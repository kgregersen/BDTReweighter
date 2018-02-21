// local includes
#include "Node.h"
#include "Branch.h"
#include "Event.h"
#include "Store.h"
#include "DecisionTree.h"

// stl includes
#include <map>
#include <iomanip>
#include <string>
#include <limits>
#include <cmath>
#include <algorithm>

// ROOT includes
#include "TTree.h"
#include "TRandom3.h"



Node::Node(const Store * store, Method::TYPE method, Branch * input) :
  m_status(NEW),
  m_input(input),
  m_output1(0),
  m_output2(0),
  m_method(method),
  m_weight(0.),
  m_weightIsSet(false),
  m_sumInitial(-1),
  m_sumTarget(-1),
  m_log("Node"),
  m_store(store)
{

  // set this node as output node of input branch, and get sum of events from input branch
  if ( input ) {

    input->SetOutputNode(this);
    m_sumInitial = input->SumInitial();
    m_sumTarget  = input->SumTarget();

    // set to FINAL if there are fewer than twice the min number of events on the node since then it can't be split
    // (there still exist other cases where it can't be split, but we have to fill the histograms to identify these...)
    int minEvents = 0;
    m_store->getif<int>("MinEventsNode", minEvents);
    if ( m_sumTarget < 2.*minEvents || m_sumInitial < 2.*minEvents) {
      m_status = FINAL;
    }

  }

  // set log level
  std::string str_level;
  store->getif<std::string>("PrintLevel", str_level);
  if (str_level.length() > 0) {
    Log::LEVEL level = Log::StringToLEVEL(str_level);
    m_log.SetLevel(level);
  }
  
}


Node::~Node()
{

  delete m_input;
  m_input = 0;

  for (unsigned int i = 0; i < m_histSetInitial.size(); ++i) {
    delete m_histSetInitial.at(i);
    m_histSetInitial.at(i) = 0;
  }

  for (unsigned int i = 0; i < m_histSetTarget.size(); ++i) {
    delete m_histSetTarget.at(i);
    m_histSetTarget.at(i) = 0;
  }

}


void Node::Initialize(const HistDefs * histDefs)
{

  // check if this node was already initialized
  if ( m_histSetInitial.size() || m_histSetTarget.size() ) {
    m_log << Log::ERROR << "Initialize() : Histograms already initialized!" << Log::endl();
    throw(0);
  }
  
  // get variables used for splitting the tree
  const std::vector<HistDefs::Entry> & histDefEntries = histDefs->GetEntries();
  std::vector<unsigned int> indices;
  if (m_method == Method::RF || m_method == Method::ET) {
    
    // Random Forest and ExtraTrees use "feature sampling", only using random subset of the variables to grow the decision tree
    static int samplingFractionSeed = m_store->get<float>("SamplingFractionSeed");
    static float featSamplingFraction = m_store->get<float>("FeatureSamplingFraction");
    static TRandom3 ran( samplingFractionSeed );
    for (unsigned int index = 0; index < histDefEntries.size(); ++index) indices.push_back( index );
    for (unsigned int index = 0; index < histDefEntries.size(); ++index) std::swap(indices[ index ], indices[static_cast<int>(ran.Rndm()*(static_cast<float>(indices.size()) - std::numeric_limits<float>::epsilon()))] );
    indices.resize(featSamplingFraction*histDefEntries.size());
    
  }
  else {
    
    // use all variables
    for (unsigned int index = 0; index < histDefEntries.size(); ++index) {
      indices.push_back( index );
    }
    
  }
  
  // declare target and initial histograms for each variable
  for (unsigned int index : indices) {
    const HistDefs::Entry & histDef = histDefEntries.at(index);
    m_histSetInitial.push_back( new Hist(histDef) );
    m_histSetTarget.push_back( new Hist(histDef) ); 
  }
  
  
}


Node::STATUS Node::Status() const
{

  return m_status;

}
  

void Node::SetStatus(Node::STATUS status)
{

  m_status = status;

}
  

const Branch * Node::InputBranch() const
{

  return m_input;

}


const Branch * Node::OutputBranch() const
{

  if      ( m_output1 && m_output1->Pass() ) return m_output1;
  else if ( m_output2 && m_output2->Pass() ) return m_output2;
  else if ( m_output1 && m_output2 ) {
    m_log << Log::ERROR << "OutputBranch() : Output branches are not both null, but the event doesn't pass one of them!" << Log::endl();
    throw(0);
  }
  
  return nullptr;

}


const Branch * Node::OutputBranch(bool isGreater) const
{

  if (isGreater) return m_output2;
  return m_output1;
  
}


void Node::SetOutputBranch(const Branch * branch, bool isGreater)
{
  
  if (isGreater) m_output2 = branch;
  else m_output1 = branch;
  
}


void Node::Build(Branch *& b1, Branch *& b2)
{

  // get number of histograms
  int nhist = m_histSetInitial.size();
  if ( nhist == 0 ) {
    m_log << Log::ERROR << "Build() : No histograms!" << Log::endl();
    throw(0);
  }
  
  // summary for chosen node
  Summary * nodeSummary = 0;

  // vector to hold node summaries
  std::vector<Summary *> nodeSummaryVec;
  
  // determine variable to cut on and the cut value
  static int minEvents = m_store->get<int>("MinEventsNode");
  if (m_method == Method::ET) {

    // radnomly chose variable
    static int samplingFractionSeed = m_store->get<float>("SamplingFractionSeed");
    static TRandom3 ran( samplingFractionSeed );
    unsigned int ranIndex = static_cast<unsigned int>(ran.Rndm()*(static_cast<float>(nhist) - std::numeric_limits<float>::epsilon()));
    
    Hist * histTarg = m_histSetTarget.at( ranIndex );
    Hist * histInit = m_histSetInitial.at( ranIndex );
    
    const TH1F * rootHistTarg = histTarg->ROOTHist();
    const TH1F * rootHistInit = histInit->ROOTHist();

    bool keepTrying = true;
    while ( keepTrying ) {

      keepTrying = false;
      
      if ( ! m_input && ! nodeSummary ) keepTrying = true;
      
      // get integrals above and below
      unsigned int xbin = static_cast<unsigned int>(ran.Rndm()*(static_cast<float>(rootHistInit->GetNbinsX()) - std::numeric_limits<float>::epsilon()));
      Double_t sumInitLowErr  = 0;
      Double_t sumTargLowErr  = 0;
      Double_t sumInitHighErr = 0;
      Double_t sumTargHighErr = 0;
      Double_t sumInitLow  = rootHistInit->IntegralAndError(0       , xbin, sumInitLowErr );
      Double_t sumTargLow  = rootHistTarg->IntegralAndError(0       , xbin, sumTargLowErr );
      Double_t sumInitHigh = rootHistInit->IntegralAndError(xbin + 1, -1  , sumInitHighErr);
      Double_t sumTargHigh = rootHistTarg->IntegralAndError(xbin + 1, -1  , sumTargHighErr);
      
      m_log << Log::DEBUG << "xbin = " << xbin << "  sumInitLow = " << sumInitLow << "  sumTargLow = " << sumTargLow << "  sumInitHigh = " << sumInitHigh << "  sumTargHigh = " << sumTargHigh << Log::endl();
      
      // check min events on potential sub-nodes
      if (sumInitLow >= minEvents && sumInitHigh >= minEvents && sumTargLow >= minEvents && sumTargHigh >= minEvents ) {
	
	// calculate chisquare and get cut value
	float chisquare = pow(sumInitLow - sumTargLow, 2)/(pow(sumInitLowErr, 2) + pow(sumTargLowErr, 2)) + pow(sumInitHigh - sumTargHigh, 2)/(pow(sumInitHighErr, 2) + pow(sumTargHighErr, 2));
	float cutValue  = rootHistTarg->GetBinLowEdge(xbin + 1);
	
	// set node summary
	nodeSummary = new Summary(histInit, histTarg, cutValue, chisquare, sumInitLow, sumTargLow, sumInitHigh, sumTargHigh);
	
	m_log << Log::DEBUG << "Node::Summary set! xbin = " << xbin << "  sumInitLow = " << sumInitLow << "  sumTargLow = " << sumTargLow << "  sumInitHigh = " << sumInitHigh << "  sumTargHigh = " << sumTargHigh << Log::endl();

      }

    }
    
  }
  else {
    
    // calculate cut-values and chisquares
    for (int i = 0; i < nhist; ++i) {
      
      Hist * histTarg = m_histSetTarget.at(i);
      Hist * histInit = m_histSetInitial.at(i);
      
      const TH1F * rootHistTarg = histTarg->ROOTHist();
      const TH1F * rootHistInit = histInit->ROOTHist();
      
      m_log << Log::DEBUG << "targ integral = " << rootHistTarg->Integral(0,-1) << "  init integral = " << rootHistInit->Integral(0,-1) << Log::endl();
      
      float maxChisquare = 0;
      float cutValue = std::numeric_limits<float>::max();
      float sumTargetLow   = 0;
      float sumTargetHigh  = 0;
      float sumInitialLow  = 0;
      float sumInitialHigh = 0;
      
      // loop over bins in histogram
      for (int xbin = 1; xbin < rootHistTarg->GetNbinsX(); ++xbin) {
	
	// get integrals above and below
	Double_t sumInitLowErr  = 0;
	Double_t sumTargLowErr  = 0;
	Double_t sumInitHighErr = 0;
	Double_t sumTargHighErr = 0;
	Double_t sumInitLow  = rootHistInit->IntegralAndError(0       , xbin, sumInitLowErr );
	Double_t sumTargLow  = rootHistTarg->IntegralAndError(0       , xbin, sumTargLowErr );
	Double_t sumInitHigh = rootHistInit->IntegralAndError(xbin + 1, -1  , sumInitHighErr);
	Double_t sumTargHigh = rootHistTarg->IntegralAndError(xbin + 1, -1  , sumTargHighErr);
	
	m_log << Log::DEBUG << "sumInitLow = " << sumInitLow << "  sumTargLow = " << sumTargLow << "  sumInitHigh = " << sumInitHigh << "  sumTargHigh = " << sumTargHigh << Log::endl();
	
	// check min events on potential sub-nodes
	if (sumInitLow < minEvents || sumInitHigh < minEvents || sumTargLow < minEvents || sumTargHigh < minEvents ) continue;
	
	// calculate chisquare and update best candidate
	float chisquare = pow(sumInitLow - sumTargLow, 2)/(pow(sumInitLowErr, 2) + pow(sumTargLowErr, 2)) + pow(sumInitHigh - sumTargHigh, 2)/(pow(sumInitHighErr, 2) + pow(sumTargHighErr, 2));
	if (chisquare > maxChisquare) {
	  maxChisquare   = chisquare;
	  cutValue       = rootHistTarg->GetBinLowEdge(xbin + 1);
	  sumInitialLow  = sumInitLow;
	  sumInitialHigh = sumInitHigh;
	  sumTargetLow   = sumTargLow;
	  sumTargetHigh  = sumTargHigh;
	}
	
      }
      
      m_log << Log::DEBUG << "sumInitialLow = " << sumInitialLow << "  sumTargetLow = " << sumTargetLow << "  sumInitialHigh = " << sumInitialHigh << "  sumTargetHigh = " << sumTargetHigh << Log::endl();
      
      // store info for this variable
      if (maxChisquare > 0) {
	nodeSummaryVec.push_back( new Summary(histInit, histTarg, cutValue, maxChisquare, sumInitialLow, sumTargetLow, sumInitialHigh, sumTargetHigh) );
      }
      
    }

    // find highest chisquare
    for (Summary * s : nodeSummaryVec) {
      if (s->Chisquare() > nodeSummary->Chisquare()) {
	nodeSummary = s;
      }
    }

  }
  
  // sanity check
  if ( m_input == 0 && nodeSummary == 0 ) {
    m_log << Log::ERROR << "Build() : This is the first node in the tree (input branch is null), but there is no node summary - we can't build the friggin tree?!?!" << Log::endl();
    throw(0);
  }
  
  // set node status 
  if ( nodeSummary == 0 ) {
    m_status = FINAL;
  }
  else if ( m_input == 0 ) {
    m_status = FIRST;
    m_sumTarget  = nodeSummary->TargetHist() ->ROOTHist()->Integral(0, -1);
    m_sumInitial = nodeSummary->InitialHist()->ROOTHist()->Integral(0, -1);
  } 
  else {
    m_status = INTERMEDIATE;
  }
  
  // set outgoing branches
  if ( m_status == FINAL ) {
    b1 = 0;
    b2 = 0;
  }
  else {
    b1 = new Branch(m_store, this, nodeSummary->Name(), nodeSummary->CutValue(), false, nodeSummary->SumInitialLow() , nodeSummary->SumTargetLow() );
    b2 = new Branch(m_store, this, nodeSummary->Name(), nodeSummary->CutValue(), true , nodeSummary->SumInitialHigh(), nodeSummary->SumTargetHigh());
  }

  // set output branches for this node
  m_output1 = b1;
  m_output2 = b2;

  // print info
  m_log << Log::VERBOSE << "Build() : -----------> INFO <-----------" << Log::endl();
  m_log << Log::VERBOSE << "Build() : Status      : " << StatusStr() << Log::endl(); 
  m_log << Log::VERBOSE << "Build() : Sum Target  : " << m_sumTarget << Log::endl(); 
  m_log << Log::VERBOSE << "Build() : Sum Initial : " << m_sumInitial << Log::endl(); 
  m_log << Log::VERBOSE << "Build() : Cuts        : ";
  const Branch * b = InputBranch();
  while ( b ) {
    
    // get cut object
    const Branch::Cut * cut = b->CutObject();
    const Branch::Smaller * lt = dynamic_cast<const Branch::Smaller *>(cut);
    const Branch::Greater * gt = dynamic_cast<const Branch::Greater *>(cut);    
    
    // check if valid
    if ( (lt && gt) || (!lt && !gt) ) {
      m_log << Log::endl();
      m_log << Log::ERROR << "Build() : Couldn't determine if greater or smaller!" << Log::endl();
      throw(0);
    }
    
    // print cut
    m_log << cut->GetVariable()->Name() << (lt ? "<" : ">") << cut->CutValue() << "|";

    // update branch
    b = b->InputNode()->InputBranch();
    
  }
  m_log << Log::endl();
  m_log << Log::VERBOSE << "Build() : ------------------------------" << Log::endl();

  // clean up
  for (Summary * nodeSum : nodeSummaryVec) {
    delete nodeSum;
    nodeSum = 0;
  }
  for (Hist * hist : m_histSetInitial) {
    delete hist;
    hist = 0;
  }
  for (Hist * hist : m_histSetTarget) {
    delete hist;
    hist = 0;
  }
  
}


void Node::FillInitial(float weight)
{

  // fill histograms
  for (Hist * hist : m_histSetInitial) {
    hist->Fill(weight);
  }

}


void Node::FillTarget(float weight)
{

  // fill histograms
  for (Hist * hist : m_histSetTarget) {
    hist->Fill(weight);
  }

}


float Node::SumInitial() const
{

  return m_sumInitial;

}


float Node::SumTarget() const
{

  return m_sumTarget;

}


void Node::SetAndLockWeight(float weight) const
{
  
  if (m_weightIsSet == false) {
    m_weight = weight;
    m_weightIsSet = true;
  }

}

float Node::GetWeight() const
{

  return m_weight;

}

 
std::string Node::StatusStr() const
{

  if (m_status == NEW         ) return "NEW";
  if (m_status == FIRST       ) return "FIRST";
  if (m_status == INTERMEDIATE) return "INTERMEDIATE";
  if (m_status == FINAL       ) return "FINAL";
  return "NOSTATUS";
  
}
