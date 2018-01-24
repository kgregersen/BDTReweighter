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

// ROOT includes
#include "TTree.h"



Node::Node(const Store * store, Branch * input) :
  m_status(NEW),
  m_input(input),
  m_output1(0),
  m_output2(0),
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

  m_log << Log::INFO << "~Node() : Called" << Log::endl();
  
  delete m_input;
  m_input = 0;
  
}

Node::STATUS Node::Status() const
{

  return m_status;

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


void Node::Build(TTree * initial, TTree * target, const HistDefs & histDefs, Branch *& b1, Branch *& b2, const std::vector<const DecisionTree *> & decisionTrees)
{

  // check if this node can be build
  static int maxLayers = m_store->get<int>("MaxTreeLayers");
  if ( NumberOfLayers() >= maxLayers ) {
    m_log << Log::VERBOSE << "Build() : Too many layers before node to build it - finalizing it instead" << Log::endl();
    m_status = FINAL;
    return;
  }
  
  // declare target and initial histograms
  std::vector<Hist *> histSetInitial;
  std::vector<Hist *> histSetTarget;
  const std::vector<HistDefs::Entry> & histDefEntries = histDefs.GetEntries();
  for (const HistDefs::Entry & histDef : histDefEntries) {
    histSetInitial.push_back( new Hist(histDef) );
    histSetTarget.push_back( new Hist(histDef) );
  }
  
  // fill histograms
  FillHistograms(target , histSetTarget);
  FillHistograms(initial, histSetInitial, &decisionTrees);
  
  // calculate cut-values and chisquares
  std::vector<Summary *> nodeSummaryVec;
  static int minEvents = m_store->get<int>("MinEventsNode");
  for (unsigned int i = 0; i < histDefEntries.size(); ++i) {

    Hist * histTarg = histSetTarget.at(i);
    Hist * histInit = histSetInitial.at(i);

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
  Summary * nodeSummary = nodeSummaryVec.size() > 0 ? nodeSummaryVec[0] : 0;
  for (Summary * s : nodeSummaryVec) {
    if (s->Chisquare() > nodeSummary->Chisquare()) {
      nodeSummary = s;
    }
  }

  // sanity check
  if ( m_input == 0 && nodeSummary == 0 ) {
    m_log << Log::ERROR << "Build() : This is the first node in the tree (input branch is null), but there is no node summary...?" << Log::endl();
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
  
}


void Node::FillHistograms(TTree * tree, std::vector<Hist *> histSet, const std::vector<const DecisionTree *> * decisionTrees) const
{

  // prepare for loop over tree entries
  long maxEvent = tree->GetEntries();
  long reportFrac = maxEvent/(maxEvent > 100000 ? 10 : 1) + 1;
  m_log << Log::VERBOSE << "FillHistograms() : Looping over events (" << tree->GetName() << ") : "  << maxEvent << Log::endl();
  std::clock_t start = std::clock();

  // Loop over ree entries
  for (long ievent = 0; ievent < maxEvent; ++ievent) {

    // print progress
    if( ievent > 0 && ievent % reportFrac == 0 ) {
      double duration     = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
      double frequency    = static_cast<double>(ievent) / duration;
      double timeEstimate = static_cast<double>(maxEvent - ievent) / frequency;
      m_log << Log::VERBOSE << "FillHistograms() : ---> processed : " << std::setw(4) << 100*ievent/maxEvent << "\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time : " << std::setw(4) << static_cast<int>(timeEstimate) << " sec"<< Log::endl(); 
    }
    
    // get event
    tree->GetEntry( ievent );

    // apply cuts
    bool pass = true;
    const Branch * b = InputBranch();
    while ( b ) {
      if ( ! b->Pass() ) {
	pass = false;
	break;
      }
      b = b->InputNode()->InputBranch();
    }
    if ( ! pass ) continue;

    // get event weight
    static const std::string & eventWeightName = m_store->get<std::string>("EventWeightVariableName");
    static float & eventWeight = Event::Instance().GetVar<float>(eventWeightName);
    
    // get weights from previous trees
    float dtreeWeight = 1.;
    if ( decisionTrees ) {
      for (const DecisionTree * dtree : *decisionTrees) {
	dtreeWeight *= dtree->GetWeight();
      }
    }
    
    // fill histograms
    for (Hist * hist : histSet) {
      hist->Fill(eventWeight*dtreeWeight);
    }
    
  }

  // print out
  double duration  = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
  double frequency = static_cast<double>(maxEvent) / duration;
  m_log << Log::VERBOSE<< "FillHistograms() : ---> processed :  100\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time :    0 sec"<< Log::endl(); 

  
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


int Node::NumberOfLayers() const
{

  // initialize return value
  int layers = 0;

  // propagate up the tree starting from the input node and count layers before reaching the first node
  const Branch * b = InputBranch();
  while ( b ) {
    ++layers;
    b = b->InputNode()->InputBranch();
  }

  // return nmber of layers
  return layers;

}

 
std::string Node::StatusStr() const
{

  if (m_status == NEW         ) return "NEW";
  if (m_status == FIRST       ) return "FIRST";
  if (m_status == INTERMEDIATE) return "INTERMEDIATE";
  if (m_status == FINAL       ) return "FINAL";
  return "NOSTATUS";
  
}
