// local includes
#include "DecisionTree.h"
#include "Branch.h"
#include "Node.h"
#include "Store.h"
#include "Event.h"
#include "HistDefs.h"

// stl includes
#include <vector>
#include <algorithm>
#include <limits>

// ROOT includes
#include "TTree.h"
#include "TRandom3.h"



DecisionTree::DecisionTree(TTree * initial, TTree * target, Method::TYPE method, const Store * store, const HistDefs * histDefs) :
  m_initial(initial),
  m_target(target),
  m_indicesInitial(0),
  m_indicesTarget(0),
  m_method(method),
  m_histDefs(histDefs),
  m_log("DecisionTree"),
  m_store(store)
{
  
  std::string str_level;
  m_store->getif<std::string>("PrintLevel", str_level);
  if (str_level.length() > 0) {
    Log::LEVEL level = Log::StringToLEVEL(str_level);
    m_log.SetLevel(level);
  }

  m_log << Log::INFO << "DecisionTree() : Preparing sorted vector of event indices" << Log::endl();

  long maxEventInit = m_initial->GetEntries();
  long maxEventTarg = m_target->GetEntries();
  static float samplingFraction   = m_store->get<float>("SamplingFraction");
  static int samplingFractionSeed = m_store->get<float>("SamplingFractionSeed");
  static TRandom3 ran( samplingFractionSeed );
  static bool bagging = false;
  m_store->getif<bool>("Bagging", bagging); 
  
  if (m_method == Method::BDT && ! bagging ) {

    std::vector<long> indices;
    indices.reserve(maxEventInit);
    // get unique vector of indices to subset of events
    // ---> initial
    for (long ievent = 0; ievent < maxEventInit; ++ievent) indices.push_back(ievent);
    for (long ievent = 0; ievent < maxEventInit; ++ievent) std::swap(indices[ievent], indices[static_cast<int>(ran.Rndm()*(static_cast<float>(indices.size()) - std::numeric_limits<float>::epsilon()))] );
    m_indicesInitial = new std::vector<long>(indices.begin(), indices.begin() + samplingFraction*maxEventInit);
    // ---> target
    indices.clear();
    indices.reserve(maxEventTarg);
    for (long ievent = 0; ievent < maxEventTarg; ++ievent) indices.push_back(ievent);
    for (long ievent = 0; ievent < maxEventTarg; ++ievent) std::swap(indices[ievent], indices[static_cast<int>(ran.Rndm()*(static_cast<float>(indices.size()) - std::numeric_limits<float>::epsilon()))] );
    m_indicesTarget = new std::vector<long>(indices.begin(), indices.begin() + samplingFraction*maxEventTarg);

  }
  else {
    
    // use same sample of data for all trees in BDT when bagging
    if (m_method == Method::BDT && bagging) {
      m_log << Log::INFO << "DecisionTree() : Method = BDT, 'Bagging' activated" << Log::endl();
      ran.SetSeed( samplingFractionSeed );  
    }
    // get non-unique vector of indices to subset of events (Random Forest, Extremely Randomised Trees)
    // ---> initial
    m_indicesInitial = new std::vector<long>();
    m_indicesInitial->reserve(samplingFraction*maxEventInit);
    while ( m_indicesInitial->size() < samplingFraction*maxEventInit ) m_indicesInitial->push_back( static_cast<long>(ran.Rndm()*maxEventInit) );
    // ---> target
    m_indicesTarget = new std::vector<long>();
    m_indicesTarget->reserve(samplingFraction*maxEventInit);
    while ( m_indicesTarget->size() < samplingFraction*maxEventTarg ) m_indicesTarget->push_back( static_cast<long>(ran.Rndm()*maxEventTarg) );

  }

  // need to sort to optimise reading of TTree (TTree::GetEntry(index) reads in chunks of sequential data, so we don't want to jump around in indices...)
  std::sort(m_indicesInitial->begin(), m_indicesInitial->end());
  std::sort(m_indicesTarget->begin(), m_indicesTarget->end());

  m_log << Log::INFO << "DecisionTree() : Sorted vector of indices created!" << Log::endl();
  
}


DecisionTree::DecisionTree(const std::vector<std::pair<float, std::vector<const Branch::Cut *> > > tree, const Store * store) :
  m_initial(0),
  m_target(0),
  m_indicesInitial(0),
  m_indicesTarget(0),
  m_method(Method::NONE),
  m_histDefs(0),
  m_log("DecisionTree"),
  m_store(store)
{

  std::string str_level;
  m_store->getif<std::string>("PrintLevel", str_level);
  if (str_level.length() > 0) {
    Log::LEVEL level = Log::StringToLEVEL(str_level);
    m_log.SetLevel(level);
  }
  
  // declare first node
  Node * firstNode = new Node(m_store, m_method, 0); 
  firstNode->SetStatus(Node::FIRST);   
  AddNodeToTree(firstNode);
  
  // loop over entries (weight, cuts) and add nodes and branches
  for (unsigned int inode = 0; inode < tree.size(); ++inode) {

    // get final node weight
    float weight = tree.at(inode).first;

    // get cuts leading to the final node
    const std::vector<const Branch::Cut *> & cuts = tree.at(inode).second;

    // now add nodes and branches as needed to reconstruct this part of the decision tree
    Node * node = firstNode;
    for (const Branch::Cut * cut : cuts) {

      // try to fetch output branch corresponding to the current cut
      bool isGreater = dynamic_cast<const Branch::Greater *>(cut);
      Branch * branch = const_cast<Branch *>(node->OutputBranch(isGreater));
      
      // if branch doesn't exist, then create it (and its output node), otherwise fetch its output node
      if ( ! branch ) {
	Branch * b = new Branch(m_store, node, cut->GetVariable()->Name(), cut->CutValue(), isGreater, 0, 0);
	node->SetOutputBranch(b, isGreater);
	node = new Node(m_store, m_method, b);
	node->SetStatus(Node::INTERMEDIATE);
	AddNodeToTree(node);
      }
      else {
	node = const_cast<Node *>(branch->OutputNode());
      }

    }

    // set weight and status of final node
    node->SetAndLockWeight(weight);
    node->SetStatus(Node::FINAL);
    
  }

  // print tree to screen
  const std::vector<const Node *> finalNodes = FinalNodes();
  m_log << Log::VERBOSE << "DecisionTree() : ----------------> VERBOSE <----------------" << Log::endl();
  m_log << Log::VERBOSE << "DecisionTree() : Final nodes : " << finalNodes.size() << " (out of " << m_nodes.size() << ")" << Log::endl();
  for (unsigned int i = 0; i < finalNodes.size(); ++i) {   
    m_log << Log::VERBOSE << "DecisionTree() : ---> Weight = " << std::setw(10) << std::left << finalNodes.at(i)->GetWeight() << " Cuts = ";
    const Branch * b = finalNodes.at(i)->InputBranch();
    while ( b ) {     
      // get cut object
      const Branch::Cut * cut = b->CutObject();
      const Branch::Smaller * lt = dynamic_cast<const Branch::Smaller *>(cut);
      const Branch::Greater * gt = dynamic_cast<const Branch::Greater *>(cut);    
      // check if valid
      if ( (lt && gt) || (!lt && !gt) ) {
  	m_log << Log::endl();
  	m_log << Log::ERROR << "DecisionTree() : Couldn't determine if cut is greater or smaller!" << Log::endl();
  	throw(0);
      }
      // print cut
      m_log << cut->GetVariable()->Name() << (lt ? "<" : ">") << cut->CutValue() << "|";
      // update branch
      b = b->InputNode()->InputBranch();
    }
    m_log << Log::endl();
  }
  m_log << Log::VERBOSE << "DecisionTree() : ----------------------------------------" << Log::endl();

  
}


DecisionTree::~DecisionTree()
{

  for (unsigned int i = 0; i < m_nodes.size(); ++i) {
    delete m_nodes.at(i);
    m_nodes.at(i) = 0;
  }
  
}


void DecisionTree::GrowTree(const std::vector<const DecisionTree *> & decisionTrees)
{

  // print info
  static int counter = 1;
  m_log << Log::INFO << "GrowTree() : Decision Tree " << counter++ << ", method = " << Method::String(m_method) << Log::endl();

  // keep track of time
  std::clock_t start = std::clock();

  // declare first node
  Node * node = new Node(m_store, m_method, 0);

  // declare vector to hold nodes to be build
  std::vector<Node *> layer;
  layer.push_back(node);
  
  // add nodes to tree
  int nlayers = 0;
  while ( layer.size() > 0 ) {

    // check number of layers
    static int maxLayers = m_store->get<int>("MaxTreeLayers");
    if ( nlayers >= maxLayers ) {

      // print verbose message
      m_log << Log::VERBOSE << "GrowTree() : Max layers reached - finalizing nodes!" << Log::endl();

      // set status of nodes to FINAL and add to decision tree
      for (Node * node : layer) {
	node->SetStatus( Node::FINAL );
	AddNodeToTree( node );
      }

      // break out of loop (no more layers to grow)
      break;
      
    }

    // initialise histograms on nodes
    for (Node * node : layer) {
      node->Initialize(m_histDefs);
    }
      
    // fill nodes (first target, then initial)
    FillNodes(layer, 0);
    FillNodes(layer, 1, &decisionTrees); 
    
    // prepare vector for next layer of nodes
    std::vector<Node *> nextLayer;
        
    // build nodes
    for (Node * node : layer) {
      
      // declare output branches
      Branch * b1 = 0;
      Branch * b2 = 0;
      
      // build node
      node->Build(b1, b2);

      // add to decision tree nodes
      AddNodeToTree(node);
      
      // create sub-nodes 
      if ( b1 ) CreateNode(b1, nextLayer);
      if ( b2 ) CreateNode(b2, nextLayer);
      
    }

    // set next layer
    layer = nextLayer;
    
    // increment layer counter
    ++nlayers;
    
  }

  // calculate and set weights on final nodes
  std::vector<const Node *> finalNodes = FinalNodes();
  std::vector<double> weights(finalNodes.size());
  double sumOfWeights = 0;
  double sumTarget    = 0;
  for (unsigned int i = 0; i < finalNodes.size(); ++i) {
    const Node * node = finalNodes[i];
    double target  = node->SumTarget();
    double initial = node->SumInitial();
    if ( ! (initial > 0) ) {
      m_log << Log::ERROR << "GrowTree() : initial is not positive! (initial = " << initial << ")" << Log::endl();
      throw(0);
    }
    static float learningRate = m_store->get<float>("LearningRate");
    double w = exp(learningRate*log(target/initial));
    weights.at(i)  = w;
    sumOfWeights  += w*initial;
    sumTarget     += target;
  }
  for (unsigned int i = 0; i < weights.size(); ++i) {
    finalNodes.at(i)->SetAndLockWeight( sumTarget*weights[i]/sumOfWeights );
  }
  
  // time spent on growing tree
  double duration = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
  
  // print tree to screen
  m_log << Log::INFO << "GrowTree() : ----------------> INFO <----------------" << Log::endl();
  m_log << Log::INFO << "GrowTree() : Time spent  : " << duration << " sec" << Log::endl();
  m_log << Log::INFO << "GrowTree() : Final nodes : " << finalNodes.size() << " (out of " << m_nodes.size() << ")" << Log::endl();
  for (unsigned int i = 0; i < finalNodes.size(); ++i) {   
    m_log << Log::INFO << "GrowTree() : ---> Weight = " << std::setw(10) << std::left << weights.at(i) << " Cuts = ";
    const Branch * b = finalNodes.at(i)->InputBranch();
    while ( b ) {     
      // get cut object
      const Branch::Cut * cut = b->CutObject();
      const Branch::Smaller * lt = dynamic_cast<const Branch::Smaller *>(cut);
      const Branch::Greater * gt = dynamic_cast<const Branch::Greater *>(cut);    
      // check if valid
      if ( (lt && gt) || (!lt && !gt) ) {
	m_log << Log::endl();
	m_log << Log::ERROR << "GrowTree() : Couldn't determine if cut is greater or smaller!" << Log::endl();
	throw(0);
      }
      // print cut
      m_log << cut->GetVariable()->Name() << (lt ? "<" : ">") << cut->CutValue() << "|";
      // update branch
      b = b->InputNode()->InputBranch();
    }
    m_log << Log::endl();
  }
  m_log << Log::INFO << "GrowTree() : ----------------------------------------" << Log::endl();

  // clean up
  delete m_indicesInitial;
  delete m_indicesTarget;
  m_indicesInitial = 0;
  m_indicesTarget = 0;
  
}


void DecisionTree::CreateNode(Branch * input, std::vector<Node *> & nextLayer) 
{

  // declare node
  Node * node = new Node(m_store, m_method, input);

  // check if it's a FINAL node or if we can grow it further
  if (node->Status() == Node::FINAL) {
    AddNodeToTree( node );
  }
  else {
    nextLayer.push_back( node );
  }

}


void DecisionTree::FillNodes(std::vector<Node *> layer, bool treeSwitch, const std::vector<const DecisionTree *> * decisionTrees) const
{

  // get TTree
  TTree * tree = 0;
  const std::vector<long> * indices = 0;
  if (treeSwitch) {
    tree = m_initial;
    indices = m_indicesInitial;
  }
  else {
    tree = m_target;
    indices = m_indicesTarget;
  }
  
  // Loop over TTree entries
  std::clock_t start = std::clock();
  long maxEvent = indices->size();
  long reportFrac = maxEvent/(maxEvent > 100000 ? 10 : 1) + 1;
  m_log << Log::VERBOSE << "FillNodes() : Looping over events (" << tree->GetName() << ") : "  << maxEvent << Log::endl();
  for (long ievent = 0; ievent < maxEvent; ++ievent) {

    // print progress
    if( ievent > 0 && ievent % reportFrac == 0 ) {
      double duration     = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
      double frequency    = static_cast<double>(ievent) / duration;
      double timeEstimate = static_cast<double>(maxEvent - ievent) / frequency;
      m_log << Log::VERBOSE << "FillNodes() : ---> processed : " << std::setw(4) << 100*ievent/maxEvent << "\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time : " << std::setw(4) << static_cast<int>(timeEstimate) << " sec"<< Log::endl(); 
    }
    
    // get event
    tree->GetEntry( indices->at(ievent) );

    // loop over nodes
    for (Node * node : layer) {
    
      // apply cuts
      bool pass = true;
      const Branch * b = node->InputBranch();
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
      static const float & eventWeight = Event::Instance().GetVar<float>(eventWeightName);
      
      // fill nodes
      if ( treeSwitch ) {

	if ( m_method == Method::BDT ) {
	  if ( ! decisionTrees ) {
	    m_log << Log::ERROR << "FillNodes() : Previous decision trees not provided for BDT!" << Log::endl();
	    throw(0);
	  }
	  float dtreeWeight = 1.;
	  for (const DecisionTree * dtree : *decisionTrees) {
	    dtreeWeight *= dtree->GetWeight();
	  }
	  node->FillInitial(eventWeight*dtreeWeight);	  
	}
	else {
	  node->FillInitial(eventWeight);
	}
	
      }
      else {
	node->FillTarget(eventWeight);
      }
       
      // all nodes are orthogonal, so we can break the loop here
      break;

    }
    
  }

  // print out
  double duration  = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
  double frequency = static_cast<double>(maxEvent) / duration;
  m_log << Log::VERBOSE<< "FillNodes() : ---> processed :  100\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time :    0 sec"<< Log::endl(); 
  
}


float DecisionTree::GetWeight() const
{
  
  // to identify which node the current event belongs to, get the first node in the tree
  const Node * node = FirstNode();

  // and then propagate down the tree until reaching a final node
  while(node->Status() != Node::FINAL) {
    node = node->OutputBranch()->OutputNode(); 
  }

  // return the weight 
  return node->GetWeight();

}


void DecisionTree::AddNodeToTree(const Node * node)
{

  // check if this node is already in the tree
  for (const Node * n : m_nodes) {
    if (node == n) {
      m_log << Log::ERROR << "AddNodeToTree() : This node already exists in the decision tree!" << Log::endl();
      throw(0);
    }
  }

  // add node to decision tree
  m_nodes.push_back( node );
  
}


const Node * DecisionTree::FirstNode() const
{

  // check if there are any nodes
  if ( m_nodes.size() == 0 ) {
    m_log << Log::ERROR << "FirstNode() : Number of nodes = " << m_nodes.size() << Log::endl();
    throw(0);
  }

  // the first node should be the first entry
  const Node * node = *(m_nodes.begin());

  // but to make sure, we check
  const Branch * b = node->InputBranch();
  while ( b ) {
    b = b->InputNode()->InputBranch();
  }

  // check status
  if ( node->Status() != Node::FIRST ) {
    m_log << Log::ERROR << "FirstNode() : Couldn't find first node! (Number of nodes = " << m_nodes.size() << ")" << Log::endl();
    throw(0);
  }

  // return first node in tree
  return node;

}


const std::vector<const Node *> DecisionTree::FinalNodes() const
{

  // initialize set of nodes
  std::vector<const Node *> finalNodes;

  // loop through nodes and identify final nodes
  for (const Node * node : m_nodes) {
    if ( node->Status() == Node::FINAL ) {
      finalNodes.push_back(node);
    }
  }

  // check if we found any
  if ( finalNodes.size() == 0 ) {
    m_log << Log::ERROR << "FinalNodes() : Couldn't find any final nodes!" << Log::endl();
    throw(0);
  }

  // return nodes
  return finalNodes;

}


void DecisionTree::Write(std::ofstream & file) const
{

  // instantiate tree counter
  static int counter = 1;

  // initial print
  file << "# Decision Tree : " << counter << "\n"; 
  ++counter;

  // print weights and corresponding cuts
  const std::vector<const Node *> & finalNodes = FinalNodes();
  for (const Node * node : finalNodes) {

    // print weight
    file << "weight=" << node->GetWeight() << ":";

    // print cuts for each final node
    const Branch * b = node->InputBranch();
    while ( b ) {

      // get cut object
      const Branch::Cut * cut = b->CutObject();
      const Branch::Smaller * lt = dynamic_cast<const Branch::Smaller *>(cut);
      const Branch::Greater * gt = dynamic_cast<const Branch::Greater *>(cut);    
      
      // check if valid
      if ( (lt && gt) || (!lt && !gt) ) {
	m_log << Log::ERROR << "Write() : Couldn't determine if greater or smaller!" << Log::endl();
	throw(0);
      }

      // print cut
      file << cut->GetVariable()->Name() << (lt ? "<" : ">") << cut->CutValue() << "|";
	
      // update branch
      b = b->InputNode()->InputBranch();

    }

    // end line for this final node
    file << "\n";
    
  }

  // add an extra empty line before next tree is printed
  file << "\n";

  
}

