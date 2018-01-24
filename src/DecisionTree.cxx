// local includes
#include "DecisionTree.h"
#include "Branch.h"
#include "Node.h"
#include "Store.h"

// stl includes
#include <vector>

// ROOT includes
#include "TTree.h"



DecisionTree::DecisionTree(TTree * initial, TTree * target, const Store * store, const HistDefs & histDefs) :
  m_initial(initial),
  m_target(target),
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
  
}


DecisionTree::~DecisionTree()
{

}


void DecisionTree::GrowTree(const std::vector<const DecisionTree *> & decisionTrees)
{

  // print info
  static int counter = 1;
  m_log << Log::INFO << "GrowTree() : Decision Tree " << counter++ << Log::endl();

  // keep track of time
  std::clock_t start = std::clock();

  // declare first node
  Node * node = new Node(m_store, 0);
  m_nodes.push_back(node);

  // add nodes to tree
  while (node->Status() == Node::NEW) {

    // declare output branches
    Branch * b1 = 0;
    Branch * b2 = 0;

    // build node
    node->Build(m_initial, m_target, m_histDefs, b1, b2, decisionTrees);

    // insert sub-nodes
    if ( b1 ) m_nodes.push_back( new Node(m_store, b1) );
    if ( b2 ) m_nodes.push_back( new Node(m_store, b2) );
        
    // find next node to build
    for (Node * n : m_nodes) {
      if ( n->Status() == Node::NEW ) {
	node = n;
	break;
      }
    }
    
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
    //float w = target/initial; 
    weights.at(i)       = w;
    sumOfWeights       += w*initial;
    sumTarget          += target;
  }
  for (unsigned int i = 0; i < weights.size(); ++i) {
    finalNodes.at(i)->SetAndLockWeight( sumTarget*weights[i]/sumOfWeights );
    //finalNodes.at(i)->SetAndLockWeight( weights[i] );
  }
  
  // time spent on growing tree
  double duration = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
  
  // print info
  m_log << Log::INFO << "GrowTree() : ----------------> INFO <----------------" << Log::endl();
  m_log << Log::INFO << "GrowTree() : Time spent  : " << duration << " sec" << Log::endl();
  m_log << Log::INFO << "GrowTree() : Final nodes : " << finalNodes.size() << " (out of " << m_nodes.size() << ")" << Log::endl();
  for (unsigned int i = 0; i < finalNodes.size(); ++i) {

    const Node * node = finalNodes[i];

    m_log << Log::INFO << "GrowTree() : ---> Weight = " << std::setw(10) << std::left << weights.at(i) << " Cuts = ";
    const Branch * b = node->InputBranch();
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
  for (Node * node : m_nodes) {
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
