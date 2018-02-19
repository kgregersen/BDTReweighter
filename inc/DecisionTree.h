#ifndef __DECISIONTREE__
#define __DECISIONTREE__

// stl includes
#include <vector>
#include <utility>
#include <fstream>

// local includes
#include "Log.h"
#include "Method.h"
#include "Branch.h"

// forward declarations
class Node;
class HistDefs;
class Store;
class TTree;


class DecisionTree {

public:

  // constructor (calculate weights)
  DecisionTree(TTree * initial, TTree * target, Method::TYPE method, const Store * store, const HistDefs * histDefs);

  // constructor (apply weights)
  DecisionTree(const std::vector<std::pair<float, std::vector<const Branch::Cut *> > > tree, const Store * store);

  // destructor
  ~DecisionTree();

  // grow tree
  void GrowTree(const std::vector<const DecisionTree *> & decisionTrees);
    
  // get weight
  float GetWeight() const;

  // write to file
  void Write(std::ofstream & file) const;

  
private:
  
  // fill nodes
  void FillNodes(std::vector<Node *> layer, bool treeSwitch, const std::vector<const DecisionTree *> * decisionTrees = 0) const;

  // create new node
  void CreateNode(Branch * input, std::vector<Node *> & nextLayer);

  // add node to decision tree
  void AddNodeToTree(const Node * node);

  // helper functions
  const Node * FirstNode() const; 
  const std::vector<const Node *> FinalNodes() const;

  // TTrees for initial and target samples, and event indices for Random Forest
  TTree * m_initial;
  TTree * m_target;
  std::vector<long> * m_indicesInitial;
  std::vector<long> * m_indicesTarget;

  // method: BDT/RF
  Method::TYPE m_method;
  
  // histogram definitions
  const HistDefs * m_histDefs;
  
  // nodes
  std::vector<const Node *> m_nodes;

  // logger
  mutable Log m_log;

  // configuration store
  const Store * m_store;
  
};

#endif
