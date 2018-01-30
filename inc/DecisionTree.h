#ifndef __DECISIONTREE__
#define __DECISIONTREE__

// stl includes
#include <vector>
#include <utility>
#include <fstream>

// local includes
#include "Log.h"
#include "HistDefs.h"

// forward declarations
class Node;
class Branch;
class Store;
class TTree;


class DecisionTree {

public:

  // constructor
  DecisionTree(TTree * initial, TTree * target, const Store * store, const HistDefs & histDefs);

  // destructor
  ~DecisionTree();

  // grow tree
  void GrowTree(const std::vector<const DecisionTree *> & decisionTrees);

  // fill nodes
  void FillNodes(std::vector<Node *> buildNodes, const std::vector<const DecisionTree *> * decisionTrees = 0) const;

  // create new node
  void CreateNode(Branch * input, std::vector<Node *> & nextLayer);

  // add node to decision tree
  void AddNodeToTree(const Node * node);
    
  // get weight
  float GetWeight() const;

  // helper functions
  const Node * FirstNode() const; 
  const std::vector<const Node *> FinalNodes() const;

  // write to file
  void Write(std::ofstream & file) const;

  
private:

  // TTrees for initial and target samples
  TTree * m_initial;
  TTree * m_target;

  // histogram definitions
  const HistDefs & m_histDefs;
  
  // nodes
  std::vector<const Node *> m_nodes;

  // logger
  mutable Log m_log;

  // configuration store
  const Store * m_store;
  
};

#endif
