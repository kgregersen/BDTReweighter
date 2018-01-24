#ifndef __NODE__
#define __NODE__

//stl includes
#include <vector>
#include <string>
#include <iostream>

// local includes
#include "Log.h"
#include "HistDefs.h"

// ROOT includes
#include "TH1F.h"

// forward declarations
class Branch;
class Store;
class Event;
class DecisionTree;
class TTree;


class Node {

public:


  // -----------------------
  // class to hold histogram
  // -----------------------
  class Hist {

  public:

    // constructor
    Hist(const HistDefs::Entry & histDef) : m_hist(0), m_variable(histDef.GetVariable())
    {
      m_hist = new TH1F(histDef.Name().c_str(), histDef.Name().c_str(), histDef.Nbins(), histDef.Xmin(), histDef.Xmax());
      m_hist->SetDirectory(0);
    }

    // destructor
    ~Hist()
    {
      delete m_hist;
      m_hist = 0;
    }

    // fill histogram
    void Fill(float weight) { m_hist->Fill(m_variable->Value(), weight); }

    // get name
    const std::string & Name() const { return m_variable->Name(); }
    
    // get ROOT hist
    const TH1F * ROOTHist() const { return m_hist; }

    
  private:

    // histogram and variable
    TH1F * m_hist;
    const Variable * m_variable;

  };

  // --------------------------
  // class to hold hist summary
  // --------------------------
  class Summary {

  public:

    // constructor
    Summary(Hist * initial, Hist * target, float cutValue, float chisquare, float sumInitLow, float sumTargLow, float sumInitHigh, float sumTargHigh) :
      m_initial(initial), m_target(target), m_cutValue(cutValue), m_chisquare(chisquare), m_sumInitialLow(sumInitLow), m_sumTargetLow(sumTargLow), m_sumInitialHigh(sumInitHigh), m_sumTargetHigh(sumTargHigh) {}

    ~Summary()
    {
      delete m_initial;
      delete m_target;
      m_initial = 0;
      m_target = 0;
    }
    
    // get information
    const Hist * InitialHist() const { return m_initial; }
    const Hist * TargetHist() const { return m_target; }
    float CutValue() const { return m_cutValue; }
    float Chisquare() const { return m_chisquare; }
    float SumInitialLow() const { return m_sumInitialLow; }
    float SumTargetLow() const { return m_sumTargetLow; }
    float SumInitialHigh() const { return m_sumInitialHigh; }
    float SumTargetHigh() const { return m_sumTargetHigh; }
    const std::string & Name() const { return m_initial->Name(); }

    
  private:

    // summary info
    Hist * m_initial;
    Hist * m_target;
    float m_cutValue;
    float m_chisquare;
    float m_sumInitialLow;
    float m_sumTargetLow;
    float m_sumInitialHigh;
    float m_sumTargetHigh;
    
  };
  
  
  // Node type
  enum STATUS {
    NEW,
    FIRST,
    INTERMEDIATE,
    FINAL
  };

  
  // constructors
  Node(const Store * store, Branch * input = 0);

  // destructor
  ~Node();

  // get status
  STATUS Status() const;

  // string version of status
  std::string StatusStr() const;
   
  // get input branch
  const Branch * InputBranch() const;

  // get output branch
  const Branch * OutputBranch() const;

  // grow node
  void Build(TTree * initial,
	     TTree * target,
	     const HistDefs & histDefs,
	     Branch *& b1,
	     Branch *& b2,
	     const std::vector<const DecisionTree *> & decisionTrees);

  // fill histograms
  void FillHistograms(TTree * tree, std::vector<Hist *> histSet, const std::vector<const DecisionTree *> * decisionTrees = 0) const;

  // get #initial
  float SumInitial() const;
  
  // get #target
  float SumTarget() const;

  // get weight
  float GetWeight() const;

  // set weight (and lock it)
  void SetAndLockWeight(float weight) const;

  // get number of layers in front of this node
  int NumberOfLayers() const;

  
private:

  // status
  STATUS m_status;
  
  // I/O branches
  const Branch * m_input;
  const Branch * m_output1;
  const Branch * m_output2;

  // weight
  mutable float m_weight;
  mutable bool m_weightIsSet;
  
  // sum of events
  float m_sumInitial;
  float m_sumTarget;
  
  // logger
  mutable Log m_log;

  // configuration store
  const Store * m_store;
  
};

#endif
