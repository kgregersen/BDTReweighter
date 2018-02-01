#ifndef __HISTDEFS__
#define __HISTDEFS__

// stl includes
#include <string>
#include <vector>
#include <limits>
#include <iterator>

// local includes
#include "Variable.h"
#include "Log.h"


// forward declarations
class Store;
class TTree;


class HistDefs {
  
public:

  
  // class to hold single histogram definition
  class Entry {
    
  public:

    // constructor
    Entry(const Variable * variable) :
      m_variable(variable),
      m_xmin(std::numeric_limits<float>::max()),
      m_xmax(-std::numeric_limits<float>::max()),
      m_nbins(500)
    {}

    // destructor
    ~Entry() {}

    // member methods
    const std::string & Name() const { return m_variable->Name(); }
    float Xmin() const { return m_xmin; }
    float Xmax() const { return m_xmax; }
    void SetXmin(float value) { m_xmin = value; }
    void SetXmax(float value) { m_xmax = value; }
    int Nbins() const { return m_nbins; }
    const Variable * GetVariable() const { return m_variable; }

    
  private:

    // data
    const Variable * m_variable;
    float m_xmin;
    float m_xmax;
    int m_nbins;
    
  };

    
  // constructor
  HistDefs(const Store * store);

  // destructor
  ~HistDefs() {}

  // add entries to set
  void Initialize();

  // update variable ranges
  void UpdateVariableRanges(TTree * tree);
  
  // get entries
  const std::vector<Entry> & GetEntries() const;
  
  
private:

  // set of histograms entries
  std::vector<Entry> m_defs;

  // config file
  const Store * m_store;
  
  // logger
  mutable Log m_log;
  
};


#endif
