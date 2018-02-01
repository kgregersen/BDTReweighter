#ifndef __EVENT__
#define __EVENT__

// ROOT includes
#include "TString.h"

// STL includes
#include <map>
#include <typeinfo>
#include <vector>

// local includes
#include "Log.h"
#include "Store.h"

// forward declarations
class TTree;
class Store;


class Event {

public:

  // internal classes to hold information
  class VarBase {
  public:
    VarBase() {}
    virtual ~VarBase() {}
  };
  template<typename T>
  class Var : public VarBase {
  public:
    Var() : value() {}
    Var(T v) : value(v) {}
    T value;
  };

  // singleton pattern
  static Event & Instance(const Store * store = 0) {
    static Event instance(store);
    return instance;
  }
  
  // disable copy-constructor and assignment operator
  Event(const Event & other) = delete;
  void operator=(const Event & other)  = delete;

  // destructor
  ~Event() {}
  
  // connect all
  void ConnectAllVariables(TTree * tree, bool disableOtherBranches = true);

  // get variable (non-const)
  template <typename T>
  T & GetVar(const TString & name);

  // get variable (const)
  template <typename T>
  const T & GetVar(const TString & name) const;

  
private:

  // constructor
  Event(const Store * store); 

  // map of variables, store and log
  std::map<TString, VarBase *> m_varMap;
  const Store * m_store;
  mutable Log m_log;

  // store variable from tree in internal map
  template<typename T>
  void ConnectVariable(const TString & name, TTree * tree);
       
};

#include "Event.icc"

#endif
