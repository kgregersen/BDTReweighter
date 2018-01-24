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
    virtual TString Type() const = 0;
  };
  template<typename T>
  class Var : public VarBase {
  public:
    Var() : value() {}
    Var(T v) : value(v) {}
    TString Type() const { return typeid(value).name(); }
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
  
  void ConnectAllVariables(TTree * tree);
  void AddWeight(float w);
  TString GetType(const TString & name) const;

  template<typename T>
  void SetVar(const TString & name, const T & value);
  
  template <typename T>
  T & GetVar(const TString & name);

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
