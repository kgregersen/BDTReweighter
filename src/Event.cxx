// local includes
#include "Event.h"
#include "Store.h"

// STL includes
#include <algorithm>
#include <vector>
#include <string>

// preprocessor macro inserting code which connects TTree and Event
#define VARIABLE(name, type) ConnectVariable<type>(#name, tree);



Event::Event(const Store * store) :
  m_store(store),
  m_log("Event")
{

  std::string str_level;
  store->getif<std::string>("PrintLevel", str_level);
  if (str_level.length() > 0) {
    Log::LEVEL level = Log::StringToLEVEL(str_level);
    m_log.SetLevel(level);
  }
  
}


void Event::ConnectAllVariables(TTree * tree)
{

  // disable all branches
  tree->SetBranchStatus("*",0);

  // connect reweighting variables
  #include "VARIABLES"

  // connect event weight
  const std::string & eventWeightName = m_store->get<std::string>("EventWeightVariableName");
  ConnectVariable<float>(eventWeightName.c_str(), tree);

}


void Event::AddWeight(float w)
{
  
  static const std::string & eventWeightVariable = m_store->get<std::string>("EventWeightVariable");
  static float & weight = GetVar<float>(eventWeightVariable.c_str());
  weight *= w;
  
}


TString Event::GetType(const TString & name) const
{

  std::map<TString, VarBase *>::const_iterator it = m_varMap.find( name );
  if (it == m_varMap.end()) {
    m_log << Log::ERROR << "GetType() : " << name.Data() << " is not in map!" << Log::endl();
    throw(0);
  }

  const VarBase * var = it->second;
  
  return var->Type();  

}
