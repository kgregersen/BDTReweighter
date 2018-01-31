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


void Event::ConnectAllVariables(TTree * tree, bool disableOtherBranches)
{

  // disable all branches
  if (disableOtherBranches) tree->SetBranchStatus("*",0);
  else tree->SetBranchStatus("*",1);
  
  // connect reweighting variables
  #include "VARIABLES"

  // connect event weight
  const std::string & eventWeightName = m_store->get<std::string>("EventWeightVariableName");
  ConnectVariable<float>(eventWeightName.c_str(), tree);

}


