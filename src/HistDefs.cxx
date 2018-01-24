// local includes
#include "HistDefs.h"
#include "Variable.h"
#include "Variables.h"
#include "Store.h"

// ROOT includes
#include "TTree.h"


HistDefs::HistDefs(const Store * store) :
  m_store(store),
  m_log("HistDefs")
{

  // set log's print level
  std::string str_level;
  m_store->getif<std::string>("PrintLevel", str_level);
  if (str_level.length() > 0) {
    Log::LEVEL level = Log::StringToLEVEL(str_level);
    m_log.SetLevel(level);
  }

}


void HistDefs::Initialize()
{

  // get string names for float/int variables
  const std::vector<const Variable *> & variables = Variables::Get();

  // initialize histogram definitions for float variables
  m_defs.clear();
  for (const Variable * var : variables) {
    m_defs.push_back( Entry(var) );
  }

}


void HistDefs::UpdateVariableRanges(TTree * tree)
{
  
  // prepare for loop over tree entries
  long maxEvent = tree->GetEntries();
  long reportFrac = maxEvent/(maxEvent > 100000 ? 10 : 1) + 1;
  m_log << Log::INFO << "UpdateVariableRanges() : Looping over events (" << tree->GetName() << ") : "  << maxEvent << Log::endl();
  std::clock_t start = std::clock();

  // Loop over ree entries
  for (long ievent = 0; ievent < maxEvent; ++ievent) {

    // print progress
    if( ievent > 0 && ievent % reportFrac == 0 ) {
      double duration     = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
      double frequency    = static_cast<double>(ievent) / duration;
      double timeEstimate = static_cast<double>(maxEvent - ievent) / frequency;
      m_log << Log::VERBOSE << "UpdateVariableRanges() : ---> processed : " << std::setw(4) << 100*ievent/maxEvent << "\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time : " << std::setw(4) << static_cast<int>(timeEstimate) << " sec"<< Log::endl(); 
    }
    
    // load event
    tree->GetEntry( ievent );

    // update ranges
    for (Entry & entry : m_defs) {
      float value = entry.GetVariable()->Value();
      if      (value < entry.Xmin() ) entry.SetXmin(value);
      else if (value > entry.Xmax() ) entry.SetXmax(value);
    }

  }

  // print out
  double duration  = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
  double frequency = static_cast<double>(maxEvent) / duration;
  m_log << Log::VERBOSE << "UpdateVariableRanges() : ---> processed :  100\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time :    0 sec"<< Log::endl(); 

}

const std::vector<HistDefs::Entry> & HistDefs::GetEntries() const
{

  return m_defs;

}
