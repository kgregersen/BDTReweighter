// Standard Template Library includes
#include <algorithm>

// Analysis includes
#include "Log.h"


Log::Log(const std::string & name, const LEVEL & level, std::ostream & stream) :
  m_outstream(stream),
  m_name(name),
  m_printlevel(level),
  m_currentlevel(INFO)
{

}


Log::LEVEL Log::StringToLEVEL(const std::string & str_level)
{

  // convert to lower case
  std::string lowercase(str_level);
  std::transform(lowercase.begin(), lowercase.end(), lowercase.begin(), ::tolower);

  // determine level
  if      ( lowercase == "debug"   ) return Log::DEBUG;
  else if ( lowercase == "verbose" ) return Log::VERBOSE;
  else if ( lowercase == "info"    ) return Log::INFO;
  else if ( lowercase == "warning" ) return Log::WARNING;
  else if ( lowercase == "error"   ) return Log::ERROR;
  
  // not recognised
  return Log::INDENT;
  
}

