#ifndef __LOG__
#define __LOG__


// Standard Library includes
#include <iostream>
#include <iomanip>
#include <string>

// Color definitions for terminal output
#define RESET       "\033[0m"
#define BLACK       "\033[30m"             /* Black */
#define RED         "\033[31m"             /* Red */
#define GREEN       "\033[32m"             /* Green */
#define YELLOW      "\033[33m"             /* Yellow */
#define BLUE        "\033[34m"             /* Blue */
#define MAGENTA     "\033[35m"             /* Magenta */
#define CYAN        "\033[36m"             /* Cyan */
#define WHITE       "\033[37m"             /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


class Log {

public: 
  
  enum LEVEL {
    DEBUG   = 0,
    VERBOSE = 1,
    INFO    = 2,
    WARNING = 3,
    ERROR   = 4,
    FATAL   = 5,
    INDENT  = 6
  };

  struct endl { /*empty structure for template tagging*/ };

  // Constructor: User provides custom output stream, or uses default (std::cout)
  Log (const std::string & name, const LEVEL & level = INFO, std::ostream & stream = std::cout);
  
  // Templated ostream operator
  template<typename T> 
  Log& operator<<(const T & data);
 
  // Get/Set print level
  const LEVEL & GetLevel() const     { return m_printlevel;  }
  void SetLevel(const LEVEL & level) { m_printlevel = level; }
  
  // Get/Set name
  const std::string & GetName() const    { return m_name; }
  void SetName(const std::string & name) { m_name = name; }

  // convert std::string to Log::LEVEL
  static LEVEL StringToLEVEL(const std::string & str_level);


private:
  
  std::ostream & m_outstream;
  std::string    m_name;
  LEVEL          m_printlevel;
  LEVEL          m_currentlevel;

};

#include "Log.icc"

#endif
