#ifndef __VARIABLE__
#define __VARIABLE__

// STL includes
#include <string>

// local includes
#include "Log.h"


class Variable {

public:

  // constructor
  Variable(const std::string & name) : m_name(name), m_log(name)
  {
    m_log << Log::INFO << "Added variable : " << m_name << Log::endl();
  }

  // destructor
  virtual ~Variable() {}

  // get value
  virtual float Value() const = 0;
  
  // get name
  const std::string & Name() const { return m_name; }

  
protected:

  // name of variable
  const std::string m_name;

  // logger
  Log m_log;
  
};


#endif
