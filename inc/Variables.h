#ifndef __VARIABLES__
#define __VARIABLES__

// STL includes
#include <vector>
#include <string>

// forward declarations
class Variable;


class Variables {
public:

  // get variables
  static std::vector<const Variable *> & Get();
  static const Variable * Get(const std::string & name);  

  // initialize
  static void Initialize();
  
};



#endif
