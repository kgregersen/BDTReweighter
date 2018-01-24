#ifndef __VARIABLES__
#define __VARIABLES__

// STL includes
#include <vector>

// forward declarations
class Variable;


class Variables {
public:

  // get variables
  static std::vector<const Variable *> & Get();

  // initialize
  static void Initialize();
  
};



#endif
