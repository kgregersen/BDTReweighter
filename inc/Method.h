#ifndef __METHOD__
#define __METHOD__

// std includes
#include <string>


class Method {

public:
  
  enum TYPE {
    NONE,
    BDT,  // Boosted Decision Trees
    RF    // Random Forest 
  };

  static TYPE Type(const std::string & methodStr)
  {
    if      (methodStr == "BDT") return BDT;
    else if (methodStr == "RF" ) return RF;
    return NONE;
  }
  
  static std::string String(Method::TYPE method)
  {
    if      (method == BDT) return "BDT";
    else if (method == RF ) return "RF";
    return "NONE";
  }
  
};


#endif
