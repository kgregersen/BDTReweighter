#ifndef __METHOD__
#define __METHOD__

// std includes
#include <string>


class Method {

public:
  
  enum TYPE {
    NONE,
    BDT,  // Boosted Decision Trees
    RF,   // Random Forest 
    ET    // Extremely Randomised Trees (ExtraTrees)
  };

  static TYPE Type(const std::string & methodStr)
  {
    if      (methodStr == "BDT") return BDT;
    else if (methodStr == "RF" ) return RF;
    else if (methodStr == "ET" ) return ET;
    return NONE;
  }
  
  static std::string String(Method::TYPE method)
  {
    if      (method == BDT) return "BDT";
    else if (method == RF ) return "RF";
    else if (method == ET ) return "ET";
    return "NONE";
  }
  
};


#endif
