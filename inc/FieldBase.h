#ifndef Utilities_FieldBase_H
#define Utilities_FieldBase_H


class FieldBase {

 public:
  
  // constructor
  FieldBase() {};

  // destructor
  virtual ~FieldBase() {};
  
  // clone 
  virtual FieldBase* clone() const = 0;

};

#endif
