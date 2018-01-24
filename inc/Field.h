#ifndef Utilities_Field_H
#define Utilities_Field_H

// Analysis includes
#include "FieldBase.h"


template <class T> 
class Field : public FieldBase {
  
public:

  // constructor
  Field(const T& value) : m_value(value) {}
  
  // copy constructor
  template <class U>
  Field(const Field<U>& other) : m_value( U(other.value) ) { }

  // assignment operator
  template <class U>
  Field& operator=(const Field<U>& other) 
  { 
    m_value = U(other.value) ; 
    return *this;
  }
  
  // destructor
  virtual ~Field() {}
  
  // method to retrieve value of field
  const T& get() const { return m_value ; }
  
  // clone
  FieldBase* clone() const { return new Field<T>( m_value ) ; }


private:

  T m_value;

};

#endif
