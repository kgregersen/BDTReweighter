// local includes
#include "Variables.h"
#include "Variable.h"
#include "Event.h"


// preprocessor macro to add variables
#define VARIABLE(name,type)						\
  class name : public Variable {					\
  public:								\
  name(const std::string & vname) : Variable(vname) {}			\
  float Value() const {							\
    static type & value = Event::Instance().GetVar<type>(m_name);	\
    return static_cast<float>(value);					\
  }									\
  };									\
  Variables::Get().push_back( new name(#name) );



std::vector<const Variable *> & Variables::Get()
{
  
  static std::vector<const Variable *> variables;
  return variables;

}


const Variable * Variables::Get(const std::string & name)
{

  // get set of variables
  const std::vector<const Variable *> & variables = Variables::Get();
  if ( variables.size() == 0 ) return 0;
  
  // find variable
  const Variable * variable = 0;
  for (const Variable * var : variables) {
    if ( var->Name() == name ) {
      variable = var;
    }
  }

  // return variable
  return variable;
  
}



void Variables::Initialize()
{

  if ( Variables::Get().size() == 0 ) {
    #include "VARIABLES"
  }  

}
