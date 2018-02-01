#ifndef __STORE__
#define __STORE__

// STL includes
#include <map>
#include <string>
#include <vector>
#include <fstream>

// Framework includes
#include "Log.h"
#include "Store.h"

// Forward declarations
class FieldBase;


class Store {

 public:

  // constructor
  Store();

  // copy constructor
  Store(const Store & other);

  // assignment operator
  Store & operator=(const Store & other);

  // destructor
  virtual ~Store();
  
  // get data field from store
  template <class T> 
  const T & get(const std::string & key) const;

  // get data field from store (if it is there)
  template <class T> 
  void getif(const std::string & key, T & value) const;

  // put data field in store
  template <class T> 
  void put(const std::string & key, const T & value, bool overwrite = false);

  // remove data field in store
  void remove(const std::string & key);
  
  // flush store (remove all entries)
  void flush();

  // write store to file
  void write(std::ofstream & file) const;
  
  // create store from card file
  static Store * createStore(const char * filename);
  static Store * createStore(const std::string & filename) { return createStore( filename.c_str() ) ; }


 private:

  // map of data fields
  std::map<std::string, FieldBase *> m_data;
  
  // convert string field to other type, used by createStore()
  template <class T>
  static T convertField(const std::string & value);

  // convert vector of strings field to vector of other type, used by createStore()
  template <class T>
  static std::vector<T> convertVector(const std::vector<std::string> & value);

  // Log
  static Log m_log;

};

#include "Store.icc"

#endif
