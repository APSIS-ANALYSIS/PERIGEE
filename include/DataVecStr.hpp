#ifndef DATAVECSTR_HPP
#define DATAVECSTR_HPP
// ----------------------------------------------------------------------------
// DataVecStr.hpp
// ----------------------------------------------------------------------------
// DataVecStr defines a data type for a std::vector with a string holding its
// name.
//
// Author: Ju Liu
// Date: Aug. 27 2023
// ----------------------------------------------------------------------------
#include <vector>
#include <string>

template <typename T>
class DataVecStr
{
  public:
    DataVecStr( const std::vector<T> &input_data, const std::string &input_name ) : data(input_data), name(input_name) {};
 
    DataVecStr() { data = {}; name = "undefined";}

    virtual ~DataVecStr() {};

    DataVecStr<T>& operator= ( const DataVecStr<T> &input )
    {
      if( this == &input) return *this;

      data = input.get_data();
      name = input.get_name();
    
      return *this;
    }

    std::vector<T> get_data() const {return data;}

    std::string get_name() const {return name;}

  private:
    std::vector<T> data;
    std::string name;
};

#endif
