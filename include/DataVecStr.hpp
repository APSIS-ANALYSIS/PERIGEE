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
struct DataVecStr
{
  std::vector<T> data {};
  std::string name {"undefined"};

  DataVecStr( const std::vector<T> &input_data, const std::string input_name ) : data(input_data), name(input_name) {};

};

#endif
