#ifndef DATAVECSTR_HPP
#define DATAVECSTR_HPP
// ----------------------------------------------------------------------------
// DataVecStr.hpp
// ----------------------------------------------------------------------------
// DataVecStr defines a data type for a std::vector with a string holding its
// name and an AssociateObject flag identifying if it is a node-related or
// cell-related data.
//
// This class is designed for the input and output of data between the code and VTK
// file, primarily.
//
// Author: Ju Liu
// Date: Aug. 27 2023
// ----------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>

enum class AssociateObject
{
  Node,
  Cell,
};

template <typename T>
class DataVecStr
{
  public:
    DataVecStr( const std::vector<T> &input_data, const std::string &input_name, const AssociateObject &input_type ) : data(input_data), name(input_name), object(input_type) {};
 
    DataVecStr() { data = {}; name = "undefined"; object = AssociateObject::Node; }

    virtual ~DataVecStr() = default;

    DataVecStr<T>& operator= ( const DataVecStr<T> &input )
    {
      if( this == &input ) return *this;

      data   = input.get_data();
      name   = input.get_name();
      object = input.get_object();
    
      return *this;
    }

    std::vector<T> get_data() const {return data;}

    std::string get_name() const {return name;}

    AssociateObject get_object() const {return object;}

    int get_data_size() const {return static_cast<int>(data.size());}

    friend std::ostream& operator<< (std::ostream& out, const DataVecStr<T>& val)
    {
      out<<val.name<<'\t';

      if(val.object == AssociateObject::Node) out<<"Node \t";
      else if(val.object == AssociateObject::Cell) out<<"Cell \t";
      else out<<"unknown object \t";
      
      for( auto it = val.data.begin(); it != val.data.end(); ++it ) out<<*it<<'\t';
      
      return out;
    }

  private:
    // ------------------------------------------------------------------------
    // data is the container that holds the values
    // name is the name to be writen or read from the file
    // object = AssociateObject::Node or AssociateObject::Cell.
    // ------------------------------------------------------------------------
    std::vector<T> data;
    std::string name;
    AssociateObject object;
};

#endif
