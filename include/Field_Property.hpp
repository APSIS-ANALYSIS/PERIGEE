#ifndef FIELD_PROPERTY_HPP
#define FIELD_PROPERTY_HPP
// ============================================================================
// Field_Property
// ----------------------------------------------------------------------------
// Field_Property is a data type that defines the property of a field. In
// particular, it defines:
//
// 1. Index of the field ( for multifield/multiphysics problems this is
// relevant);
// 2. The associated DOF number of the field ( scalar is 1, elasticity is 3, for
// exampel);
// 3. Whether or not the field is accociated with the geometry ( this will
// affect whether or not we have the control points associated with the field );
// 4. Name of the field ( this may help us better understand the field ).
//
// Date Created: Oct. 31 2023
// ============================================================================
#include <iostream>
#include <string>

class Field_Property
{
  public:
    Field_Property() : id(0), dofNum(1), is_geo_field(true), name("unspecified") {};

    Field_Property( const int &input_id, const int &input_dof, const bool &input_geo_flag, const std::string &input_name ) : id( input_id ), dofNum( input_dof ), is_geo_field( input_geo_flag ), name( input_name ) {};

    virtual ~Field_Property() = default;

    int get_id() const {return id;}

    int get_dofNum() const {return dofNum;}

    bool get_is_geo_field() const {return is_geo_flag;}

    std::string get_name() const {return name;}

    friend std::ostream& operator<< (std::ostream& out, const Field_Property &val)
    {
      out<<"Field name :"<<val.get_name()<<'\t';
      out<<"Field ID :"<<val.get_id()<<'\t';
      out<<"Field associated number of DOFs :"<<val.get_dofNum()<<'\t';
      
      if(val.get_is_geo_field() == true)
        out<<" and it is a geo field.\n";
      else
        out<<" and it is NOT geo field.\n";
      
      return out;
    }

  private:
    const int id {0};
    const int dofNum {1};
    const bool is_geo_field {true};
    const std::string name {"unspecified"};
};

#endif
