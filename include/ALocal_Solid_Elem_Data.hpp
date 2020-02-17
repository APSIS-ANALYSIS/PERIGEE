#ifndef ALOCAL_SOLID_ELEM_DATA_HPP
#define ALOCAL_SOLID_ELEM_DATA_HPP
// ==================================================================
// ALocal_Solid_Elem_Data.hpp
//
// Analysis-use local solid element class with data. The data are 
// assigned to each local solid element, typcally associated with 
// quadrature points.
// 
// We may assign prestress tensor and fiber direction with each quad
// rature point in each solid element.
//
// Date: Feb. 13 2020
// Author: Ju Liu
// ==================================================================

class ALocal_Solid_Elem_Data
{
  public:
    // --------------------------------------------------------------
    // Constructor:
    // Allocate the prestress container. If the element is fluid,
    // prestress pointer is nullptr; if the element is solid, allocate
    // 7 x nqp to store the prestress for S and p. Initialize their
    // values by zero.
    // --------------------------------------------------------------
    ALocal_Solid_Elem_Data( const ALocal_Elem_wTag * const &in_aelem,
       const int &in_nqp );

    virtual ~ALocal_Solid_Elem_Data();

    virtual void print_info();

    // --------------------------------------------------------------
    // access the prestress data
    // ee ranges from 0 to nlocalelem 
    // ii ranges form 0 to nqp, the number of quadrature points
    // The output are the prestress values in the 2nd PK stress
    // and the pressure.
    // If the element ee belongs to a fluid element, this
    // will return everything by 0.0
    // --------------------------------------------------------------
    virtual void get_prestress( const int &ee,
        const int &ii, const double &S00, const double &S01,
        const double &S02, const double &S11,
        const double &S12, const double &S22,
        const double &pressure ) const;

  private:
    // number of quadrature points and the number of elemetn in the
    // subdomain (including fluid and solid element)
    const int nqp, nlocalelem;
    
    // number of solid element in this subdomain 
    int nelem_solid;

    double ** prestress;
};

#endif
