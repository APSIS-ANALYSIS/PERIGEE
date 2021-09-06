#ifndef PRESTRESS_SOLID_HPP
#define PRESTRESS_SOLID_HPP
// ============================================================================
// Prestress_solid.hpp
//
// This class stores the values of the pre-stresses at quadrature points. 
// Pre-stress contains 6 components for the symmetric stress tensor.
//
// Date: Oct. 10 2017
// Author: Ju Liu
// ============================================================================
#include "ALocal_Elem_wTag.hpp"

class Prestress_solid
{
  public:
    Prestress_solid( const ALocal_Elem * const &locelem, const int &in_nqp_tet,
       const std::string &in_ps_fName = "prestress" );

    virtual ~Prestress_solid();

    // ------------------------------------------------------------------------
    // Input: ee the element index
    // Output: esval is the vector holding 6 components of the prestress
    // ------------------------------------------------------------------------
    virtual std::vector<double> get_prestress(const int &ee, const int &ii ) const;

    // ------------------------------------------------------------------------
    // Input: ee the element index
    //        0 <= ii < nqp is the quadrature point index
    //        in_esval is the element's prestress value at the quadrature point.
    // Users are responsible for making sure that the vector has lenght 6*nqp
    // ------------------------------------------------------------------------
    virtual void set_prestress(const int &ee, const int &ii, const double * const &in_esval );

    // ------------------------------------------------------------------------
    // record the prestress values to a h5 file
    // ------------------------------------------------------------------------
    virtual void write_prestress_hdf5() const;

    // ------------------------------------------------------------------------
    // load the prestress values from a h5 file
    // ------------------------------------------------------------------------
    virtual void read_prestress_hdf5() const;

    virtual void print_info() const;

  private:
    // ------------------------------------------------------------------------
    // The number of local element and the number of volumetric quadrature
    // points
    // ------------------------------------------------------------------------
    const int nlocalele, nqp;

    // ------------------------------------------------------------------------
    // The file name for storing the prestress values.
    // ------------------------------------------------------------------------
    const std::string ps_fileBaseName;

    // ------------------------------------------------------------------------
    // sval[ee][nqp*ii+jj] stores the ee-th solid element's pre-stress ii-th
    // component value at the jj-th quadrature point.
    // ii = 0, ..., 5 correspond to s1 s2 s3 s4 s5 s6
    //               s1  s2  s3
    //               s2  s4  s5
    //               s3  s5  s6
    // 0 <= jj < nqp
    // ------------------------------------------------------------------------
    std::vector< std::vector<double> > sval;
};

#endif
