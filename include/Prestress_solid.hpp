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
#include "HDF5_Writer.hpp"
#include "ALocal_Elem.hpp"

class Prestress_solid
{
  public:
    Prestress_solid( const ALocal_Elem * const &locelem, const int &in_nqp_tet,
       const int &in_cpu_rank, const std::string &in_ps_fName = "prestress" );

    virtual ~Prestress_solid();

    // ------------------------------------------------------------------------
    // Input: ee the element index
    // Output: the vector holding 6 components of the prestress at all quad pts
    // ------------------------------------------------------------------------
    virtual std::vector<double> get_prestress(const int &ee ) const;

    // ------------------------------------------------------------------------
    // Input: ee the element index
    //        0 <= ii < nqp is the quadrature point index
    //        in_esval is the element's prestress value at the quadrature point.
    // Users are responsible for making sure that the vector has lenght 6*nqp
    // ------------------------------------------------------------------------
    virtual void set_prestress(const int &ee, const int &ii, 
        const double * const &in_esval );

    // ------------------------------------------------------------------------
    // record the prestress values to a h5 file
    // ------------------------------------------------------------------------
    virtual void write_prestress_hdf5() const;

    virtual void print_info() const;

  private:
    // ------------------------------------------------------------------------
    // The rank or id of the subdomain
    // ------------------------------------------------------------------------
    const int cpu_rank;
    
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
    // qua_prestress[ee][nqp*ii+jj] stores the ee-th solid element's pre-stress 
    // jj-th component in Voigt notation at the ii-th quadrature point.
    // 0 <= ee < nlocalele, 0 <= ii < nqp, 0 <= jj < 6
    // ------------------------------------------------------------------------
    std::vector< std::vector<double> > qua_prestress;
};

#endif
