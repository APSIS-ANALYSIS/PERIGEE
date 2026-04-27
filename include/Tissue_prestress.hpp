#ifndef TISSUE_PRESTRESS_HPP
#define TISSUE_PRESTRESS_HPP
// ============================================================================
// Tissue_prestress.hpp
//
// This class stores the values of the pre-stresses at quadrature points. 
// Pre-stress contains 6 components for the symmetric stress tensor.
//
// Date: Oct. 10 2017
// Author: Ju Liu
// ============================================================================
#include "HDF5_Writer.hpp"
#include "HDF5_Reader.hpp"
#include "ALocal_Elem.hpp"

class Tissue_prestress
{
  public:
    Tissue_prestress( const ALocal_Elem * const &locelem, const int &in_nqp_tet,
       const int &in_cpu_rank, const bool &load_from_file,
       const std::string &in_ps_fName = "prestress" );

    virtual ~Tissue_prestress() = default;

    // ------------------------------------------------------------------------
    // Input: ee the element index
    // Output: the vector holding 6 components of the prestress at all quad pts
    // ------------------------------------------------------------------------
    virtual std::vector<double> get_prestress( const int &ee ) const;

    virtual std::array<double,6> get_prestress( const int &ee, const int &qua ) const;

    // ------------------------------------------------------------------------
    // Input: ee the element index
    //        in_esval is the element's prestress value at the quadrature point.
    // Users are responsible for making sure that the vector has lenght 6*nqp
    // The value of in_esval is ADDED into the data
    // ------------------------------------------------------------------------
    virtual void add_prestress( const int &ee, const double * const &in_esval );

    // ------------------------------------------------------------------------
    // Input: ee the element index 0 <= ee < nlocalele
    //        qua the quadrature point index 0 <= qua < nqp
    //        in_esval the element's qua-th quadrature point's stress components
    //        with length being 6.
    // The value of in_esval is ADDED into the data
    // ------------------------------------------------------------------------
    virtual void add_prestress( const int &ee, const int &qua,
        const double * const &in_esval );

    // ------------------------------------------------------------------------
    // record the prestress values to a h5 file as a row data
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
    // qua_prestress[ee][6*ii+jj] stores the ee-th element's pre-stress 
    // jj-th component in Voigt notation at the ii-th quadrature point.
    // 0 <= ee < nlocalele, 0 <= ii < nqp, 0 <= jj < 6 if the element is a solid
    // element; otherwise, qua_prestress[ee] is cleared for a fluid element.
    // 
    // Size is nlocalele x (nqp x 6)
    // ------------------------------------------------------------------------
    std::vector< std::vector<double> > qua_prestress;
};

#endif
