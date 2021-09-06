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
    Prestress_solid( const ALocal_Elem * const &locelem, const int &in_nqp_tet );

    virtual ~Prestress_solid();

    // Input: ee the element index
    // Output: esval is the vector holding 6 components of the prestress
    virtual std::vector<double> get_sval(const int &ee, const int &ii ) const;

    // Input: ee the element index
    //        in_esval is the incremental of the element's prestress values.
    //        sval[ee][ ] += in_esval
    // Users are responsible for making sure that the vector has lenght
    // 7*nqp
    virtual void update_sval(const int &ee, const std::vector<double> &in_esval );

    virtual void print_info() const;

  private:
    const int nlocalele;
    const int nqp;

    // sval[ee][nqp*ii+jj] stores the ee-th element's pre-stress ii-th
    // component value at the jj-th quadrature point.
    // ii = 0, ..., 5 correspond to s1 s2 s3 s4 s5 s6
    //               s1  s2  s3
    //               s2  s4  s5
    //               s3  s5  s6
    // 0 <= jj < nqp
    //
    // if elem_tag[ee] == 0, then the sval has zero length because fluid
    // does not need pre-stress;
    // if elem_tag[ee] == 1, reserve nqp * 7 double slots for the element
    std::vector< std::vector<double> > sval;
};

#endif
