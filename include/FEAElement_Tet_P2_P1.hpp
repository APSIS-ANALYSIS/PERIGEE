#ifndef FEAELEMENT_TET_P2_P1_HPP
#define FEAELEMENT_TET_P2_P1_HPP
// ==================================================================
// FEAElement_Tet_P2_P1.hpp
//
// This is an implementation of the Taylor-Hoold P2-P1 element.
//
//                     u
//                     ^
//                     |
//                     3(3)
//                    /| -
//                     |   -
//                  /  |     -
//                     7        8
//               /    |           -
//                9    |             -
//                     |                - 
//              /    /0 (0)----6---------2(2)----> s
//              /    /                -
//              /   /              -
//                 4           -   
//             /  /        5    
//             / /     -
//              /  - 
//             1(1)
//            /
//           *
//           r
// There are 10-node associated with the displacement/velocity fields,
// and 4 nodes associated with the pressure field.
//
// Author: Ju Liu
// Date Created: Feb. 19 2018
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_Tet_P2_P1 : public FEAElement
{
  public:
    FEAElement_Tet_P2_P1( const int &in_nqua );
    
    virtual ~FEAElement_Tet_P2_P1();

    virtual int get_elemDim() const {return 3;}

    virtual int get_Type() const {return 532;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 14;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Return the element size, by treating the element as a 4-node
    // element, and use the algorithm for linear tet.
    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const;

    // In the get_R, get_gradR, get_R_gradR functions, we will output all
    // 14 basis functions. Users should know that the first 10 are P2 disp
    // field, and the next 4 are P1 pressure field
    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_z ) const;

    virtual void get_Jacobian(const int &quaindex,
        double * const &jac_value) const
    {
      for(int ii=0; ii<9; ++ii) jac_value[ii] = dx_dr[9*quaindex+ii];
    }

    virtual void get_invJacobian(const int &quaindex,
        double * const &jac_value) const
    {
      for(int ii=0; ii<9; ++ii) jac_value[ii] = dr_dx[9*quaindex+ii];
    }

    virtual double get_detJac(const int &quaindex) const
    {return detJac[quaindex];}

  private:
    const int numQuapts;

    // Basis function and first derivatives
    // 0 <= ii < 14 * numQuapts
    double * R, * dR_dx, * dR_dy, * dR_dz;

    // Auxiliary variable for basis function evaluation
    double dR_dr [14], dR_ds [14], dR_dt [14];

    // Geometrical info.
    // 0 <= ii < 9 numQuapts
    double * dx_dr, * dr_dx;

    // detJac : 0 <= ii < numQuapts
    double * detJac;
};

#endif
