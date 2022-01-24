#ifndef FEAELEMENT_MINI_P1_2D_HPP
#define FEAELEMENT_MINI_P1_2D_HPP
// ==================================================================
// FEAElement_MINI_P1_2D.hpp
//
// This is an implementation of the element routine for 2D MINI
// element based on P1 triangular element.
//
//   s
//   2
//   |-
//   |  -
//   |    -
//   |      -
//   |   3    -
//   |          -
// 0 ------------- 1  r
//  
//  N0 = t = 1 - r - s
//  N1 = r
//  N2 = s
//  N3 = rst = 27rs(1-r-s)
//
// Author: Ju Liu
// Date Created: Oct. 29 2017
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_MINI_P1_2D : public FEAElement
{
  public:
    FEAElement_MINI_P1_2D( const int &in_nqp );

    virtual ~FEAElement_MINI_P1_2D();

    virtual int get_elemDim() const {return 2;}

    virtual int get_Type() const {return 561;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 4;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x, const double * const &ctrl_y );

    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y ) const;

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y ) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y ) const;

    virtual void get_Jacobian(const int &quaindex,
        double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex,
        double * const &jac_value) const;

    virtual double get_detJac(const int &quaindex) const {return detJac;}

  private:
    const int numQuapts;

    // length : 4 x numQuapts
    double * R;

    // Linear basis function derivatives
    double dR_dx[3];
    double dR_dy[3];

    // Bubble enriched basis derivatives
    // length : numQuapts
    double * dB_dx;
    double * dB_dy;

    // The geometry is constant straing, hence these are constants.
    // dx_dr : 0 <= ii < 4
    // dr_dx : 4 <= ii < 8
    double Jac[8];
    
    double detJac;
};

#endif
