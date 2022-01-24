#ifndef FEAELEMENT_MINI_P1_3D_HPP
#define FEAELEMENT_MINI_P1_3D_HPP
// ==================================================================
// FEAElement_MINI_P1_3D.hpp
//
// This is an implementation of the element routine for 3D MINI
// element based on P1 tetrahedral element.
// 
// The enriched quartic bubble is node 4.
//
// Author: Ju Liu
// Date Created: Feb. 1 2018
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_MINI_P1_3D : public FEAElement
{
  public:
    FEAElement_MINI_P1_3D( const int &in_nqp );

    virtual ~FEAElement_MINI_P1_3D();

    virtual int get_elemDim() const {return 3;}

    virtual int get_Type() const {return 571;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 5;}

    virtual void print() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Return the element size
    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const;

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_z ) const;

    virtual void get_Jacobian(const int &quaindex,
        double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex,
        double * const &jac_value) const;

    virtual double get_detJac(const int &quaindex) const
    {return detJac;}

  private:
    const int numQuapts;

    // R : 0 <= ii < 5 x numQuapts
    double * R;

    // tet4 basis are linear, the first derivatives are constant
    double dR_dx[4];
    double dR_dy[4];
    double dR_dz[4];

    // Bubble enriched basis derivatives, length : numQuapts
    double * dB_dx;
    double * dB_dy; 
    double * dB_dz; 
    
    // Container for 
    // dx_dr : 0 <= ii < 9
    // dr_dx : 9 <= ii < 18
    double Jac[18];

    double detJac;
};

#endif
