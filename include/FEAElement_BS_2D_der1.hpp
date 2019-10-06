#ifndef FEAELEMENT_BS_2D_DER1_HPP
#define FEAELEMENT_BS_2D_DER1_HPP
// ============================================================================
// FEAElement_BS_2D_der1.hpp
// This is an implementation of the element routine for 2D B-splines. The basis
// function are evaluated for the function itself and 1st order derivatives.
// 
// The Jacobian and the inverse of the Jacobian are evaluated.
//
// The normal vector in four directions can be computed from the Jacobian.
//
// This design is intended for no-cache-style programming. The objective is to
// reduce the heap usage and reduce the call of malloc/new. In addition, this is
// also a design for faster element routiens.
//
// Date: Sept 24 2015
// ============================================================================

#include "FEAElement.hpp"

class FEAElement_BS_2D_der1 : public FEAElement
{
  public:
    FEAElement_BS_2D_der1( const int &in_sdeg, const int &in_tdeg,
       const int &in_nquas, const int &in_nquat );

    virtual ~FEAElement_BS_2D_der1();

    virtual int get_elemDim() const {return 2;}

    virtual int get_Type() const {return 621;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual void get_R(const int &quaindex, double * const &basis) const;

    virtual double get_detJac(const int &quaindex) const;

    virtual void get_R_gradR(const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y) const;

    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex, double * const &jac_value) const;

    virtual void get_2d_normal_front( const int &quaindex,
                double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_back( const int &quaindex,
                double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_left( const int &quaindex,
                double &nx, double &ny, double &line ) const;

    virtual void get_2d_normal_right( const int &quaindex,
                double &nx, double &ny, double &line ) const;

    virtual void reset_degree(const int &new_sdeg, const int &new_tdeg);

    virtual void reset_numQua( const int &new_squa, const int &new_tqua );

    virtual void print() const;
    
    virtual double get_memory_usage() const;

    virtual void buildBasis( const double &hx, const double &hy,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ext_x,
        const double * const &ext_y );


  private:
    int num_qua_s, num_qua_t, numQuapts;
    int sdeg, tdeg, sdp1, tdp1;
    int nLocBas;
    int rlength;
    
    // container for R, dR_dx, dR_dy
    // R      : 0 <= ii < rlength
    // dR_dx  : rlength   <= ii < 2 * rlength
    // dR_dy  : 2*rlength <= ii < 3 * rlength
    double * R;

    // container for dx_ds, ds_dx, detJac
    // dx_ds :        0 <= ii < 4*numQua
    // ds_dx : 4*numQua <= ii < 8*numQua
    // detJac: 8*numQua <= ii < 9*numQua
    double * Jac;

    double * dRr;

    double * Nns, * Nnt;

    virtual void clearBasisCache();

    void BuildShape_atQua( const int &quaindex,
        const double &hx, const double &hy,
        const double * const &ctrl_x,
        const double * const &ctrl_y );

    void resize_container();
    
    void normalize_2d_vector(double &x, double &y) const
    {
      const double temp = 1.0 / sqrt(x*x+y*y);
      x = x * temp;
      y = y * temp;
    }

};


#endif
