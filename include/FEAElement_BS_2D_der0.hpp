#ifndef FEAELEMENT_BS_2D_DER0_HPP
#define FEAELEMENT_BS_2D_DER0_HPP
// ============================================================================
// FEAElement_BS_2D_der0.hpp
// This is an implementation of the element routine for 2D B-splines.
// The basis function is evaluated only for this der0 case.
//
// This is an implementation for no-cache-style programming (i.e., we do not
// cache the quadrature info in memory). The objective is to reduce the amount
// of heap usage and reduce the call of malloc/new for the memory pool purpose.
// Another design factor is to make the implementation faster than the previous
// version.
//
// Date: Sep. 15 2015
// ============================================================================
#include "FEAElement.hpp"

class FEAElement_BS_2D_der0 : public FEAElement
{
  public:
    FEAElement_BS_2D_der0( const int &in_sdeg, const int &in_tdeg,
       const int &in_nquas, const int &in_nquat );

    virtual ~FEAElement_BS_2D_der0();

    virtual int get_elemDim() const {return 2;}

    virtual int get_Type() const {return 620;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual void buildBasis( const double &hx, const double &hy,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ext_x,
        const double * const &ext_y );

    virtual void get_R(const int &quaindex, double * const &basis) const;

    virtual double get_detJac(const int &quaindex) const;

    virtual void reset_degree(const int &new_sdeg, const int &new_tdeg);

    virtual void reset_numQua( const int &new_squa, const int &new_tqua );

    virtual void print() const;

    virtual double get_memory_usage() const;

  private:
    int num_qua_s, num_qua_t, numQuapts;
    int sdeg, tdeg, sdp1, tdp1; 
    int nLocBas;
    int rlength; // length of R equals nLocBas * numQuapts

    // R and detJac are stored in the single double vector
    // R : 0<= ii < nLocBas * numQuapts
    // detJac : nLocBas * numQuapts <= ii < nLocBas*numQuapts + numQuapts
    double * R;

    // container for the B-splines, Bezier elements, and their derivatives
    // dR_ds: 0 <= ii < nLocBas
    // dR_dt: nLocBas <= ii < 2*nLocBas
    double * dRr;
    
    // Nns : 0 <= ii < sdp1
    // dNns_ds : sdp1 <= ii < 2*sdp1
    // Nnt : 0 <= ii < tdp1
    // dNnt_dt : tdp1 <= ii < 2*tdp1
    double * Nns, * Nnt;

    void BuildShape_atQua( const int &quaindex,
        const double &hx, const double &hy,
        const double * const &ctrl_x,
        const double * const &ctrl_y );

    virtual void clearBasisCache();

    virtual void resize_container();
};

#endif
