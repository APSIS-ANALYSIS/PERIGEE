#ifndef FEAELEMENT_DIVB_2D_DER1_HPP
#define FEAELEMENT_DIVB_2D_DER1_HPP
// ==================================================================
// FEAElement_DivB_2D_der1.hpp
// This is an implementation for 2D Divergence-conforming B-Splines.
// ==================================================================

#include "FEAElement.hpp"

class FEAElement_DivB_2D_der1 : public FEAElement
{
  public:
    // Constructor
    // input: sdeg and tdeg are the highest polynomial order in s/t
    //        direction.
    FEAElement_DivB_2D_der1( const int &in_sdeg, const int &in_tdeg,
        const int &in_nquas, const int &in_nquat );

    virtual ~FEAElement_DivB_2D_der1();

  private:
    int num_qua_s, num_qua_t, numQuapts;
    int sdeg, tdeg, sdp1, tdp1;
    int nLocBas;
    int rlength;

    double * Rs, * Rt;

    double * Jac;
};

#endif
