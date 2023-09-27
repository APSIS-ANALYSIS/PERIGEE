#ifndef FEAELEMENT_LINE3_3D_DER0_HPP
#define FEAELEMENT_LINE3_3D_DER0_HPP
// ==================================================================
// FEAElement_Line3_3D_der0.hpp
//
// This is an implementation of the element routine for 1D P2 element
// with basis function value evaluated only. This routine is designed
// for the Neumann boundary condition implementation.
//
// Node index :  0 ----- 2 ----- 1
// Ref. coor  : 0.0     0.5     1.0
// 
// N0 = (2r-1)(r-1)
// N1 = r(2r-1)
// N2 = 4r(1-r)
//
// Author: Ju Liu
// Date Created: Nov. 28 2017
// ==================================================================
#include "FEAElement.hpp"
#include "FE_Tools.hpp"

class FEAElement_Line3_3D_der0 : public FEAElement
{
  public:
    FEAElement_Line3_3D_der0( const int &in_nqua );

    virtual ~FEAElement_Line3_3D_der0();

    virtual int get_elemDim() const {return 1;}

    virtual int get_Type() const {return 511;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 3;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x, const double * const &ctrl_y,
        const double * const &ctrl_z );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual Vector_3 get_normal_out( const int &quaindex,
        const std::vector<Vector_3> &ctrl_pt,
        const Vector_3 &int_pt, double &len ) const;

    virtual double get_detJac(const int &quaindex) const
    {return detJac[quaindex];};

  private:
    const int numQuapts;

    // length is nLocBas x numQuapts = 3 x numQuapts
    double * R;

    // length is numQuapts
    std::vector<Vector_3> dx_dr;
    double * detJac;
};

#endif
