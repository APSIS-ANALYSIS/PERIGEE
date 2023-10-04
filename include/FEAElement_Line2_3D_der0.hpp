#ifndef FEAELEMENT_LINE2_3D_DER0_HPP
#define FEAELEMENT_LINE2_3D_DER0_HPP
// ==================================================================
// FEAElement_Line2_3D_der0.hpp
//
// This is an implementation of the element routine for 1D P1 
// element with basis function value evaluated only. This routine is
// designed for the evaulation of Neumann boundary conditions.
// 
// Node index :  0 --------------- 1
//
// Author: Ju Liu
// Date created: Nov. 27 2017
// ==================================================================
#include "FEAElement.hpp"
#include "FE_Tools.hpp"

class FEAElement_Line2_3D_der0 : public FEAElement
{
  public:
    FEAElement_Line2_3D_der0( const int &in_nqua );

    virtual ~FEAElement_Line2_3D_der0();

    virtual int get_elemDim() const {return 1;}

    virtual int get_Type() const {return 510;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 2;}

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
    {return detJac;};

  private:
    const int numQuapts;

    // length nLocBas x numQuapts = 2 x numQuapts
    double * R;

    // length should be numQuapts, 
    // here for linear element, they are all constant
    Vector_3 dx_dr;
    double detJac;
};

#endif
