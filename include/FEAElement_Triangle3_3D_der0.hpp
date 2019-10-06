#ifndef FEAELEMENT_TRIANGLE3_3D_DER0_HPP
#define FEAELEMENT_TRIANGLE3_3D_DER0_HPP
// ==================================================================
// FEAElement_Triangle3_3D_der0.hpp
// This is an implementation of the element routine for linear
// triangle element in three-dimensional space, with basis function 
// value evaluated only. 
// 
// The Triangle3 represents 3-node triangle; _3D means that this
// element has coordinates in three-dimensions; _der0 means that only
// the function value is evaluated (this is a function for Natural BC).
//
// This class is designed for the purpose of boundary integration for
// natural boundary conditions in 3D problems. Although the element
// is a 2D element, the coordinates for the geometry has x-y-z 
// components. The Jacobian of the area is calculated by the cross
// product of the two in-plane vectors. Therefore, this is not the
// classical constatn-strain triangle element.
//
// Date Created: Jan. 19 2017
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_Triangle3_3D_der0 : public FEAElement
{
  public:
    FEAElement_Triangle3_3D_der0( const int &in_nqua );

    virtual ~FEAElement_Triangle3_3D_der0();

    virtual int get_elemDim() const {return 2;}

    // element type : 521
    // 5: simplicial element
    // 2: 2D element
    virtual int get_Type() const {return 521;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual void reset_numQuapts( const int &new_num_qua );

    virtual void print() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    // Assuming the triangle points has been orgainzed so that the outward
    // direction is given by dx_dr x dx_ds
    virtual void get_2d_normal_out( const int &quaindex,
        double &nx, double &ny, double &nz, double &area ) const;

    // Assuming the triangle points are NOT organized. Use an interior point
    // to define the outward direction.
    virtual void get_normal_out( const int &quaindex,
        const double &sur_pt_x, const double &sur_pt_y, const double &sur_pt_z,
        const double &intpt_x, const double &intpt_y, const double &intpt_z,
        double &nx, double &ny, double &nz, double &len ) const;

    virtual double get_detJac(const int &quaindex) const
    {return detJac;}

  private:
    const int nLocBas;

    int numQuapts;

    // container for R0 = 1 - r - s, R1 = r, R2 = s :
    // 0 <= ii < 3 x numQuapts
    double * R;

    // container for dx_dr, dx_ds, dy_dr, dy_ds, dz_dr, dz_ds :
    double dx_dr, dx_ds, dy_dr, dy_ds, dz_dr, dz_ds;

    // unit outward normal vector
    double unx, uny, unz;

    // detJac 
    double detJac;

    // Free the space
    virtual void clearBasisCache();

    // resize the dynamic arrays R, dx_ds, detJac if the quadrature rule
    // is changed
    void resize_container(); 
};

#endif
