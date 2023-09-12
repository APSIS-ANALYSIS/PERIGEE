#ifndef FEAELEMENT_TRIANGLE3_3D_DER0_HPP
#define FEAELEMENT_TRIANGLE3_3D_DER0_HPP
// ==================================================================
// FEAElement_Triangle3_3D_der0.hpp
// Element routine for the linear triangular element in
// three-dimensional space, with evaluation of the basis functions 
// only (no derivatives). 
// 
// Triangle3 means 3-node triangle; _3D means the element has
// coordinates in three-dimensions; _der0 means that only the function
// itself is evaluated.
//
// This class is designed for boundary integrations for elemental/
// natural boundary conditions in 3D problems. Although the element
// is a 2D element, the geometry has x-y-z coordinates.
// The Jacobian determinant is calculated as the norm of the cross
// product of the two in-plane vectors. Therefore, this is not the
// classical constant-strain triangle element.
//
// Date Created: Jan. 19 2017
// ==================================================================
#include "FEAElement.hpp"
#include "Math_Tools.hpp"

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

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 3;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual std::vector<double> get_R( const int &quaindex ) const;

    // Assumes the triangle nodes are arranged such that the outward
    // direction is given by dx_dr x dx_ds
    virtual Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const;

    // If the triangle nodes are NOT arranged in any particular order,
    // use an interior node to define the outward direction.
    virtual Vector_3 get_normal_out( const int &quaindex,
        const double &sur_pt_x, const double &sur_pt_y, const double &sur_pt_z,
        const double &intpt_x, const double &intpt_y, const double &intpt_z,
        double &area ) const;

    virtual double get_detJac(const int &quaindex) const {return detJac;}

  private:
    const int numQuapts;

    // container for R0 = 1 - r - s, R1 = r, R2 = s :
    // 0 <= ii < 3 x numQuapts
    double * R;

    // unit outward normal vector
    Vector_3 un;
    
    // Jacobian determinant 
    double detJac;
};

#endif
