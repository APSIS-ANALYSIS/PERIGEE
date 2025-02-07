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

class FEAElement_Triangle3_3D_der0 final : public FEAElement
{
  public:
    FEAElement_Triangle3_3D_der0( const int &in_nqua );

    ~FEAElement_Triangle3_3D_der0() override = default;

    int get_elemDim() const override {return 2;}

    FEType get_Type() const override {return FEType::Tri3_der0;}

    int get_numQuapts() const override {return numQuapts;}

    int get_nLocBas() const override {return nLocBas;}

    void print_info() const override;

    void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) override;

    void get_R( const int &quaindex, double * const &basis ) const override;

    std::vector<double> get_R( const int &quaindex ) const override;

    // Assumes the triangle nodes are arranged such that the outward
    // direction is given by dx_dr x dx_ds
    Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const override;

    // If the triangle nodes are NOT arranged in any particular order,
    // use an interior node to define the outward direction.
    Vector_3 get_normal_out( const int &quaindex, const Vector_3 &sur_pt,
        const Vector_3 &int_pt, double &area ) const override;

    double get_detJac(const int &quaindex) const override {return detJac;}

    // Return the derivatives of the physical coordinates with respect to the
    // element reference coordinate.
    // These functions are needed in the FE_T::search_closest_point function,
    // which is called inside the sliding interface formulation.
    Vector_3 get_dx_dr( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const override
    {
      return Vector_3( - ctrl_x[0] + ctrl_x[1],
                       - ctrl_y[0] + ctrl_y[1],
                       - ctrl_z[0] + ctrl_z[1] );
    }
    
    Vector_3 get_dx_ds( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const override
    {
      return Vector_3( - ctrl_x[0] + ctrl_x[2],
                       - ctrl_y[0] + ctrl_y[2],
                       - ctrl_z[0] + ctrl_z[2] );
    }

    Vector_3 get_d2x_drr( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const override
    {return Vector_3(0.0, 0.0, 0.0);}

    Vector_3 get_d2x_dss( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const override
    {return Vector_3(0.0, 0.0, 0.0);}

    Vector_3 get_d2x_drs( const int &quaindex,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const override
    {return Vector_3(0.0, 0.0, 0.0);}

  private:
    static constexpr int nLocBas = 3;
    
    const int numQuapts;

    // container for R0 = 1 - r - s, R1 = r, R2 = s :
    // 0 <= ii < 3 x numQuapts
    std::vector<double> R {};

    // unit outward normal vector
    Vector_3 un;
    
    // Jacobian determinant 
    double detJac {};
};

#endif
