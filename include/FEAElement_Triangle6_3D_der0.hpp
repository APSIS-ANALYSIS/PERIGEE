#ifndef FEAELEMENT_TRIANGLE6_3D_DER0_HPP
#define FEAELEMENT_TRIANGLE6_3D_DER0_HPP
// ==================================================================
// FEAElement_Triangle6_3D_der0.hpp
// Element routine for the quadratic triangular element in 
// three-dimensional space, with evaluation of the basis functions
// only (no derivatives).
//
//     s
//     |
//     2
//     | -
//     |   -
//     |     -
//     |       -      t = 1 - r - s
//     5         4
//     |           -
//     |             -
//     |               -
//     0---------3------- 1 -- r
//
//
// Triangle6 means 6-node triangle; _3D means the element has
// coordinates in three-dimensions; _der0 means that only the function
// itself is evaluated.
//
// This class is designed for boundary integrations for elemental/
// natural boundary conditions in 3D problems.
//
// Author: Ju Liu
// Date Created: Feb. 17 2018.
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_Triangle6_3D_der0 final : public FEAElement
{
  public:
    FEAElement_Triangle6_3D_der0( const int &in_nqua );

    ~FEAElement_Triangle6_3D_der0() override = default;

    int get_elemDim() const override {return 2;}

    FEType get_Type() const override {return FEType::Tri6_der0;}

    int get_numQuapts() const override {return numQuapts;}

    int get_nLocBas() const override {return 6;}

    void print_info() const override;

    void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) override;

    void get_R( const int &quaindex, double * const &basis ) const override;

    std::vector<double> get_R( const int &quaindex ) const override;

    // Assuming the triangle nodes are arranged such that the outward
    // direction is given by dx_dr x dx_ds
    Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const override;
    
    // If the triangle nodes are NOT arranged in any particular order,
    // use an interior node to define the outward direction. 
    Vector_3 get_normal_out( const int &quaindex, const Vector_3 &sur_pt,
        const Vector_3 &int_pt, double &len ) const override;

    double get_detJac(const int &quaindex) const override {return detJac[quaindex];}

  private:
    const int numQuapts;

    // Container for R0, R1, R2, R3, R4, R5
    // 0 <= ii < 6 x numQuapts
    std::vector<double> R {};

    // unit normal vector components, each of length numQuapts
    std::vector<Vector_3> un;

    // Jacobian determinant, length numQuapts
    std::vector<double> detJac {};
};

#endif
