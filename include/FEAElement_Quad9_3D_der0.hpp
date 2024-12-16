#ifndef FEAELEMENT_QUAD9_3D_DER0_HPP
#define FEAELEMENT_QUAD9_3D_DER0_HPP
// ==================================================================
// FEAElement_Quad9_3D_der0.hpp
// Element routine for the biquadratic quadrilateral element in 
// three-dimensional space, with evaluation of the basis functions
// only (no derivatives).
//
//     s
//     |
//     |
//     3--------6---------2
//     |                  |
//     |                  |
//     |                  |
//     7        8         5
//     |                  |
//     |                  |
//     |                  |
//     0--------4-------- 1 -- r
//
//
// Quad4 means 9-node quad; _3D means the element has
// coordinates in three-dimensions; _der0 means that only the function
// itself is evaluated.
//
// This class is designed for boundary integrations for elemental/
// natural boundary conditions in 3D problems.
//
// Author: Ju Liu
// Date Created: Sep. 12 2023.
// ==================================================================
#include "FEAElement.hpp"
#include "Math_Tools.hpp"
#include "Vector_3.hpp"

class FEAElement_Quad9_3D_der0 final : public FEAElement
{
  public:
    FEAElement_Quad9_3D_der0( const int &in_nqua );

    virtual ~FEAElement_Quad9_3D_der0();

    virtual int get_elemDim() const {return 2;}

    virtual FEType get_Type() const {return FEType::Quad9_der0;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 9;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual std::vector<double> get_R( const int &quaindex ) const;

    // Assuming the quad nodes are arranged such that the outward
    // direction is given by dx_dr x dx_ds
    virtual Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const;
    
    // If the quad nodes are NOT arranged in any particular order,
    // use an interior node to define the outward direction. 
    virtual Vector_3 get_normal_out( const int &quaindex, const Vector_3 &sur_pt,
        const Vector_3 &int_pt, double &len ) const;

    virtual double get_detJac(const int &quaindex) const {return detJac[quaindex];}

  private:
    const int numQuapts;

    // Container for R
    // 0 <= ii < 9 x numQuapts
    double * R;

    // unit normal vector components, each of length numQuapts
    std::vector<Vector_3> un;

    // Jacobian determinant, length numQuapts
    double * detJac; 
};

#endif
