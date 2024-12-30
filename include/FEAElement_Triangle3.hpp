#ifndef FEAELEMENT_TRIANGLE3_HPP
#define FEAELEMENT_TRIANGLE3_HPP
// ==================================================================
// FEAElement_Triangle3.hpp
// This is an implementation of the element routine for linear 
// triangle element in 2D.
//
// This class is designed mainly for the 2D FEM assembly.
//
// Date created: Nov. 22 2017
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_Triangle3 final : public FEAElement
{
  public:
    FEAElement_Triangle3( const int &in_nqua );

    ~FEAElement_Triangle3() override = default;

    int get_elemDim() const override {return 2;}

    FEType get_Type() const override {return FEType::Tri3;}

    int get_numQuapts() const override {return numQuapts;}

    int get_nLocBas() const override {return nLocBas;}

    void print_info() const override;

    void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y ) override;

    double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y ) const override;

    void get_R( const int &quaindex, double * const &basis ) const override;
    
    std::vector<double> get_R( const int &quaindex ) const override;

    void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y ) const override;

    void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y ) const override;

    void get_2D_R_dR_d2R( const int &quaindex, 
        double * const &basis, 
        double * const &basis_x, double * const &basis_y, 
        double * const &basis_xx, double * const &basis_yy, 
        double * const &basis_xy ) const override;

    std::array<double,4> get_Jacobian_2D(const int &quaindex) const override;

    std::array<double,4> get_invJacobian_2D(const int &quaindex) const override;

    double get_detJac(const int &quaindex) const override {return detJac;}

  private:
    static constexpr int nLocBas = 3;
    const int numQuapts;

    std::vector<double> R {};

    // tri3 is linear element, hence the derivatives are constant
    std::array<double, 3> dR_dx, dR_dy;

    // Container for 
    // dx_dr : 0 <= ii < 4
    // dr_dx : 4 <= ii < 8
    std::array<double, 8> Jac;

    double detJac;
};

#endif
