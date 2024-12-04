#ifndef FEAELEMENT_QUAD9_HPP
#define FEAELEMENT_QUAD9_HPP
// ==================================================================
// FEAElement_Quad9.hpp
// This is an implementation of the element routine for biquadratic
// quadrilateral element in 2D.
//
// This class is designed mainly for the 2D FEM assembly.
//
// Date created: Sep. 2023
// ==================================================================
#include "FEAElement.hpp"
#include "FE_Tools.hpp"

class FEAElement_Quad9 : public FEAElement
{
  public:
    FEAElement_Quad9( const int &in_nqua );

    ~FEAElement_Quad9() override {}

    int get_elemDim() const override {return 2;}

    FEType get_Type() const override {return FEType::Quad9;}

    int get_numQuapts() const override {return numQuapts;}

    int get_nLocBas() const override {return 9;}

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

    std::vector<double> get_dR_dx( const int &quaindex ) const override;

    std::vector<double> get_dR_dy( const int &quaindex ) const override;

    void get_2D_R_dR_d2R( const int &quaindex,
        double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy,
        double * const &basis_xy ) const override;

    void get_Jacobian(const int &quaindex,
        double * const &jac_value) const override;

    void get_invJacobian(const int &quaindex,
        double * const &jac_value) const override;

    double get_detJac(const int &quaindex) const override
    {return Jac[8*numQuapts + quaindex];}

  private:
    const int numQuapts;

    // length 9 x numQuapts
    std::vector<double> R {};
    std::vector<double> dR_dx {};
    std::vector<double> dR_dy {};
    std::vector<double> d2R_dxx {};
    std::vector<double> d2R_dyy {};
    std::vector<double> d2R_dxy {};

    // length 9 x numQuapts
    // dx_ds : 0             <= ii < 4 * numQuapts
    // ds_dx : 4 * numQuapts <= ii < 8 * numQuapts
    // detJac: 8 * numQuapts <= ii < 9 * numQuapts
    std::vector<double> Jac {};
};

#endif
