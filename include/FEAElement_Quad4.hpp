#ifndef FEAELEMENT_QUAD4_HPP
#define FEAELEMENT_QUAD4_HPP
// ==================================================================
// FEAElement_Quad4.hpp
// This is an implementation of the element routine for bilinear
// quadrilateral element in 2D.
//
// This class is designed mainly for the 2D FEM assembly.
//
// Date created: Sep. 2023
// ==================================================================
#include "FEAElement.hpp"
#include "FE_Tools.hpp"

class FEAElement_Quad4 final : public FEAElement
{
  public:
    FEAElement_Quad4( const int &in_nqua );

    virtual ~FEAElement_Quad4();

    virtual int get_elemDim() const {return 2;}

    virtual FEType get_Type() const {return FEType::Quad4;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 4;}

    virtual void print_info() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y );

    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y ) const;

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual std::vector<double> get_R( const int &quaindex ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y ) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y ) const;

    virtual std::vector<double> get_dR_dx( const int &quaindex ) const;

    virtual std::vector<double> get_dR_dy( const int &quaindex ) const;

    virtual void get_2D_R_dR_d2R( const int &quaindex,
        double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy,
        double * const &basis_xy ) const;

    virtual std::vector<double> get_d2R_dxx( const int &quaindex ) const;

    virtual std::vector<double> get_d2R_dyy( const int &quaindex ) const;

    virtual std::vector<double> get_d2R_dxy( const int &quaindex ) const;

    virtual void get_Jacobian(const int &quaindex,
        double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex,
        double * const &jac_value) const;

    virtual double get_detJac(const int &quaindex) const
    {return Jac[8*numQuapts + quaindex];}

  private:
    const int numQuapts;

    // length 4 x numQuapts
    double * R, * dR_dx, * dR_dy;
    double * d2R_dxx, * d2R_dyy, * d2R_dxy;

    // length 9 x numQuapts
    // dx_ds : 0             <= ii < 4 * numQuapts
    // ds_dx : 4 * numQuapts <= ii < 8 * numQuapts
    // detJac: 8 * numQuapts <= ii < 9 * numQuapts
    double * Jac;
};

#endif
