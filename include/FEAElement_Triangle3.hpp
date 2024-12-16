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

    virtual ~FEAElement_Triangle3();

    virtual int get_elemDim() const {return 2;}

    virtual FEType get_Type() const {return FEType::Tri3;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 3;}

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

    virtual std::vector<double> get_d2R_dxx( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3::get_d2R_dxx function error.\n" );
      return { 0.0, 0.0, 0.0 };
    }

    virtual std::vector<double> get_d2R_dyy( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3::get_d2R_dyy function error.\n" );
      return { 0.0, 0.0, 0.0 };
    }

    virtual std::vector<double> get_d2R_dxy( const int &quaindex ) const
    {
      ASSERT( quaindex >= 0 && quaindex < numQuapts, "FEAElement_Triangle3::get_d2R_dxy function error.\n" );
      return { 0.0, 0.0, 0.0 };
    }

    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex, double * const &jac_value) const;

    virtual double get_detJac(const int &quaindex) const {return detJac;}

  private:
    const int numQuapts;

    double * R;

    // tri3 is linear element, hence the derivatives are constant
    double dR_dx[3], dR_dy[3];

    // Container for 
    // dx_dr : 0 <= ii < 4
    // dr_dx : 4 <= ii < 8
    double Jac[8];

    double detJac;
};

#endif
