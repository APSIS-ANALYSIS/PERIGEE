#ifndef FEAELEMENT_TET4_HPP
#define FEAELEMENT_TET4_HPP
// ==================================================================
// FEAElement_Tet4.hpp
// This is an implementation of the element routine for linear
// tetrahedral element, with derivatives evaluated.
// 
// Tet4 means 4-node element, aka, linear tets.
//
// This class is designed for the volume integration in model assembly.
//
// Date Created: Jan 19 2017
// ==================================================================
#include "FEAElement.hpp"
#include "Matrix_double_3by3_Array.hpp"

class FEAElement_Tet4 : public FEAElement
{
  public:
    FEAElement_Tet4( const int &in_nqua );

    virtual ~FEAElement_Tet4();

    virtual int get_elemDim() const {return 3;}

    virtual int get_Type() const {return 501;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return 4;}

    virtual void reset_numQuapts( const int &new_num_qua );

    virtual void print() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z );

    // Returns the element size.
    // For the linear tets element, we calculate the DIAMETER of the
    // circumscribing sphere
    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const;

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z ) const;

    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_z ) const;

    virtual void get_3D_R_dR_d2R( const int &quaindex, 
        double * const &basis, double * const &basis_x, 
        double * const &basis_y, double * const &basis_z,
        double * const &basis_xx, double * const &basis_yy, 
        double * const &basis_zz, double * const &basis_xy, 
        double * const &basis_xz, double * const &basis_yz ) const;

    virtual void get_3D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, 
        double * const &basis_y, double * const &basis_z, 
        double * const &basis_xx, double * const &basis_yy, 
        double * const &basis_zz ) const;

    virtual void get_Jacobian(const int &quaindex, 
        double * const &jac_value) const;

    virtual void get_invJacobian(const int &quaindex, 
        double * const &jac_value) const;

    virtual double get_detJac(const int &quaindex) const
    {return detJac;}

  private:
    int numQuapts;

    // R : 0 <= ii < 4 x numQuapts
    double * R;

    // tet4 is linear, the first-order derivatives are constant
    double dR_dx[4];
    double dR_dy[4];
    double dR_dz[4];

    // Container for
    // dx_dr : 0 <= ii < 9
    // dr_dx : 9 <= ii < 18
    double Jac[18]; 

    double detJac;

    // Free the space of dynamic array
    virtual void clearBasisCache();

    // resize the R array if the quadrature rule is changed
    void resize_container();
};


#endif
