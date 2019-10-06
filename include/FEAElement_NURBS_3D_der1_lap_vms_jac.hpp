#ifndef FEAELEMENT_NURBS_3D_DER1_LAP_VMS_JAC_HPP
#define FEAELEMENT_NURBS_3D_DER1_LAP_VMS_JAC_HPP
// ==================================================================
// FEAElement_NURBS_3D_der1_lap_vms_jac.hpp
// It is a finite element implementation of the physical 3D NURBS 
// element, with derivative up to laplacian. mixed second order
// derivaties are NOT evaluated! The dxi/dx and dx_dxi quantities are
// computed and cached. They are saved in ds_dx_qua, ds_dy_qua, ...
// du_dz_qua, dx_ds_qua, ... , dz_du_qua, which are vectors with 
// length numQuapts.
// The dxi/dx's will be used to calculate G_ij and g_i, which are
// further used to calculate the stabilization parameters tau_m and 
// tau_c.
//
// The dx_dxi's will be used to calculate the element surface normal
// vectors.
//
// Reference: CMAME 197 (2007) pp. 173-201.
//
// _lap: the Laplacian of the basis function will be calculated and 
//       cached.
// 
// _vms: version of variational multiscalue analysis. The element routine
//       is a direct inheritance of the _v3 version, which explores fast 
//       element implementation: This routine is based on univariate 
//       bernstein basis instead of parentElement.
// 
// _jac: the Jacobian matrix component will be calculated and cached.
//       the Jacobian matrix can be used to calculate the normal and 
//       tangential vector as well as the surface Jacobian determinant.
//
// Date:
// June 24 2015
// ==================================================================
#include <cmath>
#include <vector>
#include "Sys_Tools.hpp"
#include "Matrix_double_6by6_Array.hpp"
#include "ALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "FEANode.hpp"
#include "IAExtractor.hpp"
#include "BernsteinBasis_Array.hpp"
#include "FEAElement.hpp"

class FEAElement_NURBS_3D_der1_lap_vms_jac : public FEAElement
{
  public:
    FEAElement_NURBS_3D_der1_lap_vms_jac( const int &in_eIndex,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor,
        const ALocal_IEN * const &locIEN );
    
    
    FEAElement_NURBS_3D_der1_lap_vms_jac( const int &in_eIndex,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z,
        const double * const &ctrl_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );


    virtual ~FEAElement_NURBS_3D_der1_lap_vms_jac();

    virtual int get_elemIndex() const {return elem_index;}
    virtual int get_elemDim() const {return 3;}

    // Type = 304: 3D NURBS H1 conforming basis with up to 
    //             laplacian derivatives, dxi_dx, and dx_dxi.
    virtual int get_Type() const {return 304;}

    // numType = 1: This is a basis for scalar variable.
    virtual int get_numType() const {return 1;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual bool is_sizeNonzero() const {return is_sNonzero;}

    virtual void clearBasisCache();

    virtual void buildBasis( 
        const int &ssdp1, const int &ttdp1, const int &uudp1,
        const int &num_qs, const int &num_qt, const int &num_qu,
        const IALocal_meshSize * const &mSize, 
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor, 
        const ALocal_IEN * const &locIEN );

    // buildBasis : evaluate basis function values at quadrature points.
    // This function only relies on the local element-wise infomation.
    // \para bs/t/u : Bernstein polynomial value at quadrature points
    // \para ctrl_x/y/z/w : control points and weights
    // \para ext_x/y/z : Bezier extraction operator in x/y/z direction
    virtual void buildBasis( 
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z,
        const double * const &ctrl_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z );

    virtual int get_numQuapts() const {return numQuapts;}

    // Get function value R.
    virtual void get_R(const int &quaindex, double * const &basis) const;

    // Get grad R: dR_dx, dR_dy, dR_dz
    virtual void get_gradR(const int &quaindex, 
        double * const &basis_x, double * const &basis_y, 
        double * const &basis_z ) const;

    // Get d2R_dxx, d2R_dyy, d2R_dzz
    virtual void get_3D_R_gradR_LaplacianR(const int &quaindex,
        double * const &basis, double * const &basis_x, double * const &basis_y,
        double * const &basis_z, double * const &basis_xx, double * const &basis_yy,
        double * const &basis_zz ) const;

    // Get detJac
    virtual double get_detJac(const int &quaindex) const
    {return detJac[quaindex];}

    // Get Jacobian matrix
    virtual void get_Jacobian( const int &quaindex, double * const &jac_value ) const;

    // Get invJac
    virtual void get_invJacobian( const int &quaindex, double * const &dxi_dx ) const;

    // get unit normal vector together with surface element
    virtual void get_3d_normal_bottom( const int &quaindex, double &nx, double &ny, 
        double &nz, double &surf_area) const;
    
    virtual void get_3d_normal_top( const int &quaindex, double &nx, double &ny, 
        double &nz, double &surf_area) const;
    
    virtual void get_3d_normal_left( const int &quaindex, double &nx, double &ny, 
        double &nz, double &surf_area) const;
    
    virtual void get_3d_normal_right( const int &quaindex, double &nx, double &ny, 
        double &nz, double &surf_area) const;
    
    virtual void get_3d_normal_front( const int &quaindex, double &nx, double &ny, 
        double &nz, double &surf_area) const;
    
    virtual void get_3d_normal_back( const int &quaindex, double &nx, double &ny, 
        double &nz, double &surf_area) const;

    // Obtain element information
    virtual void print() const;

    virtual double get_memory_usage() const;

  private:
    int elem_index;
    int numQuapts;
    int nLocBas;
    bool is_sNonzero;

    double *R;
    double *dR_dx, *dR_dy, *dR_dz;
    double *d2R_dxx, *d2R_dyy,*d2R_dzz;
    double *detJac;

    // The dxi_dx value at quadrature points
    // Note: xi in [0,1].
    double * ds_dx_qua;

    // The array dx_ds_qua contains dx_ds, dx_dt, dx_du, dy_ds, ... , dz_du
    // at each quadrature points
    // Note: here in ds, s coordinate is in the knot domain. 
    double * dx_ds_qua; 

    // This is a private funtion that evaluates shape functions at
    // given quadrature point index.
    void BuildShape_atQua( 
        const int &ssdp1, const int &ttdp1, const int &uudp1,
        const int &quaindex,
        const int &quaindex_s, const int &quaindex_t, const int &quaindex_u,
        const double &hx, const double &hy, const double &hz,
        const double &invhx, const double &invhy, const double &invhz,
        const double * const &ctrl_x, 
        const double * const &ctrl_y, 
        const double * const &ctrl_z, 
        const double * const &ctrl_w,
        double * &N, double * &dN_ds, double * &dN_dt, double * &dN_du,
        double * &d2N_dss, double * &d2N_dtt, double * &d2N_duu,
        double * &d2N_dst, double * &d2N_dsu, double * &d2N_dtu,
        double * &dR_ds, double * &dR_dt, double * &dR_du,
        double * &d2R_dss, double * &d2R_dtt, double * &d2R_duu,
        double * &d2R_dst, double * &d2R_dsu, double * &d2R_dtu,
        double * &RHS, double * &sol,
        const double * const &Ns,
        const double * const &Nt,
        const double * const &Nu,
        const double * const &dNs_ds,
        const double * const &dNt_dt, 
        const double * const &dNu_du,
        const double * const &d2Ns_dss,
        const double * const &d2Nt_dtt,
        const double * const &d2Nu_duu );

    // ------------------------------------------------------------------------
    // normalize_3d_vector: given the component of an arbitrary 3d vector
    // this function normalize it to make sure the new vector has l2-norm 1.
    // ------------------------------------------------------------------------
    void normalize_3d_vector(double &x, double &y, double &z) const
    {
      const double temp = 1.0 / sqrt(x*x + y*y + z*z);
      x = x * temp;
      y = y * temp;
      z = z * temp;
    }

    
    void normalize_3d_vector(double &x, double &y, double &z, double &vec_length) const
    {
      vec_length = sqrt(x*x+y*y+z*z);
      const double temp = 1.0 / vec_length;
      x = x * temp;
      y = y * temp;
      z = z * temp;
    }

}; 

#endif
