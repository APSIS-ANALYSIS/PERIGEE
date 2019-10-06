#ifndef FEAELEMENT_NURBS_3D_DER2_VMS_HPP
#define FEAELEMENT_NURBS_3D_DER2_VMS_HPP
// ==================================================================
// FEAElement_NURBS_3D_der2_vms.hpp
// It is a finite element implementation of the physical 3D NURBS 
// element, with derivative up to 2nd order.
// 
// _vms: version vms, this class can be viewed as a public inheritance
//       of the FEAElement_NURBS_3D_der2.hpp. In addition to the
//       classical functions, I have the dxi/dx values at all quadrature
//       points calculated and cached in this class. 
//       The dxi/dx's are stored in ds_dx_qua, ds_dy_qua, ..., du_dz_qua 
//       vectors with length numQuapts. These values will be used to
//       calculate the G_ij, and g_i, which will further be used to 
//       calculate tau_m, tau_c, the stabilization parameters. The 
//       definitions of G_ij and g_i can be found on CMAME 197 (2007),
//       pp 173-201. 
//
// Date:
// Feb. 11th 2015
// ==================================================================
#include <vector>
#include "Sys_Tools.hpp"
#include "Matrix_double_6by6_Array.hpp"
#include "ALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "FEANode.hpp"
#include "IAExtractor.hpp"
#include "BernsteinBasis_Array.hpp"
#include "FEAElement.hpp"


class FEAElement_NURBS_3D_der2_vms : public FEAElement
{
  public:
    FEAElement_NURBS_3D_der2_vms( const int &in_eIndex,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const FEANode * const &feaNode,
        const IAExtractor * const &extractor,
        const ALocal_IEN * const &locIEN );

    virtual ~FEAElement_NURBS_3D_der2_vms();

    virtual int get_elemIndex() const {return elem_index;}
    virtual int get_elemDim() const {return 3;}

    // Type = 301: 3D NURBS H1 conforming basis with up to 
    //             2nd order derivatives, together with dxi/dx.
    virtual int get_Type() const {return 301;}

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

    // Get invJac
    virtual void get_invJacobian(const int &quaindex, double * const &dxi_dx) const;

    virtual void print() const;

    virtual double get_memory_usage() const;

  private:
    int elem_index;
    int numQuapts;
    int nLocBas;
    bool is_sNonzero;

    double *R;
    double *dR_dx, *dR_dy, *dR_dz;
    double *d2R_dxx, *d2R_dyy,*d2R_dzz, *d2R_dxy, *d2R_dxz, *d2R_dyz;
    double *detJac;

    // This is the dxi_dx that will be needed when calculating the stabilization
    // parameter in VMS. See the definition of tau_m & tau_c on pp 181, CMAME
    // 197 (2007)
    double *ds_dx_qua, *ds_dy_qua, *ds_dz_qua;
    double *dt_dx_qua, *dt_dy_qua, *dt_dz_qua;
    double *du_dx_qua, *du_dy_qua, *du_dz_qua;

    // This is a private funtion that evaluates shape functions at
    // given quadrature point index.
    void BuildShape_atQua( 
        const int &ssdp1, const int &ttdp1, const int &uudp1,
        const int &quaindex,
        const int &quaindex_s, const int &quaindex_t, const int &quaindex_u,
        const double &hx, const double &hy, const double &hz,
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
}; 
#endif
