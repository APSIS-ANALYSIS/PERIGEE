#ifndef PLOCASSEM_VMS_NS_HERK_HPP
#define PLOCASSEM_VMS_NS_HERK_HPP
// ==================================================================
// PLocAssem_Block_MS_NS_HERK.hpp
// 
// Parallel Local Assembly routine for VMS and HERK based NS
// solver into 5 sub blocks.
//
// Author: Yujie Sun
// Date: JMar. 5 2025
// ==================================================================
#include "ITimeMethod_RungeKutta.hpp"
#include "SymmTensor2_3D.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

class PLocAssem_Block_VMS_NS_HERK
{
  public:
    PetscScalar * Tangent0; // A00
    PetscScalar * Tangent1; // A01
    PetscScalar * Tangent2; // A10
    PetscScalar * Tangent3; // A11
    PetscScalar * Tangent4; // tilde{A}11
        
    PetscScalar * Residual0; // R0
    PetscScalar * Residual1; // R1

    // PetscScalar * sur_Residual0; // sur_R0
    PetscScalar * sur_Residual1; // sur_R1

    PLocAssem_Block_VMS_NS_HERK(
        const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s,
        const ITimeMethod_RungeKutta * const &tm_RK, const double &in_rho, 
        const double &in_vis_mu, const double &in_L0,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0, 
        const double &in_cu = 2.0, const double &in_cp = 2.0);

    ~PLocAssem_Block_VMS_NS_HERK();

    int get_dof() const {return 4;}

    int get_dof_mat() const {return 4;}

    int get_dof_mat_0() const {return 1;}
    
    int get_dof_mat_1() const {return 3;}

    int get_nLocBas() const {return nLocBas;}

    int get_snLocBas() const {return snLocBas;}

    void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size_p; ++ii) Residual0[ii] = 0.0;
      for(int ii=0; ii<vec_size_v; ++ii) Residual1[ii] = 0.0;
      for(int ii=0; ii<vec_size_p*vec_size_p; ++ii) Tangent0[ii] = 0.0;
      for(int ii=0; ii<vec_size_p*vec_size_v; ++ii) Tangent1[ii] = 0.0;
      for(int ii=0; ii<vec_size_v*vec_size_p; ++ii) Tangent2[ii] = 0.0;
      for(int ii=0; ii<vec_size_v*vec_size_v; ++ii) Tangent3[ii] = 0.0;
      for(int ii=0; ii<vec_size_v*vec_size_v; ++ii) Tangent4[ii] = 0.0;
    }

    void Zero_sur_Residual()
    {
      // for(int ii=0; ii<sur_size_v; ++ii) sur_Residual0[ii] = 0.0;
      for(int ii=0; ii<sur_size_v; ++ii) sur_Residual1[ii] = 0.0;
    }

    void Zero_Residual()
    {
      for(int ii=0; ii<vec_size_p; ++ii) Residual0[ii] = 0.0;
      for(int ii=0; ii<vec_size_v; ++ii) Residual1[ii] = 0.0;
    }

    void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size_p*vec_size_p; ++ii) Tangent0[ii] = 1.0;
      for(int ii=0; ii<vec_size_p*vec_size_v; ++ii) Tangent1[ii] = 1.0;
      for(int ii=0; ii<vec_size_v*vec_size_p; ++ii) Tangent2[ii] = 1.0;
      for(int ii=0; ii<vec_size_v*vec_size_v; ++ii) Tangent3[ii] = 1.0;
      for(int ii=0; ii<vec_size_v*vec_size_v; ++ii) Tangent4[ii] = 1.0;
    }

    void Assem_Residual_EBC_HERK_Sub(
        const int &ebc_id,
        const double &time, const double &dt,
        const int &subindex,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Residual_EBC_HERK_Final(
        const int &ebc_id,
        const double &time, const double &dt,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Residual_EBC_HERK_Pressure(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Tangent_Matrix(
      const double &time, const double &dt,
      const ITimeMethod_RungeKutta * const &tm_RK_ptr,
      const double * const &eleCtrlPts_x,
      const double * const &eleCtrlPts_y,
      const double * const &eleCtrlPts_z );

    void Assem_Residual_Sub(
        const double &time, const double &dt,
        const int &subindex,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const std::vector<std::vector<double>>& cur_velo_sols,
        const std::vector<std::vector<double>>& cur_pres_sols,
        const std::vector<std::vector<double>>& pre_velo_sols,
        const std::vector<std::vector<double>>& pre_pres_sols,
        const std::vector<double>& pre_velo,
        const std::vector<double>& pre_velo_before,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Residual_Final(
        const double &time, const double &dt,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const std::vector<std::vector<double>>& cur_velo_sols,
        const std::vector<double>& cur_velo,
        const std::vector<std::vector<double>>& cur_pres_sols,
        const std::vector<std::vector<double>>& pre_velo_sols,
        const std::vector<double>& pre_velo,
        const std::vector<std::vector<double>>& pre_pres_sols,
        const std::vector<double>& pre_velo_before,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Residual_Pressure(
        const double &time, const double &dt,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const std::vector<double>& cur_dot_velo,
        const std::vector<std::vector<double>>& cur_velo_sols,
        const std::vector<double>& cur_velo,
        const std::vector<std::vector<double>>& cur_pres_sols,
        const std::vector<double>& pre_velo,
        const std::vector<double>& cur_pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Tangent_Residual_Sub(
        const double &time, const double &dt,
        const int &subindex,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const std::vector<std::vector<double>>& cur_velo_sols,
        const std::vector<std::vector<double>>& cur_pres_sols,
        const std::vector<std::vector<double>>& pre_velo_sols,
        const std::vector<std::vector<double>>& pre_pres_sols,
        const std::vector<double>& pre_velo,
        const std::vector<double>& pre_velo_before,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Tangent_Residual_Final(
        const double &time, const double &dt,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const std::vector<std::vector<double>>& cur_velo_sols,
        const std::vector<double>& cur_velo,
        const std::vector<std::vector<double>>& cur_pres_sols,
        const std::vector<std::vector<double>>& pre_velo_sols,
        const std::vector<double>& pre_velo,
        const std::vector<std::vector<double>>& pre_pres_sols,
        const std::vector<double>& pre_velo_before,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    void Assem_Tangent_Residual_Pressure(
        const double &time, const double &dt,
        const ITimeMethod_RungeKutta * const &tm_RK_ptr,
        const std::vector<double>& cur_dot_velo,
        const std::vector<std::vector<double>>& cur_velo_sols,
        const std::vector<double>& cur_velo,
        const std::vector<std::vector<double>>& cur_pres_sols,
        const std::vector<double>& pre_velo,
        const std::vector<double>& cur_pres,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

  protected:
    // Private data
    const FEType elemType;
    
    const int nqpv, nqps;

    const std::unique_ptr<FEAElement> elementv, elements;

    const std::unique_ptr<IQuadPts> quadv, quads;

    const double rho0, vis_mu;
    
    const double CI, CT; // Constants for stabilization parameters
    
    const double Ctauc; // Constant scaling factor for tau_C
    
    const double L0, cu, cp; // Stabilization parameters for Darcy problem

    const int nLocBas, snLocBas, vec_size_v, vec_size_p, sur_size_v;

    // M matrix for tau_m
    //             mm[0], mm[1], mm[2]
    // M = coef *  mm[3], mm[4], mm[5]
    //             mm[6], mm[7], mm[8]
    const double coef;
    const std::array<double, 9> mm; 

    // Private functions
    void print_info() const;

    SymmTensor2_3D get_metric( const std::array<double, 9> &dxi_dx ) const;

    // Return tau_m and tau_c in RB-VMS
    std::array<double, 2> get_tau( const double &dt, 
        const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    std::array<double, 2> get_tau_Darcy( const double &dt, 
        const std::array<double, 9> &dxi_dx ) const;

    // Return tau_m_dot and tau_c_dot in RB-VMS
    std::array<double, 2> get_tau_dot( const double &dt, 
        const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return tau_bar := (v' G v')^-0.5 x rho0, 
    //        which scales like Time x Density
    // Users can refer to Int. J. Numer. Meth. Fluids 2001; 35: 93–116 
    // for more details
    double get_DC( const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    Vector_3 get_f(const Vector_3 &pt, const double &tt) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    Vector_3 get_H1(const Vector_3 &pt, const double &tt, 
        const Vector_3 &n_out ) const
    {
      const double p0 = 0.0;
      return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
    }

    typedef Vector_3 ( PLocAssem_Block_VMS_NS_HERK::*locassem_vms_ns_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_vms_ns_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, 
        const double &tt, const Vector_3 &n_out ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }
};

#endif
