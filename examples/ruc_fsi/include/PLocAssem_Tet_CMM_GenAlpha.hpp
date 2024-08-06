#ifndef PLOCASSEM_TET_CMM_GENALPHA_HPP
#define PLOCASSEM_TET_CMM_GENALPHA_HPP
// ==================================================================
// PLocAssem_Tet_CMM_GenAlpha.hpp
// 
// Parallel Local Assembly routine for VMS and Gen-alpha based
// solver for the CMM type FSI problem.
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "SymmTensor2_3D.hpp"

class PLocAssem_Tet_CMM_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Tet_CMM_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nqp, const int &in_face_nqp, 
        const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const double &in_wall_rho, const double &in_nu,
        const double &in_kappa, const double &in_ctauc = 1.0,
        const int &elemtype = 501 );

    virtual ~PLocAssem_Tet_CMM_GenAlpha();

    virtual int get_dof() const {return 4;}

    virtual int get_dof_mat() const {return 4;}

    virtual double get_model_para_1() const {return alpha_f;}

    virtual double get_model_para_2() const {return gamma;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
      for(int ii=0; ii<sur_size*sur_size; ++ii) sur_Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Zero_sur_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 1.0;
    }

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Mass_Residual(
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual double get_flowrate( const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void get_pressure_area( const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        double &pres, double &area );

    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id, const double &val,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_BackFlowStab(
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC_Wall(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &sol_wall_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &ele_thickness,
        const double * const &ele_youngsmod,
        const double * const &ele_springconst,
        const double * const &ele_dampingconst,
        const double * const &qua_prestress,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual_EBC_Wall(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &sol_wall_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &ele_thickness,
        const double * const &ele_youngsmod,
        const double * const &ele_springconst,
        const double * const &ele_dampingconst,
        const double * const &qua_prestress,
        const IQuadPts * const &quad );

    virtual void get_Wall_CauchyStress(
        const double * const &sol_wall_disp,
        const FEAElement * const &element,
        const double * const &ele_youngsmod,
        std::vector<Tensor2_3D> &stress );

  private:
    // Private data
    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta;

    // wall properties: density, Poisson ratio, shear correction factor (kappa)
    const double rho_w, nu_w, kappa_w; 

    const int nqp, face_nqp; // number of quadrature points

    double CI, CT; // Constants for stabilization parameters

    const double Ctauc; // Constant scaling factor for tau_C

    int nLocBas, snLocBas, vec_size, sur_size;

    // Private functions
    void print_info() const;

    SymmTensor2_3D get_metric( const std::array<double, 9> &dxi_dx ) const;

    // Return tau_m and tau_c in RB-VMS
    std::array<double, 2> get_tau( 
        const double &dt, const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return tau_bar := (v' G v')^-0.5 x rho0, 
    //        which scales like Time x Density
    // Users can refer to Int. J. Numer. Meth. Fluids 2001; 35: 93–116 for more details
    double get_DC( const std::array<double, 9> &dxi_dx,
        const double &u, const double &v, const double &w ) const;

    // Return body force acting on the fluid domain
    Vector_3 get_f( const Vector_3 &pt, const double &tt ) const
    {
      return Vector_3(0.0, 0.0, 0.0);
    }

    // Return body force acting on the wall domain
    Vector_3 get_fw( const Vector_3 &pt, const double &tt ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    Vector_3 get_H1( const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const
    {
      const double p0 = 0.0;
      return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
    }

    typedef Vector_3 ( PLocAssem_Tet_CMM_GenAlpha::*locassem_tet_cmm_funs )(
      const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_tet_cmm_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id,
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }
};

#endif
