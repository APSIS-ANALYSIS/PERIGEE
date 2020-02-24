#ifndef PLOCASSEM_TET_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_TET_VMS_NS_GENALPHA_HPP
// ==================================================================
// PLocAssem_Tet_VMS_NS_GenAlpha.hpp
// 
// Parallel Local Assembly routine for VMS and Gen-alpha based NS
// solver.
//
// Date: Feb. 10 2020
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_Tet_VMS_NS_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Tet_VMS_NS_GenAlpha(
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const int &elemtype = 501 );

    virtual ~PLocAssem_Tet_VMS_NS_GenAlpha();

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
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual double get_flowrate( const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void get_pressure_area( const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        double &pres, double &area );

    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id,
        const double &val,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_BackFlowStab(
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  private:
      // Private data
      const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta;

      const int nqp;

      double CI, CT;
      
      int nLocBas, snLocBas, vec_size, sur_size;

      std::vector<double> R, dR_dx, dR_dy, dR_dz;
      
      std::vector<double> d2R_dxx, d2R_dyy, d2R_dzz;

      double dxi_dx[9];

      std::vector< std::vector<double> > Sub_Tan, Sub_sur_Tan;

      // Private functions
      void print_info() const;

      void get_metric( const double * const &dxi_dx,
          double &G11, double &G12, double &G13,
          double &G22, double &G23, double &G33 ) const;

      // Return tau_m and tau_c in RB-VMS
      void get_tau( double &tau_m_qua, double &tau_c_qua,
          const double &dt, const double * const &dxi_dx,
          const double &u, const double &v, const double &w ) const;

      // Return tau_bar := (v' G v')^-0.5, which scales like a time variable
      void get_DC( double &dc_tau, const double * const &dxi_dx,
          const double &u, const double &v, const double &w ) const;

      void get_f(const double &x, const double &y, const double &z,
          const double &t, double &fx, double &fy, double &fz ) const
      {
        fx = 0.0; fy = 0.0; fz = 0.0;
      }

      void get_H1(const double &x, const double &y, const double &z,
          const double &t, const double &nx, const double &ny,
          const double &nz, double &gx, double &gy, double &gz ) const
      {
        const double p0 = 0.0;
        gx = p0*nx; gy = p0*ny; gz = p0*nz;
      }

      typedef void ( PLocAssem_Tet_VMS_NS_GenAlpha::*locassem_tet_vms_ns_funs )( const double &x, const double &y, const double &z,
          const double &t, const double &nx, const double &ny,
          const double &nz, double &gx, double &gy, double &gz ) const;

      locassem_tet_vms_ns_funs * flist;

      void get_ebc_fun( const int &ebc_id,
          const double &x, const double &y, const double &z,
          const double &t, const double &nx, const double &ny,
          const double &nz, double &gx, double &gy, double &gz ) const
      {
        return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz,gx,gy,gz);
      }

      void Zero_Sub_Tan()
      {
        for(int ii=0; ii<16; ++ii)
        {
          for(int jj=0; jj<nLocBas * nLocBas; ++jj) Sub_Tan[ii][jj] = 0.0;
        }
      }

      void Zero_Sub_sur_Tan()
      {
        for(int ii=0; ii<16; ++ii)
        {
          for(int jj=0; jj<snLocBas * snLocBas; ++jj) Sub_sur_Tan[ii][jj] = 0.0;
        }
      }
};

#endif
