#ifndef PLOCASSEM_VMS_NS_SEMIBDF1_HPP
#define PLOCASSEM_VMS_NS_SEMIBDF1_HPP
// ==================================================================
// PLocAssem_VMS_NS_GenAlpha.hpp
// 
// Parallel Local Assembly routine for VMS and Semi-BDF1 based NS
// solver.
//
// Author : Chi Ding
// Date   : March 4, 2024
// ==================================================================
#include "IPLocAssem.hpp"
#include "SymmTensor2_3D.hpp"
#include "Math_Tools.hpp"

class PLocAssem_VMS_NS_SemiBDF1 : public IPLocAssem
{
  public:
    PLocAssem_VMS_NS_SemiBDF1(
        const int &in_nlocbas, const int &in_nqp,
        const int &in_snlocbas, const double &in_rho, 
        const double &in_vis_mu, const double &in_beta,
        const int &elemtype,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0 );

    virtual ~PLocAssem_VMS_NS_SemiBDF1();

    virtual int get_dof() const {return 4;}

    virtual int get_dof_mat() const {return 4;}

    //virtual double get_model_para_1() const {return alpha_f;}

    //virtual double get_model_para_2() const {return gamma;}

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
        const double * const &sol_0,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &sol_0,
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
        const double * const &sol_0,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  protected:
    // Private data
    const double rho0, vis_mu, beta;

    const double CI, CT; // Constants for stabilization parameters
    
    const double Ctauc; // Constant scaling factor for tau_C

    const int nqp; // number of quadrature points
    
    const int nLocBas, snLocBas, vec_size, sur_size;

    // M matrix for tau_m
    //             mm[0], mm[1], mm[2]
    // M = coef *  mm[3], mm[4], mm[5]
    //             mm[6], mm[7], mm[8]
    const double coef;
    const std::array<double, 9> mm; 

    // Private functions
    virtual void print_info() const;

    SymmTensor2_3D get_metric( const std::array<double, 9> &dxi_dx ) const;

    // Return tau_m and tau_c in RB-VMS
    std::array<double, 2> get_tau( const double &dt, 
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

      const double RR   = 1.0;
      const double RR_3 = RR * RR * RR;
      const double Area = MATH_T::PI * RR * RR;

      const double yy   = pt.y(), zz = pt.z();
      const double rr   = std::sqrt(yy * yy + zz * zz);
      const double rr_3 = rr * rr * rr;
     
      const double fr    = 5.0/(3.0 * Area) * (1 - rr_3/RR_3 );
      //const double d2f_dxx = 0.0;
      //const double d2f_dyy = 5.0/(3.0 * Area) * ( -3.0 * (yy * yy / rr + rr) / RR_3 );
      //const double d2f_dzz = 5.0/(3.0 * Area) * ( -3.0 * (zz * zz / rr + rr) / RR_3 );

      const double Q_thd = 1.0;
      const double T_thd = 1.0;
      
      double Qt = 0.0, dQ_dt = 0.0;

      // if-else statement
      if (tt < T_thd)
    {
        Qt    = Q_thd * 0.5 * (1-std::cos(MATH_T::PI*tt/T_thd));
        dQ_dt = Q_thd * 0.5 * (MATH_T::PI/T_thd * std::sin(MATH_T::PI*tt/T_thd));
    }
     else
    { 
        Qt    = Q_thd;
        dQ_dt = 0.0;
    }

      //const double d2u_dxx = Qt * d2f_dxx;
      //const double d2u_dyy = Qt * d2f_dyy;
      //const double d2u_dzz = Qt * d2f_dzz;

      const double du_dt  = dQ_dt * fr;

      const double fl_mu = 4.0e-2, fl_density = 1.0;

      const double f_acc = fl_density * du_dt ; 
      const double dp_dx = -100.0;
      const double f_vis = fl_mu * Qt * 5.0/(3.0 * Area) * ( -3.0 * (3.0 * rr) / RR_3 );

      double fx =  (f_acc + dp_dx - f_vis)/ fl_density;
      
      return Vector_3( fx, 0.0, 0.0 );
    }

    Vector_3 get_H1(const Vector_3 &pt, const double &tt, 
        const Vector_3 &n_out ) const
    {
      const double p0 = 0.0;
      return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
    }

    typedef Vector_3 ( PLocAssem_VMS_NS_SemiBDF1::*locassem_vms_ns_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_vms_ns_funs * flist;

    Vector_3 get_cubic_velo_traction(const Vector_3 &pt, const double &tt, const Vector_3 &n_out) const
    {
      const double RR   = 1.0;
      const double RR_3 = RR * RR * RR;
      const double Area = MATH_T::PI * RR * RR;

      const double yy   = pt.y(), zz = pt.z();
      const double rr   = std::sqrt(yy * yy + zz * zz);

      const double df_dx = 0.0;
      const double df_dy = 5.0/(3.0 * Area) * ( -3.0 * rr * yy / RR_3 );
      const double df_dz = 5.0/(3.0 * Area) * ( -3.0 * rr * zz / RR_3 );

      const double Q_thd = 1.0;
      const double T_thd = 1.0;
      
      // if-else statement
      double Qt = 0.0;
      if (tt < T_thd)  Qt = Q_thd * 0.5 * (1-std::cos(MATH_T::PI*tt/T_thd));
      else  Qt  = Q_thd;
      
      const double du_dx = Qt * df_dx; 
      const double du_dy = Qt * df_dy;
      const double du_dz = Qt * df_dz;

      const double fl_mu = 4.0e-2;
      
      Tensor2_3D velo_grad (du_dx, du_dy, du_dz,
                            0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0);

      Tensor2_3D velo_grad_T = velo_grad;
      velo_grad_T.transpose();

      Tensor2_3D SS = velo_grad + velo_grad_T;

      SS *= fl_mu;

      return SS * n_out;
    }

    Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, 
        const double &tt, const Vector_3 &n_out ) const
    {
      //return Vector_3( 0.0, 0.0, 0.0 ); 
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }
};

#endif