#ifndef PLOCASSEM_VMS_NS_GENALPHA_HPP
#define PLOCASSEM_VMS_NS_GENALPHA_HPP
// ==================================================================
// PLocAssem_VMS_NS_GenAlpha.hpp
// 
// Parallel Local Assembly routine for VMS and Gen-alpha based NS
// solver.
//
// Author: Ju Liu
// Date: Feb. 10 2020
// ==================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "SymmTensor2_3D.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"

class PLocAssem_VMS_NS_GenAlpha : public IPLocAssem
{
  public:
    // PLocAssem_VMS_NS_GenAlpha(
    //     const TimeMethod_GenAlpha * const &tm_gAlpha,
    //     const int &in_nlocbas, const int &in_nqp,
    //     const int &in_snlocbas, const double &in_rho, 
    //     const double &in_vis_mu, const double &in_beta,
    //     const FEType &elemtype,
    //     const double &in_ct = 4.0, const double &in_ctauc = 1.0 );
    PLocAssem_VMS_NS_GenAlpha(
        const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s,
        const TimeMethod_GenAlpha * const &tm_gAlpha, const double &in_rho,
        const double &in_vis_mu, const double &in_beta,
        const double &in_ct = 4.0, const double &in_ctauc = 1.0 );

    virtual ~PLocAssem_VMS_NS_GenAlpha();

    virtual int get_dof() const {return 4;}

    virtual int get_dof_mat() const {return 4;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual int get_snLocBas() const {return snLocBas;}

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

    // virtual void Assem_Residual(
    //     const double &time, const double &dt,
    //     const double * const &dot_sol,
    //     const double * const &sol,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad );
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // virtual void Assem_Tangent_Residual(
    //     const double &time, const double &dt,
    //     const double * const &dot_sol,
    //     const double * const &sol,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad );
    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // virtual void Assem_Mass_Residual(
    //     const double * const &sol,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad );
    virtual void Assem_Mass_Residual(
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // virtual void Assem_Residual_EBC(
    //     const int &ebc_id,
    //     const double &time, const double &dt,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad );
    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // virtual double get_flowrate( const double * const &sol,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad );
    virtual double get_flowrate( const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z )

    // virtual void get_pressure_area( const double * const &sol,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad,
    //     double &pres, double &area );
    virtual void get_pressure_area( const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        double &pres, double &area );

    // virtual void Assem_Residual_EBC_Resistance(
    //     const int &ebc_id, const double &val,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad );
    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id, const double &val,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // virtual void Assem_Residual_BackFlowStab(
    //     const double * const &sol,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad );
    virtual void Assem_Residual_BackFlowStab(
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

    // virtual void Assem_Tangent_Residual_BackFlowStab(
    //     const double &dt,
    //     const double * const &sol,
    //     FEAElement * const &element,
    //     const double * const &eleCtrlPts_x,
    //     const double * const &eleCtrlPts_y,
    //     const double * const &eleCtrlPts_z,
    //     const IQuadPts * const &quad );
    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &sol,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z );

  protected:
    // Private data
    const FEType elemType;
    
    const int nqpv, npqs;

    const std::unique_ptr<FEAElement> elementv, elements, elementvs;

    const std::unique_ptr<IQuadPts> quadv, quads;

    const double rho0, vis_mu, alpha_f, alpha_m, gamma, beta;

    const double CI, CT; // Constants for stabilization parameters
    
    const double Ctauc; // Constant scaling factor for tau_C

    // const int nqp; // number of quadrature points
    
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
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    Vector_3 get_H1(const Vector_3 &pt, const double &tt, 
        const Vector_3 &n_out ) const
    {
      const double p0 = 0.0;
      return Vector_3( p0*n_out.x(), p0*n_out.y(), p0*n_out.z() );
    }

    typedef Vector_3 ( PLocAssem_VMS_NS_GenAlpha::*locassem_vms_ns_funs )( 
        const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const;

    locassem_vms_ns_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, 
        const double &tt, const Vector_3 &n_out ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt, n_out);
    }
};

#endif
