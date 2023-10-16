#ifndef PLOCASSEM_ELASTODYNAMICS_GENALPHA_HPP
#define PLOCASSEM_ELASTODYNAMICS_GENALPHA_HPP
// ============================================================================
// PLocAssem_Elastodynamics_GenAlpha.hpp
//
// Date: Jan 21 2022
// ============================================================================
#include "IPLocAssem.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_Elastodynamics_GenAlpha : public IPLocAssem
{
  public:
    PLocAssem_Elastodynamics_GenAlpha( 
        const double &in_rho, const double &in_module_E, const double &in_nu,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_snlocbas,
        const int &in_num_ebc_fun );

    virtual ~PLocAssem_Elastodynamics_GenAlpha();

    virtual int get_dof() const {return 6;}

    virtual int get_dof_mat() const {return 3;}

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Zero_sur_Tangent_Residual()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Residual[ii] = 0.0;
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
        const double * const &dot_sol_velo,
        const double * const &dot_sol_disp,
        const double * const &sol_velo,
        const double * const &sol_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_sol_velo,
        const double * const &dot_sol_disp,
        const double * const &sol_velo
        const double * const &sol_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Mass_Residual(
        const double * const &sol_velo,
        const double * const &sol_disp,
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

  private:
    // Private data
    const double rho, module_E, nu, lambda, mu;
    const double alpha_f, alpha_m, gamma;
    
    const int num_ebc_fun;

    const int nLocBas, snLocBas;
    const int vec_size, sur_size;

    void print_info() const;

    Vector_3 get_f( const Vector_3 &pt, const double &tt ) const
    {
      return Vector_3(0.0, 0.0, 0.0);
    }

    typedef Vector_3 ( PLocAssem_Elastodynamics_GenAlpha::*locassem_transport_funs )
        ( const Vector_3 &pt, const double &t, const Vector_3 &n_out ) const;

    locassem_transport_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, const double &tt, const Vector_3 &n_out ) const
    {
      return Vector_3(0.0, 0.0, 0.0);
    }

    Vector_3 get_g_0( const Vector_3 &pt, const double &time, const Vector_3 &n_out ) const
    {
      return Vector_3(0.0, 0.0, 0.0);
    }

    Vector_3 get_g_1( const Vector_3 &pt, const double &time, const Vector_3 &n_out ) const
    {
      return Vector_3(0.0, 0.0, 0.0);
    }
};

#endif
