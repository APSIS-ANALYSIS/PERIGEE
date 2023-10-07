#ifndef PLOCASSEM_LINEARPDE_GENALPHA_HPP
#define PLOCASSEM_LINEARPDE_GENALPHA_HPP
// ============================================================================
// PLocAssem_LinearPDE_GenAlpha.hpp
//
// Date: Oct 6 2023
// ============================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "FEAElement.hpp"
#include "ALocal_IEN.hpp"

class PLocAssem_LinearPDE_GenAlpha : public IPLocAssem_Linear
{
  public:
    PLocAssem_LinearPDE_GenAlpha( 
        const double &in_Young_modulus, const double &in_Poisson_ratio,
        const double &in_rho,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_snlocbas,
        const int &in_num_ebc_fun, const int &in_dof,
        const int &in_dof_mat, const int &elemtype = 501 );

    ~PLocAssem_LinearPDE_GenAlpha();

    int get_dof() const {return dof;}

    int get_dof_mat() const {return dof_mat;}

    void Zero_Mass_Stiffness_Load()
    {
      for(int ii=0; ii<vec_size; ++ii) Load[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) 
      {
        Mass[ii] = 0.0;
        Stiffness[ii] = 0.0;
      }
    }

    void Zero_Stiffness()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) Stiffness[ii] = 0.0;
    }

    void Zero_Mass()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) Mass[ii] = 0.0;
    }

    void Zero_Load()
    {
      for(int ii=0; ii<vec_size; ++ii) Load[ii] = 0.0;
    }

    void Zero_sur_Mass_Load()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Load[ii] = 0.0;
      for(int ii=0; ii<sur_size*sur_size; ++ii) sur_Mass[ii] = 0.0;
    }

    void Zero_sur_Mass()
    {
      for(int ii=0; ii<sur_size*sur_size; ++ii) sur_Mass[ii] = 0.0;
    }

    void Zero_sur_Load()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Load[ii] = 0.0;
    }

    void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) 
      {
        Mass[ii] = 1.0
        Stiffness[ii] = 1.0;
      }
    }

    void Assem_Load(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    void Assem_Stiffness(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    void Assem_Mass(
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    void Assem_Load_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );
    
    void Assem_Mass_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  private:
    // Private data
    const double Young_modulus, Possion_ratio, rho;
    const double alpha_f, alpha_m, gamma;
    
    const int num_ebc_fun;
    const int dof, dof_mat;

    int nLocBas, snLocBas, vec_size, sur_size;
    double mu, lambda;

    PetscScalar * Mass;
    PetscScalar * Stiffness;
    PetscScalar * Load;

    PetscScalar * sur_Mass;
    PetscScalar * sur_Load;

    void print_info() const;

    Vector_3 get_f( const Vector_3 &pt, const double &tt ) const
    {
      return Vector_3(0.0, 0.0, 0.0); 
    }

    typedef Vector_3 ( PLocAssem_LinearPDE_GenAlpha::*locassem_load_funs )( const Vector_3 &pt, const Vector_3 &nout, const double &t ) const;

    locassem_load_funs * flist;

    typedef double (PLocAssem_LinearPDE_GenAlpha::*locassem_robin_coefficient)( const Vector_3 &pt, const double &time) const;
    
    locassem_robin_coefficient * colist;

    Vector_3 get_ebc_fun( const int &ebc_id, const Vector_3 &pt, const Vector_3 &nout, const double &tt ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, nout, tt);
    }

    double get_robin_coefficient( const int &ebc_id, const Vector_3 &pt, const double &tt) const
    {
      return ((*this).*(colist[ebc_id]))(pt, tt);
    }

    double get_g_0( const Vector_3 &pt, const double &time ) const
    {
      return 0.0;
    }

    double get_g_1( const Vector_3 &pt, const double &time ) const
    {
      return 0.0;
    }
};

#endif
