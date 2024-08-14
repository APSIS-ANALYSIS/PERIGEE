#ifndef PLOCASSEM_STRESS_RECOVERY_HPP
#define PLOCASSEM_STRESS_RECOVERY_HPP
// ============================================================================
// PLocAssem_Stress_Recovery.hpp
//
// Parallel Local Assembly routine for solution smoother
//
// Reference: Oden and Brauchli, IJNME, Vol.3, 317-325 (1971) 
// Date: Jan. 16 2024
// ============================================================================
#include "IPLocAssem.hpp"

class PLocAssem_Stress_Recovery : public IPLocAssem
{
  public: 
    PLocAssem_Stress_Recovery( const int &in_nlocbas,
       const double &in_nu, const double &in_E );
    
    virtual ~PLocAssem_Stress_Recovery();

    int get_dof() const { return 6; }

    virtual void Zero_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
    }

    virtual void Zero_Tangent_Residual()
    {
      for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 1.0;
    }

    virtual void Assem_Residual(
        const double * const &isol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );
    
    virtual void Assem_Mass_Residual(
        const double * const &sol_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );
    
  private:
    const int nLocBas, vec_size;
    const double modulus_E, nu, lambda, mu, l2mu;

    SymmTensor2_3D get_Cauchy_stress( const Tensor2_3D &F ) const;

    void print_info() const;
};

#endif
