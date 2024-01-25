#ifndef PLOCASSEM_SMOOTH_VOL_HPP
#define PLOCASSEM_SMOOTH_VOL_HPP
// ============================================================================
// PLocAssem_Smooth_Vol.hpp
//
// Parallel Local Assembly routine for solution smoother
//
// Date: Jan. 16 2024
// ============================================================================
#include "IPLocAssem.hpp"

class PLocAssem_Smooth_Vol : public IPLocAssem
{
  public: 
    PLocAssem_Smooth_Vol(
        const double &in_module_E, const double &in_nu,
        const int &in_isol_dof, const int &in_osol_dof,
        const int &in_nlocbas );
    
    virtual ~PLocAssem_Smooth_Vol();

    int get_dof() const {return isol_dof;}

    int get_dof_mat() const {return osol_dof;}

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
    const double lambda, mu;
    const int isol_dof, osol_dof;
    const int nLocBas;
    const int vec_size;

    void print_info() const;
};

#endif