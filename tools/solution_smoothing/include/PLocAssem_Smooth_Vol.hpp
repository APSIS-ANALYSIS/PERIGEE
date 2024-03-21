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
#include "IMaterialModel.hpp"

class PLocAssem_Smooth_Vol : public IPLocAssem
{
  public: 
    PLocAssem_Smooth_Vol(
        IMaterialModel * const &in_matmodel, const int &in_nlocbas );
    
    virtual ~PLocAssem_Smooth_Vol();

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
    IMaterialModel * matmodel;

    void print_info() const;
};

#endif