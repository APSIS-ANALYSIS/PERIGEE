#ifndef PLOCASSEM_2X2BLOCK_TET4_VMS_INCOMPRESSIBLE_HPP
#define PLOCASSEM_2X2BLOCK_TET4_VMS_INCOMPRESSIBLE_HPP
// ============================================================================
// PLocAssem_2x2Block_Tet4_VMS_Incompressible.hpp
//
// This is a local assembly for VMS formulation of the fully incompressible
// hyperelasticity.
//
// Author: Ju Liu
// Date: Jan 1 2022
// ============================================================================
#include "IPLocAssem_2x2Block.hpp"
#include "IMaterialModel.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_2x2Block_Tet4_VMS_Incompressible : public IPLocAssem_2x2Block
{
  public:
    PLocAssem_2x2Block_Tet4_VMS_Incompressible( 
        IMaterialModel * const &in_matmodel,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nqp );

    virtual ~PLocAssem_2x2Block_Tet4_VMS_Incompressible();
    
  private:
    const double rho0, alpha_f, alpha_m, gamma;

    // memory layout
    // dof_per_node = 7 to make it compatible with the problem setting
    // vec_size = 4 * nLocBas, which defines the local matrix/vector length
    const int nLocBas, vec_size, nqp, snLocBas;

    // useful tensors for the material model
    IMaterialModel * matmodel;

    void print_info() const;

    void get_tau( double &tau_m_qua, double &tau_c_qua,
        const double &dt, const double &Jin, const double &dx ) const;

    Vector_3 get_f(const double &x, const double &y, const double &z,
        const double &t ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    // Use pointers to the member functions to facilitate the automatic
    // treatment of ebc surface integration.
    typedef Vector_3 ( PLocAssem_2x2Block_Tet4_VMS_Incompressible::*locassem_2x2block_vms_ela_fem_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz ) const;

    locassem_2x2block_vms_ela_fem_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz);
    }

};

#endif
