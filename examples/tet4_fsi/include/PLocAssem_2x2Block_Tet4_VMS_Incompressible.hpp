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
    
    virtual int get_dof_0() const {return 3;}

    virtual int get_dof_1() const {return 1;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual(); 
    
    virtual void Assem_Estimate();

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
        const IQuadPts * const &quad );

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_disp,
        const double * const &dot_velo,
        const double * const &dot_pres,
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
        const IQuadPts * const &quad );

    virtual void Assem_Mass_Residual(
        const double * const &disp,
        const double * const &velo,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Residual_EBC(
        const double &time,
        const double * const &pres,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    // ------------------------------------------------------------------------
    // This function will calculate the Cauchy stress at every quadrature points
    // within this element. The output stress has length quad -> get_num_quadPts()
    // ------------------------------------------------------------------------
    virtual std::vector<Matrix_3x3> get_Wall_CauchyStress(
        const double * const &disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad ) const;

  private:
    const double rho0, alpha_f, alpha_m, gamma;

    // memory layout
    // vec_size = 4 * nLocBas, which defines the local matrix/vector length
    const int nLocBas, vec_size_0, vec_size_1, nqp, snLocBas;

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
