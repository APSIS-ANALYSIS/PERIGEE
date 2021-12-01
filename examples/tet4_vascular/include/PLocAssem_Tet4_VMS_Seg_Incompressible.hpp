#ifndef PLOCASSEM_TET4_VMS_SEG_INCOMPRESSIBLE_HPP
#define PLOCASSEM_TET4_VMS_SEG_INCOMPRESSIBLE_HPP
// ==================================================================
// PLocAssem_Tet4_VMS_Seg_Incompressible.hpp
//
// This is the local assembly routine for the stabilized mixed 
// formulation for the p-v implicit solver using the segregated
// algorithm for FULLY INCOMPRESSIBLE materials.
//
// To simplify the implementation, this local assembly code assumes
// linear tetrahedral element, and this file is in the tet project
// folder. Second-order derivatives will be ignored.
//
// Time integration is the Generalized-alpha method.
//
// Author: Ju Liu
// Date Created: June 12 2017
// ==================================================================
#include "IPLocAssem.hpp"
#include "IMaterialModel.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_Tet4_VMS_Seg_Incompressible : public IPLocAssem
{
  public:
    PLocAssem_Tet4_VMS_Seg_Incompressible(
        IMaterialModel * const &in_matmodel,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nqp );

    virtual ~PLocAssem_Tet4_VMS_Seg_Incompressible();

    virtual int get_dof() const {return dof_per_node;}

    virtual int get_dof_mat() const {return 4;}

    virtual int get_num_ebc_fun() const {return num_ebc_fun;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Assem_Estimate();

    // Assembly routines.
    // Note: vec_a / vec_b has 3 + 1 + 3 = 7 d.o.f.
    //       vec_a : dot(variable), i.e., velo,
    //       vec_b : variable, i.e., disp.    
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
        const IQuadPts * const &quad );
    
    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
        const IQuadPts * const &quad );
    
    virtual void Assem_Mass_Residual(
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &qua_prestress,
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

    // ------------------------------------------------------------------------
    // This function will calculate the Cauchy stress at every quadrature points
    // within this element. The output stress has length quad -> get_num_quadPts()
    // ------------------------------------------------------------------------
    virtual void get_Wall_CauchyStress(
        const double * const &disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        std::vector<Matrix_3x3> &stress ) const;

  private:
    const double rho0, alpha_f, alpha_m, gamma;

    const int num_ebc_fun;

    // memory layout
    // dof_per_node = 7 to make it compatible with the problem setting
    // vec_size = 4 * nLocBas, which defines the local matrix/vector length
    const int nLocBas, dof_per_node, vec_size, nqp, snLocBas;

    // useful tensors for the material model
    IMaterialModel * matmodel;

    void print_info() const;

    void get_tau( double &tau_m_qua, double &tau_c_qua,
        const double &dt, const double &Jin,
        const double &dx ) const;

    Vector_3 get_f(const double &x, const double &y, const double &z,
        const double &t ) const
    {
      return Vector_3( 0.0, 0.0, 0.0 );
    }

    // Use pointers to the member functions to facilitate the automatic
    // treatment of ebc surface integration.
    typedef Vector_3 ( PLocAssem_Tet4_VMS_Seg_Incompressible::*locassem_vms_seg_ela_fem_funs )( const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz ) const;

    locassem_vms_seg_ela_fem_funs * flist;

    Vector_3 get_ebc_fun( const int &ebc_id,
        const double &x, const double &y, const double &z,
        const double &t, const double &nx, const double &ny,
        const double &nz ) const
    {
      return ((*this).*(flist[ebc_id]))(x,y,z,t,nx,ny,nz);
    }

};

#endif
