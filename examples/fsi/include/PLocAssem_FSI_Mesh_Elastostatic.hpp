#ifndef PLOCASSEM_FSI_MESH_ELASTOSTATIC_HPP
#define PLOCASSEM_FSI_MESH_ELASTOSTATIC_HPP
// ==================================================================
// PLocAssem_FSI_Mesh_Elastostatic.hpp
//
// This is the local assembly routine for the mesh motion in FSI 
// using the elastostatic equations to update the fluid domain.
//
// This assembly assumes:
// 1. The elastic modulus will be divided by the element jacobian
//    to improve the mesh quality;
//
// 2. There is NO element bc for the mesh equations.
//
// 3. There is NO body force in the right-hand side of the equations.
//
// 4. We only implement the residual-tangent assembly routine, wherein
//    vec_a gives the previous disp,
//    vec_b gives the current disp.
//
// Date Created: Oct. 15 2017
// Author: Ju Liu
// ==================================================================
#include "IPLocAssem.hpp"

class PLocAssem_FSI_Mesh_Elastostatic : public IPLocAssem
{
  public:
    PLocAssem_FSI_Mesh_Elastostatic( const double &in_mat_E, 
        const double &in_mat_nu, const int &in_nlocbas );

    virtual ~PLocAssem_FSI_Mesh_Elastostatic();

    virtual int get_dof() const {return 3;}

    virtual int get_dof_mat() const {return 3;}

    virtual int get_num_ebc_fun() const {return 0;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual int get_snLocBas() const {return 3;}

    virtual void Zero_Tangent_Residual();

    virtual void Zero_Residual();

    virtual void Assem_Estimate();

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::print_fatal("Error: PLocAssem_FSI_Mesh_Laplacian::Assem_Residual is not implemented. \n");}

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Mass_Residual(
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::print_fatal("Error: PLocAssem_FSI_Mesh_Laplacian::Assem_Mass_Residual is not implemented. \n");}

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::print_fatal("Error: PLocAssem_FSI_Mesh_Laplacian::Assem_Residual_EBC is not implemented. \n");}

  private:
    const double E, nu, lambda, mu, kappa;
    const int nLocBas, vec_size;

    void print_info() const;

    void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const double * const &sol,
        double * const &cur_x,
        double * const &cur_y,
        double * const &cur_z );
};

#endif
