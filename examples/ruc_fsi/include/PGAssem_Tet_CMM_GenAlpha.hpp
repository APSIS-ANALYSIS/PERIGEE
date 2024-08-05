#ifndef PGASSEM_TET_CMM_GENALPHA_HPP
#define PGASSEM_TET_CMM_GENALPHA_HPP
// ==================================================================
// PGAssem_Tet_CMM_GenAlpha.hpp
//
// Parallel global assembly based on PETSc, using AIJ matrix format.
// The assembly routine is designed for classical C0 FEM method, which
// means we do not need extraction operators and local mesh sizes.
//
// The assembly is for the CMM-type FSI problem written in VMS formualtion
// for the NS equations. The input solution vectors contains
//  [ pressure; velocity ],
// the dot solution contains
//  [ dot pressure; dot velocity ].
//
// Author: Ju Liu 
// Date Created: Feb. 10 2020
// ==================================================================
#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"
#include "PDNSolution_NS.hpp"
#include "PDNSolution_Wall_Disp.hpp"

class PGAssem_Tet_CMM_GenAlpha : public IPGAssem
{
  public:
    // Constructor for CMM equations
    PGAssem_Tet_CMM_GenAlpha( 
        IPLocAssem * const &locassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quads,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const ALocal_NBC * const &part_nbc,
        const ALocal_RingBC * const &part_ringnbc,
        const ALocal_EBC * const &part_ebc,
        const IGenBC * const &gbc,
        const int &in_nz_estimate=60 );

    // Destructor
    virtual ~PGAssem_Tet_CMM_GenAlpha();

    // Nonzero pattern estimate for the CMM equations
    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    // Assemble mass matrix and residual vector
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part );

    // Assemble the residual vector for the CMM equations
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &sol_wall_disp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementw,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_wall_part,
        const IGenBC * const &gbc );

    // Assemble the residual vector and tangent matrix 
    // for the CMM equations
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &sol_wall_disp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementw,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_wall_part,
        const IGenBC * const &gbc );

    // Assembly routine for the surface integrals of flow rate and
    // pressure
    virtual double Assem_surface_flowrate(
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id );

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id );

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id );

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id );

    // Assembly routine for a matrix-free manner of getting the tangent
    // stiffness due to the reduced model coupling.
    // The effective m = dP/dQ will be calculated, and 
    // n = XX dot int_NA will
    // be calculated by performing a loop over surface elements. Then
    // YY will be int_NA scaled by alpha_f gamma dt m n.
    virtual void Assem_matrix_free_K( const Vec &XX,
        const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbci,
        Vec &YY );

  private:
    // Private data
    const int nLocBas, dof_sol, dof_mat, num_ebc, nlgn;

    int snLocBas;

    // Private function
    // Essential boundary condition
    void EssBC_KG( const ALocal_NBC * const &nbc_part, const int &field );

    void EssBC_G( const ALocal_NBC * const &nbc_part, const int &field );

    // Ring nodal BC: 1) clamped, or 2) in-plane motion (skew bc)
    // References:
    //   i. Griffiths DV (Computers & Structures 1990) Treatment of skew boundary
    //      conditions in finite element analysis
    //  ii. Bathe KJ (1996) Finite element procedures
    void RingBC_KG(
        const ALocal_RingBC * const &ringnbc_part,
        const int &dof, const int &nrow, const int &ncol,
        const PetscInt * const &row_index,
        const PetscInt * const &col_index,
        PetscScalar * const &Ke,
        PetscScalar * const &Ge );

    void RingBC_G(
        const ALocal_RingBC * const &ringnbc_part,
        const int &dof, const int &nrow,
        const PetscInt * const &row_index,
        PetscScalar * const &Ge );

    // Natural boundary condition
    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part );

    // Backflow integral on outlet surfaces
    void BackFlow_G( const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part );

    void BackFlow_KG( const double &dt,
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part );

    // Resistance type boundary condition on outlet surfaces
    void NatBC_Resis_G( const double &curr_time, const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    // Note: to be replaced by the SHELL approach
    void NatBC_Resis_KG( const double &curr_time, const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    // Wall integral for thin-walled linear membrane
    void WallMembrane_G( const double &curr_time,
        const double &dt, 
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &sol_wall_disp,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_w,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_wall_part );

    void WallMembrane_KG( const double &curr_time,
        const double &dt, 
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &sol_wall_disp,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_w,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_RingBC * const &ringnbc_part,
        const ALocal_EBC * const &ebc_wall_part );

    void GetLocal(const double * const &array, const int * const &IEN,
        double * const &local_array) const
    {
      for(int ii=0; ii<nLocBas; ++ii)
      {
        const int offset1 = ii * dof_sol;
        const int offset2 = IEN[ii] * dof_sol;
        for(int jj=0; jj<dof_sol; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

    void GetLocal( const double * const &array, const int * const &IEN,
        const int &in_locbas, double * const &local_array) const
    {
      for(int ii=0; ii<in_locbas; ++ii)
      {
        const int offset1 = ii * dof_sol;
        const int offset2 = IEN[ii] * dof_sol;
        for(int jj=0; jj<dof_sol; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

    void GetLocal( const double * const &array, const int * const &IEN,
        const int &in_locbas, const int &in_dof,
        double * const &local_array) const
    {
      for(int ii=0; ii<in_locbas; ++ii)
      {
        const int offset1 = ii * in_dof;
        const int offset2 = IEN[ii] * in_dof;
        for(int jj=0; jj<in_dof; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }
};

#endif
