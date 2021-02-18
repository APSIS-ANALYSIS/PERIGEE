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
        const ALocal_NodalBC * const &part_nbc,
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
        const ALocal_NodalBC * const &nbc_part,
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
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
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
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
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
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_Inflow_NodalBC * const &infbc_part );

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
        const ALocal_Inflow_NodalBC * const &infbc_part );

    // **** PRESTRESS TODO:
    // **** Loop over num_selem:
    //        - call PLocAssem_Tet_CMM_GenAlpha::get_Wall_CauchyStress()
    //        - call ALocal_Wall_Prestress::update_prestress()
    // virtual void Update_Wall_Prestress(
    //     const PDNSolution * const &sol_wall_disp,
    //     IPLocAssem * const &lassem_ptr,
    //     FEAElement * const &element_w,
    //     const IQuadPts * const &quad_s,
    //     const ALocal_EBC * const &ebc_wall_part,
    //     ALocal_Wall_Prestress * const &wall_prestress );

  private:
    // Private data
    const int nLocBas, dof_sol, dof_mat, num_ebc, nlgn;
    
    int snLocBas;

    // Private function
    // Essential boundary condition
    void EssBC_KG( const ALocal_NodalBC * const &nbc_part, const int &field );
    
    void EssBC_G( const ALocal_NodalBC * const &nbc_part, const int &field );

    // Natural boundary condition
    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    // Backflow integral on outlet surfaces
    void BackFlow_G( const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    void BackFlow_KG( const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    // Resistance type boundary condition on outlet surfaces
    void NatBC_Resis_G( const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    void NatBC_Resis_KG( const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    // Wall integral for thin-walled linear membrane
    // **** PRESTRESS TODO: additional arg ALocal_Wall_Prestress
    void WallMembrane_G( const double &curr_time,
        const double &dt, 
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol_wall_disp,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_w,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_wall_part );

    // **** PRESTRESS TODO: additional arg ALocal_Wall_Prestress
    void WallMembrane_KG( const double &curr_time,
        const double &dt, 
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol_wall_disp,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_w,
        const IQuadPts * const &quad_s,
        const ALocal_NodalBC * const &nbc_part,
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
