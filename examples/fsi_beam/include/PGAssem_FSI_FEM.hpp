#ifndef PGASSEM_FSI_FEM_HPP
#define PGASSEM_FSI_FEM_HPP
// ==================================================================
// PGAssem_FSI_FEM.hpp
//
// Parallel golbal assembly based on PETSc, using AIJ matrix format.
// The assembly routine is designed for classical C0 FEM method, which
// means we do not need extraction operators and local mesh sizes.
//
// The assembly is for the segregated solver for FSI problem. 
// The input solution vectors contains, in fluid,
//  [ mesh disp; pressure; velocity ],
// the dot solution contains
//  [ mesh velo; dot pressure; dot velcoty ];
// in solid,
//  [ displacement; pressure; velocity ],
// the dot solution contains
//  [ dot disp; dot pressure; dot velo ].
//
// Note: 1. We only assume that the EBC are applied for fluids. The
//          NatBC_G performs for fluids only for now.
//
// Author: Ju Liu
// Date Created: Aug 10 2017
// ==================================================================
#include "IPGAssem.hpp"

class PGAssem_FSI_FEM : public IPGAssem
{
  public:
    PGAssem_FSI_FEM( IPLocAssem const * const &locassem_f_ptr,
        IPLocAssem const * const &locassem_s_ptr,
        IAGlobal_Mesh_Info const * const &agmi_ptr,
        ALocal_Elem const * const &alelem_ptr,
        ALocal_IEN const * const &aien_ptr,
        APart_Node const * const &pnode_ptr,
        ALocal_NodalBC const * const &part_nbc,
        ALocal_EBC const * const &part_ebc );

    virtual ~PGAssem_FSI_FEM();

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part );

    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );


    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );


    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

  private:
    int nLocBas, snLocBas, dof_sol, dof_mat, num_ebc;

    PetscInt * row_index, * srow_index;

    double * array_a, * array_b; // length: dof_sol x nlocalghonode 

    double * local_a, * local_b; // length: dof_sol x nLocBas

    double * local_as, * local_bs; // length: dof_sol x snLocBas

    int * IEN_e, * LSIEN;

    double * ectrl_x, * ectrl_y, * ectrl_z;
    double * sctrl_x, * sctrl_y, * sctrl_z;

    void EssBC_KG( const ALocal_NodalBC * const &nbc_part, const int &field );
    void EssBC_G( const ALocal_NodalBC * const &nbc_part, const int &field );


    // For now, we only do ebc boundary integration for fluids 
    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem * const &lassem_f_ptr,
        FEAElement * const &element_s,
        const int &in_loc_dof,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const lien_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );


    void GetLocal(const double * const &array, const int * const &IEN,
        double * const &local_array) const
    {
      int offset1, offset2;
      for(int ii=0; ii<nLocBas; ++ii)
      {
        offset1 = ii * dof_sol;
        offset2 = IEN[ii] * dof_sol;
        for(int jj=0; jj<dof_sol; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

    void GetLocal( const double * const &array, const int * const &IEN,
        const int &in_locbas, double * const &local_array) const
    {
      int offset1, offset2;
      for(int ii=0; ii<in_locbas; ++ii)
      {
        offset1 = ii * dof_sol;
        offset2 = IEN[ii] * dof_sol;
        for(int jj=0; jj<dof_sol; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

    void Get_dnz_onz( const int &nlocnode,
        const int &empirical_neighbor_node_number,
        const ALocal_NodalBC * const &nbc_ptr,
        PetscInt * const &dnz, PetscInt * const &onz ) const;

};

#endif
