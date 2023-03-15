#ifndef PGASSEM_HEAT_HPP
#define PGASSEM_HEAT_HPP
// ==================================================================
// PGAssem_Heat.hpp
// Global assembly for the laplace linear operator
// ==================================================================
#include "IPGAssem.hpp"
#include "ALocal_Elem.hpp"
#include "IAGlobal_Mesh_Info.hpp"
#include "IPLocAssem.hpp"
#include "ALocal_NBC.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "ALocal_EBC.hpp"
#include "PETSc_Tools.hpp"
#include "PDNSolution.hpp"

class PGAssem_Heat : public IPGAssem
{
  public:
    Vec * intNA;
    Vec b, diag;
    KSP ksp_K;
    
    PGAssem_Heat(
        IPLocAssem * const &locassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quads,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const ALocal_NBC * const &part_nbc,
        const ALocal_EBC * const &part_ebc,
        const int &in_nz_estimate=60 );

    virtual ~PGAssem_Heat();

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual void Assem_intNA(
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    virtual int get_num_ebc() const {return num_ebc;}

  private:
    const int nLocBas, dof_sol, dof_mat, num_ebc;

    int snLocBas;

    PetscInt * row_index, * srow_index;

    double * array_a, * array_b; // length: dof_sol x nlocalghonode 

    double * local_a, * local_b; // length: dof_sol x nLocBas

    double * local_as, * local_bs; // length: dof_sol x snLocBas

    int * IEN_e, * LSIEN;

    double * ectrl_x, * ectrl_y, * ectrl_z;

    double * sctrl_x, * sctrl_y, * sctrl_z;

    void EssBC_KG( const ALocal_NBC * const &nbc_part, const int &field );

    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const int &in_loc_dof,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
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
        const ALocal_NBC * const &nbc_ptr,
        PetscInt * const &dnz, PetscInt * const &onz ) const;
};

#endif
