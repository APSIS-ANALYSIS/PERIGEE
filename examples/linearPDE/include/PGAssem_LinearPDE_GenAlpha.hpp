#ifndef PGASSEM_LINEARPDE_GENALPHA_HPP
#define PGASSEM_LINEARPDE_GENALPHA_HPP

#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"
#include "PDNSolution_Elastodynamics.hpp"
#include "PDNSolution_Transport.hpp"

class PGAssem_LinearPDE_GenAlpha : public IPGAssem
{
  public:
    // Constructor for equations
    PGAssem_LinearPDE_GenAlpha(
        IPLocAssem * const &locassem_ptr,
        const AGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const ALocal_NBC * const &part_nbc,
        const ALocal_EBC * const &part_ebc,
        const int &in_nz_estimate = 60 );

    // Destructor
    virtual ~PGAssem_LinearPDE_GenAlpha();

    // Nonzero pattern estimate
    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NBC * const &nbc_part );

    // Assembly the residual vector
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    // Assembly the residual vector and tangent matrix 
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    // Assembly the residual and mass matrix
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

  private:
    // Private data
    const int num_ebc, nLocBas, snLocBas, dof_mat, nlgn;

    // Private function
    // Essential boundary condition
    void EssBC_KG( const ALocal_NBC * const &nbc_part, const int &field );

    void EssBC_G( const ALocal_NBC * const &nbc_part, const int &field );

    // Natural boundary condition
    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    void GetLocal( const double * const &array, const int * const &IEN,
        double * const &local_array ) const
    {
      for(int ii=0; ii<nLocBas; ++ii)
      {
        const int offset1 = ii * dof_mat;
        const int offset2 = IEN[ii] * dof_mat;
        for(int jj=0; jj<dof_mat; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }
};

#endif
