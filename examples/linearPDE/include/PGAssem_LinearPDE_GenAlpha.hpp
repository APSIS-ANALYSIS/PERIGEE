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
        std::unique_ptr<ALocal_IEN> in_locien,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode,
        std::unique_ptr<ALocal_NBC> in_nbc,
        std::unique_ptr<ALocal_EBC> in_ebc,
        std::unique_ptr<IPLocAssem> in_locassem,
        const int &in_nz_estimate = 60 );

    // Destructor
    virtual ~PGAssem_LinearPDE_GenAlpha();

    // Nonzero pattern estimate
    virtual void Assem_nonzero_estimate();

    // Assembly the residual vector
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const double &curr_time,
        const double &dt );

    // Assembly the residual vector and tangent matrix 
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const double &curr_time,
        const double &dt );

    // Assembly the residual and mass matrix
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol );

  private:
    // Private data
    const std::unique_ptr<const ALocal_IEN> locien;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode;
    const std::unique_ptr<const ALocal_NBC> nbc;
    const std::unique_ptr<const ALocal_EBC> ebc;
    const std::unique_ptr<IPLocAssem> locassem;

    const int num_ebc, nLocBas, snLocBas, dof_mat, nlgn;

    // Private function
    // Essential boundary condition
    void EssBC_KG( const int &field );

    void EssBC_G( const int &field );

    // Natural boundary condition
    void NatBC_G( const double &curr_time, const double &dt );

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
