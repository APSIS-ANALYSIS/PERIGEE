#ifndef PGASSEM_MESH_HPP
#define PGASSEM_MESH_HPP
// ============================================================================
// PGAssem_Mesh.hpp
//
// Parallel Global Assembly for Mesh motion equations.
//
// Author: Ju Liu
// Date: Dec. 30 2021
// ============================================================================
#include "IPGAssem.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_NBC.hpp"
#include "PETSc_Tools.hpp"

class PGAssem_Mesh : public IPGAssem
{
  public:
    PGAssem_Mesh(
        std::unique_ptr<ALocal_IEN> in_locien_v,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode_v,
        std::unique_ptr<ALocal_NBC> in_mesh_nbc,
        std::unique_ptr<ALocal_EBC> in_mesh_ebc,
        std::unique_ptr<IPLocAssem> in_locassem,
        const int &in_nz_estimate );

    virtual ~PGAssem_Mesh();

    virtual void Assem_nonzero_estimate();

    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a );

    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt );

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt );

  private:
    const std::unique_ptr<const ALocal_IEN> locien;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode;
    const std::unique_ptr<const ALocal_NBC> mesh_nbc;
    const std::unique_ptr<const ALocal_EBC> mesh_ebc;
    const std::unique_ptr<IPLocAssem> locassem;

    const int nLocBas, snLocBas, dof, num_ebc, nlgn;
    
    void EssBC_KG();

    void EssBC_G();

    void NatBC_G( const double &curr_time, const double &dt );

    std::vector<double> GetLocal( const std::vector<double> &array,
        const std::vector<int> &IEN, const int &in_locbas ) const
    {
      std::vector<double> out( in_locbas * dof, 0.0 );
      for(int ii=0; ii<in_locbas; ++ii)
        for(int jj=0; jj<dof; ++jj)
          out[ii*dof + jj] = array[ IEN[ii] * dof + jj ];
      
      return out;
    }

};

#endif
