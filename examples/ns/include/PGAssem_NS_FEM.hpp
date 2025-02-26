#ifndef PGASSEM_NS_FEM_HPP
#define PGASSEM_NS_FEM_HPP
// ==================================================================
// PGAssem_NS_FEM.hpp
//
// Parallel golbal assembly based on PETSc, using AIJ matrix format.
// The assembly routine is designed for classical C0 FEM method, which
// means we do not need extraction operators and local mesh sizes.
//
// The assembly is for the NS equations written in VMS formualtion
// of NS equations. The input solution vectors contains
//  [ pressure; velocity ],
// the dot solution contains
//  [ dot pressure; dot velcoty ].
//
// Author: Ju Liu 
// Date Created: Feb. 10 2020
// ==================================================================
#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"
#include "PDNSolution_NS.hpp"

class PGAssem_NS_FEM : public IPGAssem
{
  public:
    // Constructor for NS equations
    PGAssem_NS_FEM( 
        const IGenBC * const &gbc,
        std::unique_ptr<ALocal_IEN> in_locien,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode,
        std::unique_ptr<ALocal_NBC> in_nbc,
        std::unique_ptr<ALocal_EBC> in_ebc,
        std::unique_ptr<ALocal_WeakBC> in_wbc,
        std::unique_ptr<IPLocAssem> in_locassem,    
        const int &in_nz_estimate=60 );

    // Destructor
    virtual ~PGAssem_NS_FEM();

    // Nonzero pattern estimate for the NS equations
    virtual void Assem_nonzero_estimate(
        const IGenBC * const &gbc );

    // Assem mass matrix and residual vector
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a );

    // Assembly the residual vector for the NS equations
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const IGenBC * const &gbc );

    // Assembly the residual vector and tangent matrix 
    // for the NS equations
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const IGenBC * const &gbc );

    // Assembly routine for the surface integrals of flow rate and
    // pressure
    virtual double Assem_surface_flowrate(
        const PDNSolution * const &sol,
        const int &ebc_id ) const;

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &sol,
        const ALocal_InflowBC * const &infbc_part,
        const int &infnbc_id ) const;

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &sol,
        const int &ebc_id ) const;

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &sol,
        const ALocal_InflowBC * const &infbc_part,
        const int &infnbc_id ) const;

  private:
    // Private data
    const std::unique_ptr<const ALocal_IEN> locien;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode;
    const std::unique_ptr<const ALocal_NBC> nbc;
    const std::unique_ptr<const ALocal_EBC> ebc;
    const std::unique_ptr<const ALocal_WeakBC> wbc;
    const std::unique_ptr<IPLocAssem> locassem;

    const int nLocBas, snLocBas, dof_sol, dof_mat, num_ebc, nlgn;

    // Private function
    // Essential boundary condition
    void EssBC_KG( const int &field );
    
    void EssBC_G( const int &field );

    // Natural boundary condition
    void NatBC_G( const double &curr_time, const double &dt );

    // Backflow integral on outlet surfaces
    void BackFlow_G( const PDNSolution * const &sol );

    void BackFlow_KG( const double &dt,
        const PDNSolution * const &sol );

    // Resistance type boundary condition on outlet surfaces
    void NatBC_Resis_G( const double &curr_time, const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const IGenBC * const &gbc );

    void NatBC_Resis_KG( const double &curr_time, const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const IGenBC * const &gbc );

    // Weak imposition of no-slip boundary condition on wall
    void Weak_EssBC_KG( const double &curr_time, const double &dt,
        const PDNSolution * const &sol );

    void Weak_EssBC_G( const double &curr_time, const double &dt,
        const PDNSolution * const &sol );

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
};

#endif
