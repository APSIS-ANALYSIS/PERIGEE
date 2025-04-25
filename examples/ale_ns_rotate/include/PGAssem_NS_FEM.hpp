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
#include "ALocal_Elem.hpp"
#include "ALocal_NBC.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "PETSc_Tools.hpp"
#include "PDNSolution_NS.hpp"
#include "FE_Tools.hpp"
#include "FEAElementFactory.hpp"
#include "QuadPtsFactory.hpp"
#include "Sliding_Interface_Tools.hpp"

class PGAssem_NS_FEM : public IPGAssem
{
  public:
    // Constructor for NS equations
    PGAssem_NS_FEM(
        const FEType &in_type,
        const double &in_nqps,
        std::unique_ptr<ALocal_IEN> in_locien,
        std::unique_ptr<ALocal_Elem> in_locelem,
        std::unique_ptr<FEANode> in_fnode,
        std::unique_ptr<APart_Node> in_pnode,
        std::unique_ptr<ALocal_NBC> in_nbc,
        std::unique_ptr<ALocal_EBC> in_ebc,
        std::unique_ptr<ALocal_WeakBC> in_wbc,
        std::unique_ptr<ALocal_Interface> in_itf,
        std::unique_ptr<IPLocAssem> in_locassem,
        std::unique_ptr<SI_T::SI_solution> in_SI_sol,
        std::unique_ptr<SI_T::SI_quad_point> in_SI_qp,
        const IGenBC * const &gbc,
        const int &in_nz_estimate=60 );

    // Destructor
    virtual ~PGAssem_NS_FEM();

    // Nonzero pattern estimate for the NS equations
    virtual void Assem_nonzero_estimate( 
        const IGenBC * const &gbc );

    // Assem mass matrix and residual vector
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &mdisp);

    // Assembly the residual vector for the NS equations
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp, 
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
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp,
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
        const int &nbc_id ) const;

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &sol,
        const int &ebc_id ) const;

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &sol,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id ) const;

    virtual void Interface_K_MF(Vec &X, Vec &Y);

    virtual void Update_SI_state(
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp )
    {
      SI_sol->update_node_sol(sol);
      SI_sol->update_node_mvelo(mvelo);
      SI_sol->update_node_mdisp(mdisp);
      SI_qp->search_all_opposite_point(itf.get(), SI_sol.get());
    }

    virtual void Update_SI_sol(
        const PDNSolution * const &sol )
    { SI_sol->update_node_sol(sol); }

    virtual const FEANode * Get_fnode()
    { return fnode.get(); }

    virtual const APart_Node * Get_pnode()
    { return pnode.get(); }

  private:
    // Private data
    const std::unique_ptr<SI_T::SI_solution> SI_sol;
    const std::unique_ptr<SI_T::SI_quad_point> SI_qp;

    const std::unique_ptr<const ALocal_IEN> locien;
    const std::unique_ptr<const ALocal_Elem> locelem;
    const std::unique_ptr<const FEANode> fnode;
    const std::unique_ptr<const APart_Node> pnode;
    const std::unique_ptr<const ALocal_NBC> nbc;
    const std::unique_ptr<const ALocal_EBC> ebc;
    const std::unique_ptr<const ALocal_WeakBC> wbc;
    const std::unique_ptr<const ALocal_Interface> itf;
    const std::unique_ptr<IPLocAssem> locassem;

    const int nLocBas, snLocBas, dof_sol, dof_mat, num_ebc, nlgn;

    // For the interface integral
    const std::unique_ptr<FEAElement> anchor_elementv;

    // Defined with only one quadrature point,
    // given by the found closest point
    const std::unique_ptr<FEAElement> opposite_elementv;

    // For the interface integral
    const std::unique_ptr<const IQuadPts> quad_s;

    // Free parametric surface quadrature point
    const std::unique_ptr<IQuadPts> free_quad;

    double time_step;

    // Private function
    // Essential boundary condition
    void EssBC_KG( const int &field );
    
    void EssBC_G( const int &field );

    // Natural boundary condition
    void NatBC_G( const double &curr_time, const double &dt );

    // Backflow integral on outlet surfaces
    void BackFlow_G( const PDNSolution * const &sol );

    void BackFlow_KG( const double &dt, const PDNSolution * const &sol );

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
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp);

    void Weak_EssBC_G( const double &curr_time, const double &dt,
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp);

    virtual void Interface_KG(
        const double &dt );

    virtual void Interface_G(
        const double &dt );

    virtual void local_MatMult_MF(
        const int &dof,
        PetscInt * &row_index,
        PetscInt * &col_index,
        PetscScalar * &Mat,
        Vec &X,
        Vec &Y );

    void GetLocal(const double * const &array, const int * const &IEN,
        double * const &local_array) const
    {
      GetLocalImpl( array, nLocBas, dof_sol, IEN, local_array );
    }

    void GetLocal( const double * const &array, const int * const &IEN,
        const int &in_locbas, double * const &local_array) const
    {
      GetLocalImpl( array, in_locbas, dof_sol, IEN, local_array );
    }

    void GetLocal(const double * const &array, const int &in_dof, 
        const int * const &IEN, double * const &local_array) const
    {
      GetLocalImpl( array, nLocBas, in_dof, IEN, local_array );
    }

    void GetLocalImpl( const double * const &array, 
        const int &in_locbas, const int &in_dof,
        const int * const &IEN, double * const &local_array ) const
    {
      for(int ii=0; ii<in_locbas; ++ii)
      {
        const int offset1 = ii * in_dof;
        const int offset2 = IEN[ii] * in_dof;
        for(int jj=0; jj<in_dof; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

    void get_currPts( const double * const &ept_x,
        const double * const &ept_y,
        const double * const &ept_z,
        const double * const &disp,
        const int &len,
        double * const &currPt_x,
        double * const &currPt_y,
        double * const &currPt_z ) const
    {
      for(int ii=0; ii<len; ++ii)
      {
        currPt_x[ii] = ept_x[ii] + disp[3*ii];
        currPt_y[ii] = ept_y[ii] + disp[3*ii+1];
        currPt_z[ii] = ept_z[ii] + disp[3*ii+2];
      }
    }
};

#endif
