#ifndef PGASSEM_2X2BLOCK_NS_FEM_HPP
#define PGASSEM_2X2BLOCK_NS_FEM_HPP
// ==================================================================
// PGAssem_2x2Block_NS_FEM.hpp
//
// Parallel global assembly using nest matrix.
// The assembly is designed for C0 Lagrangian elements, which does
// not need extraction operator or local mesh sizes.
//
// Author: Ju Liu
// Date Created: Aug. 11 2020
// ==================================================================
#include "APart_Node.hpp"
#include "ALocal_Elem.hpp"
#include "IAGlobal_Mesh_Info.hpp"
#include "IPLocAssem_2x2Block.hpp"
#include "PDNSolution_NS.hpp"
#include "ALocal_NBC.hpp"
#include "ALocal_EBC.hpp"
#include "IGenBC.hpp"
#include "PETSc_Tools.hpp"

class PGAssem_2x2Block_NS_FEM
{
  public:
    Mat K;
    Mat subK[4];
    
    Vec G;
    Vec subG[2];

    IS is[2];

    // Constructor
    PGAssem_2x2Block_NS_FEM(
        IPLocAssem_2x2Block * const &locassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quads,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const ALocal_NBC * const &part_nbc,
        const ALocal_EBC * const &part_ebc,
        const IGenBC * const &gbc,
        const int &in_nz_estimate=60 );

    // Destructor
    virtual ~PGAssem_2x2Block_NS_FEM();

    // ------------------------------------------------------------------------
    // ! Flag : Fix nonzero structure
    //          Add or insert in new location is ignored. Set after
    //         the first MatAssemblyEnd().
    // ------------------------------------------------------------------------
    void Fix_nonzero_str()
    {
      MatSetOption(subK[0], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(subK[1], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(subK[2], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(subK[3], MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
    }

    // ------------------------------------------------------------------------
    // ! Flag : New allocation error
    //          Add or insert in new locations will generate an error message.
    //          Supports AIJ and BAIJ formats.
    // ------------------------------------------------------------------------
    void Fix_nonzero_err_str()
    {
      MatSetOption(subK[0], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(subK[1], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(subK[2], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(subK[3], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    }
    
    // ------------------------------------------------------------------------
    // ! Flag : Ignore new allocation
    //          Add or insert in a new allocation will NOT generate error.
    // ------------------------------------------------------------------------
    void Release_nonzero_err_str()
    {
      MatSetOption(subK[0], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(subK[1], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(subK[2], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(subK[3], MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }

    // ------------------------------------------------------------------------
    // ! Flag : Keep nonzero pattern of the matrix K
    // ------------------------------------------------------------------------
    void Keep_nonzero_pattern()
    {
      MatSetOption(subK[0], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(subK[1], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(subK[2], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(subK[3], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    }

    // ------------------------------------------------------------------------
    // ! Clear K and G to be zero
    // ------------------------------------------------------------------------
    void Clear_KG()
    {
      MatZeroEntries(K);
      VecSet(G, 0.0);
    }

    // ------------------------------------------------------------------------
    // ! Clear G to be zero
    // ------------------------------------------------------------------------
    void Clear_G() {VecSet(G, 0.0);}

    // Nonzero pattern for the NS equation    
    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    
    // Assembly mass matrix and residual vector
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );


    // Assembly residual vector
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );


    // Assembly tangent matrix and residual vector
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );


    // Assembly routine for the surface integrals of flow rate
    virtual double Assem_surface_flowrate(
        const PDNSolution * const &sol,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id );

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &sol,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part );

    // Assembly routine for the surface averaged pressure
    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &sol,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id );

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &sol,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part );

  private:
    // Data:
    const int nLocBas, dof_sol, dof_mat_v, dof_mat_p, num_ebc, nlgn;

    int snLocBas;

    // Functions:
    void EssBC_KG( const ALocal_NBC * const &nbc_part );

    void EssBC_G( const ALocal_NBC * const &nbc_part );

    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    void BackFlow_G( const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    void BackFlow_KG( const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );

    void NatBC_Resis_G( const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

    void NatBC_Resis_KG( const double &dt,
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc );

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
