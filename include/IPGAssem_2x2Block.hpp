#ifndef IPGASSEM_2X2BLOCK_HPP
#define IPGASSEM_2X2BLOCK_HPP
// ==================================================================
// IPGAssem_2x2Block.hpp
// Interface for parallel global assembly for a 2-by-2 block 
// structured matrices.
//
// The data structure for the matrices and vectors are the PETSc
// Mat and Vec objects, which are based on the same MPI partition.
// 
// [ K_00, K_01    [x_0   = [ G_0
//   K_10, K_11 ]   x_1]      G_1 ].
//
// Author: Ju Liu
// Date: Jan 18 2018
// ==================================================================
#include "PETSc_Tools.hpp"
#include "AGlobal_Mesh_Info.hpp"
#include "IPLocAssem_2x2Block.hpp"
#include "ALocal_NBC.hpp"
#include "ALocal_NodalBC_2x2Block.hpp"
#include "ALocal_EBC.hpp"
#include "FEANode.hpp"
#include "PDNSolution.hpp"

class IPGAssem_2x2Block
{
  public:
    Mat K_00;
    Mat K_01;
    Mat K_10;
    Mat K_11;

    Vec G_0;
    Vec G_1;

    // Constructor
    IPGAssem_2x2Block(){};

    // Destructor
    virtual ~IPGAssem_2x2Block(){};

    // --------------------------------------------------------------
    // Fix nonzero structure
    // Add or insert in new location is ignored. This flag will be set
    // after the first MatAssemblyEnd.
    // --------------------------------------------------------------
    void Fix_nonzero_str()
    {
      MatSetOption(K_00, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(K_01, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(K_10, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
      MatSetOption(K_11, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
    }

    // --------------------------------------------------------------
    // New allocation error
    // Add or insert in new allocations will generate an error message.
    // --------------------------------------------------------------
    void Fix_nonzero_err_str()
    {
      MatSetOption(K_00, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(K_01, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(K_10, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
      MatSetOption(K_11, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    }

    // --------------------------------------------------------------
    // Ignore new allocation
    // Add or insert in new allocations will NOT generate error
    // --------------------------------------------------------------
    void Release_nonzero_err_str()
    {
      MatSetOption(K_00, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(K_01, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(K_10, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetOption(K_11, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }

    // --------------------------------------------------------------
    // Keep nonzero pattern of the matrices
    // --------------------------------------------------------------
    void Keep_nonzero_pattern()
    {
      MatSetOption(K_00, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(K_01, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(K_10, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(K_11, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    }

    // --------------------------------------------------------------
    // Clear K matrices and G vectors
    // --------------------------------------------------------------
    void Clear_KG()
    {
      MatZeroEntries(K_00); MatZeroEntries(K_01);
      MatZeroEntries(K_10); MatZeroEntries(K_11);
      VecSet(G_0, 0.0);     VecSet(G_1, 0.0);
    }

    // --------------------------------------------------------------
    // Clear G vectors by assign zero values
    // --------------------------------------------------------------
    void Clear_G()
    {
      VecSet(G_0, 0.0);     
      VecSet(G_1, 0.0);
    }

    
    // --------------------------------------------------------------
    // Print functions
    // --------------------------------------------------------------
    void Print_G_0() const {VecView(G_0, PETSC_VIEWER_STDOUT_WORLD);}
    
    void Print_G_1() const {VecView(G_1, PETSC_VIEWER_STDOUT_WORLD);}

    void Print_K_00() const {MatView(K_00, PETSC_VIEWER_STDOUT_WORLD);}
    
    void Print_K_01() const {MatView(K_01, PETSC_VIEWER_STDOUT_WORLD);}
    
    void Print_K_10() const {MatView(K_10, PETSC_VIEWER_STDOUT_WORLD);}
    
    void Print_K_11() const {MatView(K_11, PETSC_VIEWER_STDOUT_WORLD);}

    // --------------------------------------------------------------
    // Calculate the residual norms
    //  norm_0 = \| G0 - K01 sol_1 - K00 sol_0 \|_l2
    //  norm_1 = \| G1 - K11 sol_1 - K10 sol_0 \|_l2
    // --------------------------------------------------------------
    virtual void Get_res_norms( const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        double &norm_0, double &norm_1 ) const
    {
      SYS_T::commPrint("Warning: Get_res_norms() is not implemented.\n");
    }


    // --------------------------------------------------------------
    // Assem_nonzero_estimate : Assembly nonzero estimate for matrices
    // by inserting 1.0 to every possible nonzero locations.
    // --------------------------------------------------------------
    virtual void Assem_nonzero_estimate(
        IPLocAssem_2x2Block * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NBC * const &nbc_part )
    {
      SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented.\n");
    }

    // --------------------------------------------------------------
    // Assem_mass_residual : Assembly mass matrices and corresponding
    // residual vectors.
    // --------------------------------------------------------------
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

    
    // Interface for Segregated elastodynamic with disp, pres, and velo 
    // as solution input.
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        const PDNSolution * const &sol_2,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

    // --------------------------------------------------------------
    // Assem_residual : Assembly residual vectors.
    // --------------------------------------------------------------
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol_0,
        const PDNSolution * const &dot_sol_1,
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    
    // Interface for Segregated elastodynamic with disp, pres, and velo 
    // as solution input.
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol_0,
        const PDNSolution * const &dot_sol_1,
        const PDNSolution * const &dot_sol_2,
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        const PDNSolution * const &sol_2,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    // --------------------------------------------------------------
    // Assem_tangent_residual : Assembly tangent matrices and residual
    // vectors.
    // --------------------------------------------------------------
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol_0,
        const PDNSolution * const &dot_sol_1,
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    
    // Interface for Segregated elastodynamic with disp, pres, and velo 
    // as solution input.
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol_0,
        const PDNSolution * const &dot_sol_1,
        const PDNSolution * const &dot_sol_2,
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        const PDNSolution * const &sol_2,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}


    // Interface for NonIsoparametric FEM
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol_0,
        const PDNSolution * const &dot_sol_1,
        const PDNSolution * const &sol_0,
        const PDNSolution * const &sol_1,
        const double &curr_time,
        const double &dt,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_0_ptr,
        const ALocal_IEN * const &lien_1_ptr,
        const ALocal_IEN * const &lien_geo_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC_2x2Block * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}


};

#endif
