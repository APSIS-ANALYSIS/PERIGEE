#ifndef IPGASSEM_HPP
#define IPGASSEM_HPP
// ============================================================================
// IPGAssem.hpp
// Interface for parallel global assembly.
//
// This class provide an interface for different instantiations of global
// assembly of tangent, mass matrix and residual vector, with different BC
// imposed. 
// 
// The data structure for the matrix and the vector are based on the PETSc Mat
// and Vec objects.
//
// Date: July 31 2015 
// ============================================================================
#include "APart_Node.hpp"
#include "ALocal_Elem.hpp"
#include "AGlobal_Mesh_Info.hpp"
#include "IPLocAssem.hpp"
#include "IPLocAssem_2x2Block.hpp"
#include "FEANode.hpp"
#include "PDNSolution.hpp"
#include "ALocal_NBC.hpp"
#include "ALocal_InflowBC.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "ALocal_Interface.hpp"
#include "Sliding_Interface_Tools.hpp"
#include "IGenBC.hpp"
#include "Tissue_prestress.hpp"

class IPGAssem
{
  public:
    // ------------------------------------------------------------------------
    // Matrix K is a generic matrix object that mainly serves as a SPARSE matrix
    // that we solve with. It can be used as other purposes as well.     
    // ------------------------------------------------------------------------
    Mat K;

    // ------------------------------------------------------------------------
    // Vector G is a generic vector object.
    // ------------------------------------------------------------------------
    Vec G;

    IPGAssem(){};

    virtual ~IPGAssem(){};

    // ------------------------------------------------------------------------
    // ! Flag : Fix nonzero structure
    //          Add or insert in new location is ignored. Set after the first
    //          MatAssemblyEnd().
    // ------------------------------------------------------------------------
    void Fix_nonzero_str()
    {MatSetOption(K, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);}


    // ------------------------------------------------------------------------
    // ! Flag : New allocation error
    //          Add or insert in new locations will generate an error message.
    //          Supports AIJ and BAIJ formats.
    // ------------------------------------------------------------------------
    void Fix_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);}


    // ------------------------------------------------------------------------
    // ! Flag : Ignore new allocation
    //          Add or insert in a new allocation will NOT generate error.
    // ------------------------------------------------------------------------
    void Release_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);}


    // ------------------------------------------------------------------------
    // ! Flag : Keep nonzero pattern of the matrix K
    // ------------------------------------------------------------------------
    void Keep_nonzero_pattern()
    {MatSetOption(K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);}


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
    void Clear_G()
    {VecSet(G, 0.0);}

    // ------------------------------------------------------------------------
    // ! Print the vector G on screen
    // ------------------------------------------------------------------------
    void Print_G() const
    {VecView(G, PETSC_VIEWER_STDOUT_WORLD);}

    // ------------------------------------------------------------------------
    // ! Assem_nonzero_estimate : Assembly nonzero estimate matrix for K.
    //                            Insert 1.0 to every possible nonzero locations.
    // ------------------------------------------------------------------------
    virtual void Assem_nonzero_estimate()
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate(
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &pnode_ptr )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NBC * const &nbc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const ALocal_NBC * const &nbc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    // Nonzero pattern estimate for the NS equations
    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NBC * const &nbc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const APart_Node * const &pnode_v,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    // Nonzero pattern for FSI wall prestressing
    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    // ------------------------------------------------------------------------
    // ! Assem_mass_residual : assembly mass matrix and corresponding residual 
    //                         vector for 3D problems WITHOUT pre-existing 
    //                         cached quadrature info.
    // ------------------------------------------------------------------------
    // Assemble mass matrix and residual vector for NS equations
     virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        const IQuadPts * const &quad_v,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const APart_Node * const &pnode_ptr )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

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
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Tissue_prestress * const &ps_ptr )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

    virtual void Assem_mass_residual(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_part,
        const Tissue_prestress * const &ps_ptr )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementvs,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_WeakBC * const &wbc_part )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &mdisp,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementvs,
        FEAElement * const &elementvs_rotated,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IQuadPts * const &free_quad,
        const FEANode * const &fnode_ptr,
        const ALocal_Interface * const &itf_part,
        const SI_T::SI_solution * const &SI_sol,
        const SI_T::SI_quad_point * const &SI_qp )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

    virtual void Assem_mass_residual( const PDNSolution * const &sol )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}

    // ------------------------------------------------------------------------
    // ! Assem_residual : assembly residual vector for 3D problem WITHOUT
    //                    pre-existing cached quadrature info.
    // ------------------------------------------------------------------------
    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        const IQuadPts * const &quad_v,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const APart_Node * const &pnode_ptr )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
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
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}


    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
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
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}


    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
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
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

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
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}


    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &sol_np1,
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
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}


    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
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
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const Tissue_prestress * const &ps_ptr )
        {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}


    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Tissue_prestress * const &ps_ptr )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    virtual void Assem_Residual(
        const double &curr_time, const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const PDNSolution * const &dot_velo_np1,
        const PDNSolution * const &velo_np1,
        const PDNSolution * const &disp_np1,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const Tissue_prestress * const &ps_ptr )
        {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    // Assembly in the prestress generation
    virtual void Assem_Residual(
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p,
        const Tissue_prestress * const &ps_ptr) 
        {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    // Assembly with weak BC
    virtual void Assem_residual(
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
        FEAElement * const &elementvs,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const ALocal_WeakBC * const &wbc_part )
        {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    // Assembly with interface integral
    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementvs,
        FEAElement * const &elementvs_rotated,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IQuadPts * const &free_quad,
        const FEANode * const &fnode_ptr,
        const IGenBC * const &gbc,
        const ALocal_Interface * const &itf_part,
        const SI_T::SI_solution * const &SI_sol,
        const SI_T::SI_quad_point * const &SI_qp )
        {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    virtual void Assem_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const double &curr_time,
        const double &dt )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}

    // ------------------------------------------------------------------------
    // ! Assem_tangent_residual : assembly tangent matrix and residual vector 
    //                            for 3D problem WITHOUT pre-existing cached 
    //                            quadrature info.
    // ------------------------------------------------------------------------
    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
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
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}


    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
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
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}


    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
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
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
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
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

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
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &sol_np1,
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
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
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
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const Tissue_prestress * const &ps_ptr )
        {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const Tissue_prestress * const &ps_ptr )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    virtual void Assem_Tangent_Residual(
        const double &curr_time, const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const PDNSolution * const &dot_velo_np1,
        const PDNSolution * const &velo_np1,
        const PDNSolution * const &disp_np1,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_f_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const Tissue_prestress * const &ps_ptr )
        {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    // Assembly in prestress generation
    virtual void Assem_Tangent_Residual(
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &dot_disp,
        const PDNSolution * const &dot_velo,
        const PDNSolution * const &dot_pres,
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        const PDNSolution * const &pres,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_v,
        const ALocal_NBC * const &nbc_p,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p,
        const Tissue_prestress * const &ps_ptr )
        {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    // Assembly with weak BC
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
        FEAElement * const &elementvs,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc,
        const ALocal_WeakBC * const &wbc_part )
        {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    // Assembly with interface integral
    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const PDNSolution * const &mvelo,
        const PDNSolution * const &mdisp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementvs,
        FEAElement * const &elementvs_rotated,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        IQuadPts * const &free_quad,
        const FEANode * const &fnode_ptr,
        const IGenBC * const &gbc,
        const ALocal_Interface * const &itf_part,
        const SI_T::SI_solution * const &SI_sol,
        const SI_T::SI_quad_point * const &SI_qp )
        {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    virtual void Assem_tangent_residual(
        const PDNSolution * const &dot_sol,
        const PDNSolution * const &sol,
        const double &curr_time,
        const double &dt )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}

    // --------------------------------------------------------------
    // Assembly boundary integrals
    // --------------------------------------------------------------
    // Assem_surface_flowrate
    // Performs surface integral to calculate the flow rate on the face
    // with ebc_id.
    // Integration is performed for the ebc_part(ebc_id) domain for
    //        v_x n_x + v_y n_y + v_z n_z
    // and then a MPI_Reduce is called to collect the value over multiple
    // CPUs    
    virtual double Assem_surface_flowrate(
        const PDNSolution * const &vec,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const int &ebc_id )
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &vec,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id )
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id )
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id )
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &vec,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id )
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &vec,
        const int &ebc_id )
    {
      SYS_T::commPrint("Warning: IPGAssem::Assem_surface_flowrate is not implemented. \n");
      return 0.0;
    }

    // Assem_surface_ave_pressure
    // Performs surface integral to calculated the pressure integrated
    // over the surface as well as the surface area. Return the
    // pressure integration divided by the area to give a surface
    // averaged pressure. Before division, the values will be collected
    // from multiple CPUs by MPI_Reduce.
    //  integration p dGamma / integration 1 dGamma
    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &vec,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const int &ebc_id )
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &vec,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id )
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &vec,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id )
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &vec,
        const int &ebc_id )
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_v,
        const ALocal_EBC * const &ebc_p,
        const int &ebc_id )
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented. \n");
      return 0.0;
    }

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_InflowBC * const &infbc_part,
        const int &nbc_id )
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented. \n");
      return 0.0;
    }

    // Update wall prestress at all surface quadrature points
    virtual void Update_Wall_Prestress(
        const PDNSolution * const &sol_wall_disp,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_w,
        const IQuadPts * const &quad_s,
        ALocal_EBC * const &ebc_wall_part )
    {SYS_T::commPrint("Warning: Update_Wall_Prestress() is not implemented. \n");}

    // Update solid prestress at all volumetric quadrature points (in
    // tet4_vascular) 
    virtual void Update_Wall_Prestress(
        const PDNSolution * const &sol,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element,
        const IQuadPts * const &quad,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        Tissue_prestress * const &ps_ptr ) const
    {SYS_T::commPrint("Warning: Update_Wall_Prestress() is not implemented. \n");}

    // Update solid prestress at all volumetric quadrature points (in tet4_fsi) 
    virtual void Update_Wall_Prestress(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem_2x2Block * const &lassem_s_ptr,
        FEAElement * const &elementv,
        const IQuadPts * const &quadv,
        const ALocal_IEN * const &lien_v,
        const ALocal_IEN * const &lien_p,
        const FEANode * const &fnode_ptr,
        Tissue_prestress * const &ps_ptr ) const
    {SYS_T::commPrint("Warning: Update_Wall_Prestress() is not implemented. \n");}
};

#endif
