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
#include "IAGlobal_Mesh_Info.hpp"
#include "IALocal_BC.hpp"
#include "IPLocAssem.hpp"
#include "FEANode.hpp"
#include "PDNSolution.hpp"
#include "ALocal_meshSize_3D_NURBS.hpp"
#include "AExtractor_3D_NURBS_xyz.hpp"
#include "ALocal_NodalBC.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "ALocal_EBC.hpp"
#include "IGenBC.hpp"

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
    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const IALocal_BC * const &bc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}


    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    virtual void Assem_nonzero_estimate( 
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_f_ptr,
        IPLocAssem * const &lassem_s_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_nonzero_estimate() is not implemented. \n");}

    // ------------------------------------------------------------------------
    // ! Assem_mass_residual : assembly mass matrix and corresponding residual 
    //                         vector for 3D problems WITHOUT pre-existing 
    //                         cached quadrature info.
    // ------------------------------------------------------------------------
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr, 
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const AInt_Weight * const &wei_ptr,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const IALocal_BC * const &bc_part )
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
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part )
    {SYS_T::commPrint("Warning: Assem_mass_residual() is not implemented. \n");}


    // ------------------------------------------------------------------------
    // ! Assem_residual : assembly residual vector for 3D problem WITHOUT
    //                    pre-existing cached quadrature info.
    // ------------------------------------------------------------------------
    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr, 
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const AInt_Weight * const &wei_ptr,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const IALocal_BC * const &bc_part )
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
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
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
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_residual() is not implemented. \n");}


    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &sol_wall_disp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementw,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_wall_part,
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
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
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
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const AInt_Weight * const &wei_ptr,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const IALocal_BC * const &bc_part )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}


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
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
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
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const IGenBC * const &gbc )
    {SYS_T::commPrint("Warning: Assem_tangent_residual() is not implemented. \n");}
    

    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const PDNSolution * const &sol_wall_disp,
        const PDNSolution * const &dot_sol_np1,
        const PDNSolution * const &sol_np1,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &elementv,
        FEAElement * const &elements,
        FEAElement * const &elementw,
        const IQuadPts * const &quad_v,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
        const ALocal_EBC * const &ebc_wall_part,
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
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
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
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part,
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
        const ALocal_EBC * const &ebc_part,
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
        const ALocal_Inflow_NodalBC * const &infbc_part )
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
        const ALocal_EBC * const &ebc_part,
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
        const ALocal_Inflow_NodalBC * const &infbc_part )
    {
      SYS_T::commPrint("Warning: Assem_surface_ave_pressure is not implemented. \n");
      return 0.0;
    }
};

#endif
