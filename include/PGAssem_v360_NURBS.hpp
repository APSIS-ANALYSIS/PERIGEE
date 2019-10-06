#ifndef PGASSEM_V360_NURBS_HPP
#define PGASSEM_V360_NURBS_HPP
// ==================================================================
// PGAssem_v360_NURBS.hpp
// Parallel global assembly based on PETSc-3.6.0, using AIJ matrix format.
// The quadrature info for NURBS basis function is computed within 
// the local assembly code.
//
// This code is an improvement of the PGAssem_v360_noCache_NURBS:
//  1. The boundary condition is split into nodalBC and elemBC in input:
//         bc_part -> nbc_part and ebc_part
//  2. The imposition of essential BC sets all bc values in one setvalues
//     function call.
//  3. The imposition of natural BC is in the NatBC_G assembly routine. 
//     Inside this function, the element object is resized to account 
//     for the boundary quadrature rule. Once the boundary integraiton 
//     is finished, the element is resized back to the volumetric 
//     integration setting. 
//
// Date: Oct. 13 2015
// ==================================================================
#include <ctime>
#include "IPGAssem.hpp"

class PGAssem_v360_NURBS : public IPGAssem
{
  public:
    PGAssem_v360_NURBS(
        IPLocAssem const * const &locassem_ptr,
        IAGlobal_Mesh_Info const * const &agmi_ptr,
        ALocal_Elem const * const &alelem_ptr,
        ALocal_IEN const * const &aien_ptr,
        APart_Node const * const &pnode_ptr,
        ALocal_NodalBC const * const &part_nbc );

    virtual ~PGAssem_v360_NURBS();


    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part );

    // ------------------------------------------------------------------------
    // Assem_mass_residual : Assembly the mass matrix and the corresponding
    //                       Residual vector. The element is passed essentially
    //                       as a container. The quadrature rules are built 
    //                       inside the local assembly code.
    //
    //                       This assembly is initially designed for my TNSK 
    //                       2D & 3D based on B-spline discretization
    // ------------------------------------------------------------------------
    virtual void Assem_mass_residual(
        const PDNSolution * const &sol_a,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const AInt_Weight * const &wei_ptr,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_ElemBC * const &ebc_part );


    virtual void Assem_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const AInt_Weight * const &wei_ptr,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_ElemBC * const &ebc_part );



    virtual void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
        const double &curr_time,
        const double &dt,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const AInt_Weight * const &wei_ptr,
        const IALocal_meshSize * const &mSize,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const IAExtractor * const &extractor,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_ElemBC * const &ebc_part );


  private:
    int nLocBas, dof;

    PetscInt * row_index;

    double * array_a;
    double * array_b;

    double * local_a;
    double * local_b;

    int * IEN_e;

    double * ectrl_x;
    double * ectrl_y;
    double * ectrl_z;
    double * ectrl_w;

    
    // -------------------------------------------------------------------
    // Natural (elemental) boundary conditions
    // This integrate -int_{\partial Omega} w_i H_i dA and provide
    // update for the residual vector.
    // Note: This function only updates the residual vector G.
    //       IN this function, the element object will resize its
    //       memory allocation to prepare for the boundary integration.
    //       At the end of the funciton, the element is resized back
    //       for volumetric integration.
    // -------------------------------------------------------------------
    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element,
        const int &in_loc_dof,
        const int &nqpx, const int &nqpy, const int &nqpz,
        const ALocal_IEN * const &lien_ptr,
        const FEANode * const &fnode_ptr,
        const IALocal_meshSize * const &mSize,
        const IAExtractor * const &extractor,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_ElemBC * const &ebc_part );


    // Essential (nodal) boundary conditions
    void EssBC_KG( const ALocal_NodalBC * const &nbc_part, const int &field );
    void EssBC_G( const ALocal_NodalBC * const &nbc_part, const int &field );


    void GetLocal(const double * const &array, const int * const &IEN,
        double * const &local_array) const
    {
      int offset1, offset2;
      for(int ii=0; ii<nLocBas; ++ii)
      {
        offset1 = ii * dof;
        offset2 = IEN[ii] * dof;
        for(int jj=0; jj<dof; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }

    
    int Get_part_id( const std::vector<unsigned int> &nlist, 
        const unsigned int &val ) const
    {
      for(unsigned int ii=1; ii<nlist.size(); ++ii)
      {
        if(val < nlist[ii])
          return int(ii-1);
      }
      return -1;
    }


    // Get_dnz_onz
    // get a row-by-row estimate of the number of nonzeros per row.
    // This estimate should be accurate for regular, dirichlet, and slave nodes;
    // This estimate may overestimate the number of nz for periodic master row.
    void Get_dnz_onz( const int &nElem,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part,
        PetscInt * const &dnz, PetscInt * const &onz ) const;

    
    // Get_dnz_onz
    // This is a rather simplified and fast version of sparsity estimator.
    // For regular nodes, we assign dof * (2sdeg+1) * (2tdeg+1) * (2udeg+1) for
    // dnz and onz;
    // for dirichlet nodes, we assign 1 for dnz and onz;
    // for slave nodes, we assign 2 for dnz and onz;
    // for master nodes, we assign (1+number of
    // slaves)*dof*(2sdeg+1)*(2tdeg+1)*(2udeg+1) for dnz and onz
    void Get_dnz_onz( const int &nlocnode,
        const int &insdeg, const int &intdeg, const int &inudeg,
        const ALocal_NodalBC * const &nbc_ptr,
        PetscInt * const &dnz, PetscInt * const &onz ) const;

};


#endif
