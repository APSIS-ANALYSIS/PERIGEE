#ifndef PGASSEM_SV_NS_FEM_HPP
#define PGASSEM_SV_NS_FEM_HPP
// ==================================================================
// PGAssem_SV_NS_FEM.hpp
// 
// Parallel global assembly based on PETSc, using AIJ matrix
// format. The quadrature info for basis functions is computed within
// the local assembly routine.
// 
// This assembly routine is designed for classical C0 FEM method, (
// instead of the IGA method.)
// 
// We utilize the function pointers for the local assembly member
// function to facilitate the implementation of natural boundary
// conditions.
//
// Date Crated: June 08 2018
// ==================================================================
#include "IPGAssem.hpp"
#include "PETSc_Tools.hpp"

class PGAssem_SV_NS_FEM : public IPGAssem
{
  public:
    PGAssem_SV_NS_FEM(
        IPLocAssem * const &locassem_ptr,
        IAGlobal_Mesh_Info const * const &agmi_ptr,
        ALocal_Elem const * const &alelem_ptr,
        ALocal_IEN const * const &aien_ptr,
        APart_Node const * const &pnode_ptr,
        ALocal_NodalBC const * const &part_nbc,
        ALocal_EBC const * const &part_ebc );

    virtual ~PGAssem_SV_NS_FEM();

    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part );

    // --------------------------------------------------------------
    // Assem_mass_residual: Assembly the mass matrix and the residual
    // vector. The volumetric and surface element area passed into
    // the routine separately; the corresponding quadrature rules for
    // volume and surface integration are passed accordingly.
    // --------------------------------------------------------------
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
        const ALocal_EBC * const &ebc_part );


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
        const ALocal_EBC * const &ebc_part );


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
        const ALocal_EBC * const &ebc_part );

    virtual double Assem_outsurface_flowrate(
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &ele_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id ) const;

  private:
    int nLocBas, dof;

    int snLocBas; // surface_nLocBas

    int num_ebc; // number of ebc surface regions

    PetscInt * row_index;

    PetscInt * srow_index;

    // array_a/b stores the local processor's portion of the solution,
    // including the solution on the ghost nodes. It has length: dof x
    // nlocalghonode.
    double * array_a;
    double * array_b;

    // local_a/b stores the elemenet's nodes' solution values. It has
    // length: dof x nLocBas. 
    double * local_a;
    double * local_b;

    // local_a/b on surfaces. It has dof x surface_nLocBas
    double * local_as;
    double * local_bs;

    // IEN_e: the element's IEN array
    int * IEN_e;

    // LSIEN: local surface element's IEN array
    int * LSIEN;

    // ectrl_x/y/z: the element's nodal control points. They have length
    // nLocBas
    double * ectrl_x;
    double * ectrl_y;
    double * ectrl_z;

    // sctrl_x/y/z: the surface element's nodal control points. They have
    // length surface_nLocBas;
    double * sctrl_x;
    double * sctrl_y;
    double * sctrl_z;

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


    void GetLocal( const double * const &array, const int * const &IEN,
        const int &in_locbas, double * const &local_array) const
    {
      int offset1, offset2;
      for(int ii=0; ii<in_locbas; ++ii)
      {
        offset1 = ii * dof;
        offset2 = IEN[ii] * dof;
        for(int jj=0; jj<dof; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }


    void NatBC_G( const double &curr_time, const double &dt,
        IPLocAssem * const &lassem_ptr,
        FEAElement * const &element_s,
        const int &in_loc_dof,
        const IQuadPts * const &quad_s,
        const ALocal_IEN * const lien_ptr,
        const ALocal_NodalBC * const &nbc_part,
        const ALocal_EBC * const &ebc_part );


    // --------------------------------------------------------------
    // Get_dnz_onz: obtain the row-by-row estimate of the number of
    // nonzeros per row. This estimate should be accurate for regular,
    // dirichlet, and perioidc-slave nodes; it may overestimate the
    // number of nz for periodic-master nodes.
    //
    // Note: This estimate is "too" accurate for a working matrix. In
    // practice, we do not often need such accurate estimate since
    // there is alwasy a room left for preconditioners. So, a rough,
    // fast estimator will be sufficient for practical simulations.
    // --------------------------------------------------------------
    void Get_dnz_onz( const int &nElem,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const ALocal_NodalBC * const &nbc_part,
        PetscInt * const &dnz, PetscInt * const &onz ) const;


    // --------------------------------------------------------------
    // Get_dnz_onz: obtain the row-by-row estimate of the number of 
    // nonzeros per row by empirical estimates. This estimate will
    // be loose, since the preconditioner will take a significant 
    // amount of memeory.
    // For structural NURBS mesh, empirical_neighbor_node_number =
    // (2*sdeg+1) * (2*tdeg+1) * (2*udeg+1)
    // For linear tetrahedral mesh, we assign empirical_neighbor_node_
    // number = 25.
    // --------------------------------------------------------------
    void Get_dnz_onz( const int &nlocnode,
        const int &empirical_neighbor_node_number,
        const ALocal_NodalBC * const &nbc_ptr,
        PetscInt * const &dnz, PetscInt * const &onz ) const;


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
};


#endif
