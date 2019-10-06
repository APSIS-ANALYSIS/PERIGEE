#ifndef PGASSEM_V360_NOCACHE_NURBS_HPP
#define PGASSEM_V360_NOCACHE_NURBS_HPP
// ============================================================================
// PGAssem_v360_noCache_NURBS.hpp
// Parallel global assembly based on PETSc-3.6.0, with no pre-cached quadrature
// info for NURBS basis functions.
// 
// The Dirichlet boundary condition is STRONGLY enforced;
// The Periodic boundary condition is STRONGLY enforced;
// The Neumann boundary condition is enforced weakly.
//
// Date: July 31 2015
// ============================================================================
#include "IPGAssem.hpp"

class PGAssem_v360_noCache_NURBS : public IPGAssem
{
  public:
    // ------------------------------------------------------------------------
    // ! Constructor
    // ------------------------------------------------------------------------
    PGAssem_v360_noCache_NURBS( 
        const IPLocAssem * const &locassem_ptr,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const IALocal_BC * const &part_bc );

    // ------------------------------------------------------------------------
    // ! Destructor
    // ------------------------------------------------------------------------
    virtual ~PGAssem_v360_noCache_NURBS();


    // ------------------------------------------------------------------------
    // ! Assem_nonzero_estimate
    //   Assembly the nonzero estimated matrix for K (add 1.0 to every possible
    //   slot.)
    //   \para ALocal_Elem
    //   \para IPLocAssem
    //   \para ALocal_IEN
    //   \para APart_Node
    //   \para IALocal_BC
    // ------------------------------------------------------------------------
    virtual void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const IALocal_BC * const &bc_part );



    // ------------------------------------------------------------------------
    // ! Assembly mass matrix and its residual vector for 3D problems
    //   without cached quadrature info.
    //    -------------------
    //   The following input are needed for element quadrature evaluation
    //   \para IALocal_meshSize : object saving the mesh size
    //   \para BernsteinBasis : pre-evaluated Bernstein polynomial at Gauss 
    //                          quadrature points for volumetric integration
    //   \para IAExtractor : Bezier extraction operator
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
        const IALocal_BC * const &bc_part );



    // ------------------------------------------------------------------------
    // ! Assem residual vector for 3D problem without cached quadrature
    //   info.
    //   -----------------
    //   The following input are needed for element quadrature evaluation
    //   IALocal_meshSize : object saving the mesh size
    //   BernsteinBasis : pre-evaluated Bernstein polynomial at Gauss 
    //                    quadrature points for volumetric integration
    //   IAExtractor : Bezier extraction operator
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
        const IALocal_BC * const &bc_part );


    // ------------------------------------------------------------------------
    // ! Assembly tangent matrix and residual vector together.
    //     The quadrature info is computed within the assembly routine.
    //     This function is for 3D problems
    //     sol_a : disp or sol_{n+1}
    //     sol_b : velo or sol_{n}
    //   -----------------
    //   The following input are needed for element quadrature evaluation
    //     IALocal_meshSize : object saving the mesh size
    //     BernsteinBasis : pre-evaluated Bernstein polynomial at Gauss 
    //                      quadrature points for volumetric integration
    //     IAExtractor : Bezier extraction operator
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
        const IALocal_BC * const &bc_part );



  private:
    // ------------------------------------------------------------------------
    // Disallow constructor with no parameter
    // ------------------------------------------------------------------------
    PGAssem_v360_noCache_NURBS()
      : nLocBas(0), dof(0)
    {};

    // ------------------------------------------------------------------------
    // ! Initialize MPIAIJ matrix using PETSc-3.6.0
    //   \para dnz : number of nonzero in diagonal portion for each row
    //   \para onz : number of nonzero in off-diagonal portion for each row
    //   \para num_loc_row : number of rows for this subdomain (partition)
    // ------------------------------------------------------------------------
    void Init_petsc_360(const PetscInt * const &dnz,
        const PetscInt * const &onz, const int &num_loc_row );


    // ------------------------------------------------------------------------
    // ! Get_dnz_onz
    //   This function gives an estimate for the number of the diagonal
    //   nonzeros and the off-diagonal nonzeros. The results are
    //   returned as two arrays with length dof * nlocalnode.
    //   This function is used to give an accurate matrix preallocation for 
    //   the PETSc matrix Mat.
    //   Output : 
    //   \para dnz : the number of nonzeros in diagonal portion
    //   \para onz : the number of nonzeros in off-diagonal portion
    // ------------------------------------------------------------------------
    void Get_dnz_onz( const int &nElem,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const IALocal_BC * const &bc_part,
        PetscInt * const &dnz, PetscInt * const &onz ) const;


    // ------------------------------------------------------------------------
    // ! Get_part_id
    //   Get the partition id for given node index from the given list of nodal
    //   partition list. nlist has to start as nlist[0] = 0;
    //   nlist = {0, 3, 5, 8} means that partition assigns 3 cpus.
    //   cpu 0 owns node 0 to node 2.
    //   cpu 1 owns node 3 to node 4.
    //   cpu 2 owns node 5 to node 8.
    //   if val = 1, returns 0.
    //   if val > *(nlist.end()), return -1.
    // -------------------------------------------------------------------------
    int Get_part_id( const std::vector<unsigned int> &nlist, const unsigned int &val ) const
    {
      for(unsigned int ii=1; ii<nlist.size(); ++ii)
      {
        if(val < nlist[ii])
          return int(ii-1);
      }
      return -1;
    }


    // ------------------------------------------------------------------------
    // ! Get local element's array.
    //   the long array should be listed compatible with the LIEN
    //   structure, for example, local_to_global() ~ LIEN[][].
    //   user is responsible for creating and deleting local_array pointer
    // ------------------------------------------------------------------------
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


    // ------------------------------------------------------------------------
    // ! EssBC_KG : Essential boundary condition enforcement for Mat K and Vec
    //              G. This function ADD_VALUES for the matrix K, which
    //              implicityly implies that the previous assembly should use
    //              LID in global assembly. Only one MatAssemblyBegin/End needs
    //              to be called.
    // ------------------------------------------------------------------------
    void EssBC_KG( const IALocal_BC * const &bc_part, const int &field );


    // ------------------------------------------------------------------------
    // ! EssBC_G : update the residual vector for essential boundary conditions.
    // ------------------------------------------------------------------------
    void EssBC_G( const IALocal_BC * const &bc_part, const int &field );


    // ------------------------------------------------------------------------
    // ! PRIVATE DATA STRUCTURES
    // ------------------------------------------------------------------------
    // nLocBas : number of basis functions in one NURBS element
    const int nLocBas;

    // dof : degree of freedom of this problem
    const int dof;

    // row_index & col_index : cache the element's LID indicies in assembly
    PetscInt * row_index;
    //PetscInt * col_index;

    // array_a & array_b : store the solution vector in this subdomain
    double * array_a;
    double * array_b;

    // local_a & local_b : store the solution vector in one element
    double * local_a;
    double * local_b;

    // IEN_e : LIEN indices in one element
    int * IEN_e;

    // ectrl_x/y/z/w : control points xyz-coordinates and weights in one element
    double * ectrl_x;
    double * ectrl_y;
    double * ectrl_z;
    double * ectrl_w;

};


#endif
