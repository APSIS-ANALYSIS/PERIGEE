#ifndef MATRIX_PETSC_HPP
#define MATRIX_PETSC_HPP
// ============================================================================
// Matrix_PETSc.hpp
// This is a class that implements the assembly of PETSc sparse matrices. The
// main purpose is to implement sparse matrix that serve simple matrix
// manipulations or preconditioners.
//
// Date: Feb. 14 2016
// Author: Ju Liu
// ============================================================================
#include "PETSc_Tools.hpp"
#include "ALocal_NBC.hpp"
#include "PDNSolution.hpp"
#include "petscmat.h"

class Matrix_PETSc
{
  public:
    // ------------------------------------------------------------------------
    // Constructor of PETSc Matrix object with dimension: loc_row x loc_col
    // in the local cpu portion. The diangonal and off-diagonal nonzero entries
    // are assumbed to be 1 respectively. The purpose is to create a sparse
    // permutation matrix.
    // ------------------------------------------------------------------------
    Matrix_PETSc( const int &loc_row, const int &loc_col, const int &dnz = 1,
       const int &onz = 1 );

    // loc_col = loc_row
    Matrix_PETSc( const int &loc_row, const int &dnz = 1, const int &onz = 1 );

    // ------------------------------------------------------------------------
    // Constructor: Generate a sparse square matrix with size defined
    // locally. Local row/column num = pnode_ptr->get_nlocalnode() 
    // * pnode_ptr->get_dof();
    // ------------------------------------------------------------------------
    Matrix_PETSc( const APart_Node * const &pnode_ptr, const int &dnz = 1,
       const int &onz = 1 );

    // ------------------------------------------------------------------------
    // Constructor: Generate a sparse square matrix with size defined
    // locally. Local row/column num = pnode_ptr->get_nlocalnode()
    // * bc_part->get_dof_LID();
    // ------------------------------------------------------------------------
    Matrix_PETSc( const APart_Node * const &pnode_ptr,
       const ALocal_NBC * const &bc_part,
       const int &dnz = 1, const int &onz = 1 );

    // ------------------------------------------------------------------------
    // Destroyer
    // ------------------------------------------------------------------------
    virtual ~Matrix_PETSc();

    // ------------------------------------------------------------------------
    // ! Flag : Fix nonzero structure
    //          Adding or inserting in new locations will be ignored. 
    //          Set after a correct nonzero structure has been assembled.
    // ------------------------------------------------------------------------
    virtual void Fix_nonzero_str()
    {MatSetOption(K, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);}

    // ------------------------------------------------------------------------
    // ! Flag : New allocation error
    //          Adding or inserting in new locations will generate an error.
    //          Supports AIJ and BAIJ formats.
    //          Typically used for debugging.
    // ------------------------------------------------------------------------
    virtual void Fix_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);}

    // ------------------------------------------------------------------------
    // ! Flag : Ignore new allocation
    //          Adding or inserting in new locations will NOT generate an error.
    //          Typically used if the estimate of the nonzero structure is an
    //          underestimate. This flag is called for a first practical
    //          assembly.
    // ------------------------------------------------------------------------
    virtual void Release_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);}

    // ------------------------------------------------------------------------
    // ! Flag : Keep nonzero pattern of the matrix K
    // ------------------------------------------------------------------------
    virtual void Keep_nonzero_pattern() 
    {MatSetOption(K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);}

    // ------------------------------------------------------------------------
    // MatView : print the matrix on screen
    // ------------------------------------------------------------------------
    virtual void print_info() const {MatView(K, PETSC_VIEWER_STDOUT_WORLD);}

    // ------------------------------------------------------------------------
    // Clear K
    // ------------------------------------------------------------------------
    virtual void Clear() {MatZeroEntries(K);}

    // ------------------------------------------------------------------------
    // gen_id : Generate an identity matrix that is compatible with the
    // mesh partition.
    // ------------------------------------------------------------------------
    virtual void gen_id( const APart_Node * const &pnode_ptr );

    // ------------------------------------------------------------------------
    // gen_perm_bc : Generate a permutation matrix accounting for the essential
    // boundary conditions.
    // ------------------------------------------------------------------------
    virtual void gen_perm_bc( const APart_Node * const &pnode_ptr,
        const ALocal_NBC * const &bc_part );

    // ------------------------------------------------------------------------
    // gen_perm_bc : Generate a permutation matrix accounting for the essential
    // boundary conditions, following the rule of arranging velocity first, then pressure.
    // ------------------------------------------------------------------------
    virtual void gen_perm_bc_block( const APart_Node * const &pnode_ptr,
      const ALocal_NBC * const &bc_part );

    // ------------------------------------------------------------------------
    // gen_perm_bc : Generate a permutation matrix for essential boundary
    // conditions for Multi-Field problems.
    // NOTE: LID from the ALocal_NBC here should give the matrix row/col id
    // ------------------------------------------------------------------------
    virtual void gen_perm_bc( const std::vector<APart_Node *> &pnode_list,
        const std::vector<ALocal_NBC *> &bc_part_list,
        const std::vector<int> &start_idx );

    // ------------------------------------------------------------------------
    // Gen_extractor_for_Dirichlet_nodes : Generate a matrix that is zero for
    // all regular rows. On each row corresponding to a Dirichlet dof, set
    // the diagonal entry in that row to 1. When applied to the
    // initial solution vector with a correct Dirichlet description, this matrix
    // will extract the Dirichlet nodal values; the remaining values are zero.
    // 
    //           [ 1 0 0       [ 3        [ 3
    //             0 0 0    x    2    =     0
    //             0 0 0 ]       1 ]        0 ]
    // 
    // ------------------------------------------------------------------------
    virtual void gen_extractor_for_Dirichlet_nodes( const APart_Node * const &pnode_ptr,
        const ALocal_NBC * const &bc_part );

    // ------------------------------------------------------------------------
    // MatMultSol : perform a matrix-vector multiplication : sol = K sol
    //              If the matrix partition is incompatible with the vector's,
    //              a PETSc error from PETSc::MatMult will be thrown.
    // ------------------------------------------------------------------------
    virtual void MatMultSol( PDNSolution * const &sol ) const;
    
    virtual void MatMultSol( Vec &sol ) const;

  protected:
    Mat K;

    // Global dimension of K matrix
    int gm, gn;

    // Local dimension of K matrix
    int lm, ln;

    // flag indicating whether the matrix is assembled 
    bool is_set;
};

#endif
