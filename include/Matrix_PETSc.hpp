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
#include "ALocal_NodalBC.hpp"
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
    Matrix_PETSc( const int &loc_row, const int &loc_col );

    
    // ------------------------------------------------------------------------
    // Constructor: Generate a sparse square matrix with size defined
    // locally. Local row/column num = pnode_ptr->get_nlocalnode() 
    // * pnode_ptr->get_dof();
    // ------------------------------------------------------------------------
    Matrix_PETSc( const APart_Node * const &pnode_ptr );

    
    // ------------------------------------------------------------------------
    // Constructor: Generate a sparse square matrix with size defined
    // locally. Local row/column num = pnode_ptr->get_nlocalnode()
    // * bc_part->get_dofMat();
    // ------------------------------------------------------------------------
    Matrix_PETSc( const APart_Node * const &pnode_ptr,
       const ALocal_NodalBC * const &bc_part );


    // ------------------------------------------------------------------------
    // Destroyer
    // ------------------------------------------------------------------------
    virtual ~Matrix_PETSc();


    // ------------------------------------------------------------------------
    // ! Flag : Fix nonzero structure
    //          Add or insert in new location is ignored. 
    //          Set after a correct nonzero structure assembled.
    // ------------------------------------------------------------------------
    void Fix_nonzero_str()
    {MatSetOption(K, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);}


    // ------------------------------------------------------------------------
    // ! Flag : New allocation error
    //          Add or insert in new locations will generate an error message.
    //          Supports AIJ and BAIJ formats.
    //          Typically used for debugging.
    // ------------------------------------------------------------------------
    void Fix_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);}


    // ------------------------------------------------------------------------
    // ! Flag : Ignore new allocation
    //          Add or insert in a new allocation will NOT generate error.
    //          Typically used if the estimate of nonzero structure
    //          underestimates, this flag is called for a first practical
    //          assembly.
    // ------------------------------------------------------------------------
    void Release_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);}


    // ------------------------------------------------------------------------
    // ! Flag : Keep nonzero pattern of the matrix K
    // ------------------------------------------------------------------------
    void Keep_nonzero_pattern() 
    {MatSetOption(K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);}


    // ------------------------------------------------------------------------
    // MatView : print the matrix on screen
    // ------------------------------------------------------------------------
    void print_matrix() const {MatView(K, PETSC_VIEWER_STDOUT_WORLD);}


    // ------------------------------------------------------------------------
    // Clear K
    // ------------------------------------------------------------------------
    void Clear() {MatZeroEntries(K);}


    // ------------------------------------------------------------------------
    // Gen_id : Generate an idnetity matrix that is compatible with the
    // mesh partition.
    // ------------------------------------------------------------------------
    void gen_id( const APart_Node * const &pnode_ptr );


    // ------------------------------------------------------------------------
    // Gen_perm_bc : Generate a permutation matrix accounting for the essential
    // boundary conditions.
    // ------------------------------------------------------------------------
    void gen_perm_bc( const APart_Node * const &pnode_ptr , 
       const ALocal_NodalBC * const &bc_part );

    
    // ------------------------------------------------------------------------
    // Gen_extractor_for_Dirichlet_nodes : Generate a matrix that is zero for
    // all regular rows. On the row corresponding to the Dirichlet dof, make
    // the diagonal entry in that row to be 1. Applying this matrix to the
    // initial solution vector with correct Dirichlet description, this matrix
    // will extract the Dirichlet nodal values; the rest values are zero.
    // 
    //           [ 1 0 0       [ 3        [ 3
    //             0 0 0    x    2    =     0
    //             0 0 0 ]       1 ]        0 ]
    // 
    // ------------------------------------------------------------------------
    void gen_extractor_for_Dirichlet_nodes( const APart_Node * const &pnode_ptr ,
        const ALocal_NodalBC * const &bc_part );


    // ------------------------------------------------------------------------
    // MatMultSol : perform a matrix-vector multiplication : sol = K sol
    //              If the matrix partition is incompatible with the vector's,
    //              a PETSc error from PETSc::MatMult will throw.
    // ------------------------------------------------------------------------
    void MatMultSol( PDNSolution * const &sol ) const;

  private:
    Mat K;

    // Global dimension of K matrix
    int m, n;

    // Local dimension of K matrix
    int lm, ln;

    // flag indicating if the matrix is assembled    
    bool is_set;
};

#endif
