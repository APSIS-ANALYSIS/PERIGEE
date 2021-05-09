#ifndef MATRIX_PETSC_CMM_HPP
#define MATRIX_PETSC_CMM_HPP
// ============================================================================
// Matrix_PETSc_CMM.hpp
//
// This is the Matrix that we use to enforce the essential BC after within the
// multi-corrector of nonlinear solvers.
//
// This particular derived matrix object is designed for CMM type problems, in
// which there are 4 dofs defined for each node and the 1st dof is pressure.
//
// Author: Ju Liu
// Date Created: April 28 2021
// ============================================================================
#include "Matrix_PETSc.hpp"
#include "ALocal_Ring_NodalBC.hpp"

class Matrix_PETSc_CMM : public Matrix_PETSc
{
  public:
    Matrix_PETSc_CMM( const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &bc_part,
        const ALocal_Ring_NodalBC * const &ring_bc_part );

    virtual ~Matrix_PETSc_CMM();

  private:
    // ------------------------------------------------------------------------
    // Type 1 : Generate a matrix accounting for the essential
    // boundary conditions that describe the in plane motion for a 4-dof system
    // like the NS equations.
    // Assumption: the matrix system has 4 degrees-of-freedom per node and the
    //             1st dof is pressure, the next 3 dofs are velocity.
    // 1. For ring nodes that belong to the ALocal_Ring_NodalBC class, the
    // dominant component's row will be modified. Let the dominant component's
    // corresponding entry in the outward normal be Nd, and the remaining two
    // non-dominant entries be Na, Nb. Then in the dominant component's row, the
    // diagonal entry will be 0, and entries in the non-dominant components'
    // columns will be -Na / Nd and - Nb / Nd. For the two non-dominant
    // components' rows, diagonal entries will be 1.
    // 2. For the remaining essential BC nodes, we assign 0 to all entries in all
    // components' rows.
    // 3. For all remaining nodes, we assign 1 to the diagonal entries.
    // ------------------------------------------------------------------------
    void gen_ring_inplane_bc( const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &bc_part,
        const ALocal_Ring_NodalBC * const &ring_bc_part );

    // ------------------------------------------------------------------------
    // Type 2: Generate a matrix accounting for the essential boundary
    // conditions that enforce purely radial motion for a 4-dof system.
    // Assumption: the matrix system has 4 degrees-of-freedom per node and the
    //             1st dof is pressure, the next 3 dofs are velocity.
    // The constraint equations are
    //             nx ux + ny uy + nz uz = 0
    //             tx ux + ty uy + tz uz = 0
    // Let a and b denote the dominant comp of n and t vectors respectively.
    //             na ua + nb ub = -nc uc
    //             ta ua + tb ub = -tc uc
    // Solve that equation gives
    //             ua = val1 uc
    //             ub = val2 uc
    //             val1 = (nb tc - nc tb) / (na tb - nb ta)
    //             val2 = (ta nc - na tc) / (na tb - nb ta)         
    // 1. For ring nodes that belong to the ALocal_Ring_NodalBC class, the rows
    //    for ua and ub with columns at uc are inserted with value val1 & val2,
    //    respectively.
    // 2. For the remaining essential BC nodes, we assign 0 to all entries in
    //    all components' rows.
    // 3. For all remaining nodes, we assign 1 to the diagonal entries. 
    // ------------------------------------------------------------------------
    void gen_ring_radial_motion_bc( const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &bc_part,
        const ALocal_Ring_NodalBC * const &ring_bc_part );

    // ------------------------------------------------------------------------
    // Type 3 : Generate a matrix accounting for the essential
    // boundary conditions that describe inlet clamping and outlet in plane motion
    // for a 4-dof system like the NS equations.
    // Assumption: the matrix system has 4 degrees-of-freedom per node and the
    //             1st dof is pressure, the next 3 dofs are velocity.
    // 1. For outlet ring nodes that belong to the ALocal_Ring_NodalBC class, the
    // dominant component's row will be modified. Let the dominant component's
    // corresponding entry in the outward normal be Nd, and the remaining two
    // non-dominant entries be Na, Nb. Then in the dominant component's row, the
    // diagonal entry will be 0, and entries in the non-dominant components'
    // columns will be -Na / Nd and - Nb / Nd. For the two non-dominant
    // components' rows, diagonal entries will be 1.
    // 2. For the remaining essential BC nodes, we assign 0 to all entries in all
    // components' rows.
    // 3. For all remaining nodes, we assign 1 to the diagonal entries.
    // ------------------------------------------------------------------------
    void gen_outlet_ring_inplane_bc( const APart_Node * const &pnode_ptr,
        const ALocal_NodalBC * const &bc_part,
        const ALocal_Ring_NodalBC * const &ring_bc_part );

};

#endif
