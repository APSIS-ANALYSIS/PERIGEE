#ifndef PGASSEM_SEGREGATED_HPP
#define PGASSEM_SEGREGATED_HPP
// ============================================================================
// PGAssem_Segregated.hpp
// Parallel global assembly routine based on PETSc-3.6, using AIJ matrix format.
// The quadrature info for finite element basis funcitons are computed within
// the local assembly code.
//
// This assembly code is designed to release the fully implicit solver
// assumption. The local_array will be allowed to have more flexible memory size
// so that some parts of the equation can be solved explicitly, or in a
// segregated fashion.
//
// This assembly code is designed for FSI or segregated solids solver.
//
// Date: Dec. 10 2016
// Author: Ju Liu
// ============================================================================

#include <ctime>
#include "IPGAssem.hpp"

class PGAssem_Segregated : public IPGAssem
{
  public:
    PGAssem_Segregated(
        IPLocAssem const * const &locassem_ptr,
        IAGlobal_Mesh_Info const * const &agmi_ptr,
        ALocal_Elem const * const &alelem_ptr,
        ALocal_IEN const * const &aien_ptr,
        APart_Node const * const &pnode_ptr,
        ALocal_NodalBC const * const &part_nbc );

    virtual ~PGAssem_Segregated();

  private:
    int nLocBas, dof_a, dof_b;

    PetscInt * row_index;

    double * array_a, * array_da, * array_b, * array_db;

    double * local_a, * local_da, * local_b, * local_db;

    int * IEN_e;
};


#endif
