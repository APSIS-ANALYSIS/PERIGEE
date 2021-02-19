// Test the new PDNSolution functions

#include "PDNSolution.hpp"

int main( int argc, char *argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const PetscMPIInt rank = SYS_T::get_MPI_rank();

  std::string part_file("part");

  APart_Node * pNode = new APart_Node(part_file, rank);

  PDNSolution * sol_a = new PDNSolution( pNode );
  
  delete sol_a;
  delete pNode; 

  PetscFinalize();
  return 0;
}

//EOF
