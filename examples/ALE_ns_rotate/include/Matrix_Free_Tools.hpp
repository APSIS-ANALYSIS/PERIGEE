#ifndef MATRIX_FREE_TOOLS_HPP
#define MATRIX_FREE_TOOLS_HPP

#include "PGAssem_NS_FEM.hpp"

namespace MF_T
{
  PetscErrorCode MF_MatMult(Mat shell, Vec X, Vec Y)
  {
    void * ptr;
    PGAssem_NS_FEM * user;
    MatShellGetContext(shell, &ptr);
    user = (PGAssem_NS_FEM*) ptr;

    MatMult(user->K, X, Y);

    user->Interface_K_MF(X, Y);

    return 0;
  }
}

#endif