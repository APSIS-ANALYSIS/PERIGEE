#include "SymmTensor4_3D.hpp"

SymmTensor4_3D::SymmTensor4_3D()
{
  gen_zero();

}

SymmTensor4_3D::~SymmTensor4_3D()
{}

void SymmTensor4_3D::gen_rand()
{}

void SymmTensor4_3D::gen_zero()
{
  ten = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0};
}

// EOF
