#include "MaterialModel_StVenant_Kirchhoff.hpp"

MaterialModel_StVenant_Kirchhoff::MaterialModel_StVenant_Kirchhoff( 
    const double &in_E, const double &in_nu )
: E( in_E ), nu( in_nu ), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{}

MaterialModel_StVenant_Kirchhoff::~MaterialModel_StVenant_Kirchhoff()
{}


void MaterialModel_StVenant_Kirchhoff::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  MaterialModel_StVenant_Kirchhoff: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Bulk modulus kappa = %e \n", kappa);
}

void MaterialModel_StVenant_Kirchhoff_M94_Mixed::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar("E", E);
    h5w -> write_doubleScalar("nu", nu);
    h5w -> write_doubleScalar("lambda", lambda);
    h5w -> write_doubleScalar("mu", mu);
    h5w -> write_doubleScalar("kappa", kappa);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}


void MaterialModel_StVenant_Kirchhoff::get_PK( const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S )
{
  G.MatMultTransposeLeft(F); G.AXPY(-1.0, I); G.scale(0.5);
  S.gen_id(); S.scale(lambda * G.tr());
  S.AXPY(2.0 * mu, G); P.MatMult(F, S);
}


void MaterialModel_StVenant_Kirchhoff::get_PK_Stiffness( const Matrix_3x3 &F, Matrix_3x3 &P,
    Matrix_3x3 &S, Tensor4_3D &CC)
{
  G.MatMultTransposeLeft(F); G.AXPY(-1.0, I);  G.scale(0.5);
  S.gen_id(); S.scale(lambda * G.tr()); S.AXPY(2.0 * mu, G);
  P.MatMult(F, S);
  CC.gen_symm_id(); CC.scale(2.0 * mu);
  CC.add_OutProduct(lambda, I, I);
}


double MaterialModel_StVenant_Kirchhoff::get_strain_energy( const Matrix_3x3 &F )
{
  G.MatMultTransposeLeft(F); G.AXPY(-1.0, I); G.scale(0.5);
  const double trG = G.tr();
  G.MatMult(G,G);
  const double trG2 = G.tr();
  return 0.5 * lambda * trG + mu * trG2;
}


// EOF
