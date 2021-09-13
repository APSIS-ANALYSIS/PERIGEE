#include "MaterialModel_GOH06_Incompressible_Mixed.hpp"

MaterialModel_GOH06_Incompressible_Mixed::MaterialModel_GOH06_Incompressible_Mixed(
    const double &in_rho, const double &in_mu, 
    const double &in_f1the, const double &in_f1phi,
    const double &in_f2the, const double &in_f2phi,
    const double &in_fk1, const double &in_fk2, 
    const double &in_fkd )
: pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ), pi( MATH_T::PI ),
  rho0( in_rho ), E(3.0 * in_mu), nu(0.5), mu( in_mu ),
  f1_the( in_f1the*pi/180.0 ), f1_phi( in_f1phi*pi/180.0 ),
  f2_the( in_f2the*pi/180.0 ), f2_phi( in_f2phi*pi/180.0 ),
  fk1(in_fk1), fk2(in_fk2), fkd(in_fkd),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{
  a1[0] = sin(f1_the) * cos(f1_phi);
  a1[1] = sin(f1_the) * sin(f1_phi);
  a1[2] = cos(f1_the);

  a2[0] = sin(f2_the) * cos(f2_phi);
  a2[1] = sin(f2_the) * sin(f2_phi);
  a2[2] = cos(f2_the);

  a1xa1.gen_outprod(a1);
  a2xa2.gen_outprod(a2);
}


MaterialModel_GOH06_Incompressible_Mixed::~MaterialModel_GOH06_Incompressible_Mixed()
{}


void MaterialModel_GOH06_Incompressible_Mixed::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  MaterialModel_GOH06_Incompressible_Mixed: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Density rho  = %e \n", rho0);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Ground Matrix Neo-Hookean: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Fibre Fung: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Angle theta_1 (deg)= %e \n", f1_the*180/pi);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Angle phi_1 (deg)  = %e \n", f1_phi*180/pi);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Angle theta_1 (rad)= %e \n", f1_the);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Angle phi_1 (rad)  = %e \n", f1_phi);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Angle theta_2 (deg)= %e \n", f2_the*180/pi);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Angle phi_2 (deg)  = %e \n", f2_phi*180/pi);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Angle theta_2 (rad)= %e \n", f2_the);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Angle phi_2 (rad)  = %e \n", f2_phi);
  PetscPrintf(PETSC_COMM_WORLD, "\t  a1: [ %e , %e , %e ] \n", a1[0], a1[1], a1[2]);
  PetscPrintf(PETSC_COMM_WORLD, "\t  a2: [ %e , %e , %e ] \n", a2[0], a2[1], a2[2]);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Fibre k1   = %e \n", fk1);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Fibre k2   = %e \n", fk2);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Fibre k_dispersion = %e \n", fkd);
}

void MaterialModel_GOH06_Incompressible_Mixed::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar( "rho0", rho0 );
    h5w -> write_doubleScalar("E", E);
    h5w -> write_doubleScalar("nu", nu);
    h5w -> write_doubleScalar("lambda", lambda);
    h5w -> write_doubleScalar("mu", mu);
    h5w -> write_doubleScalar("kappa", kappa);
    h5w -> write_doubleScalar("f1_the", f1_the);
    h5w -> write_doubleScalar("f1_phi", f1_phi);
    h5w -> write_doubleScalar("f2_the", f2_the);
    h5w -> write_doubleScalar("f2_phi", f2_phi);
    h5w -> write_doubleScalar("fk1", fk1);
    h5w -> write_doubleScalar("fk2", fk2);
    h5w -> write_doubleScalar("k_dispersion" fkd);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

void MaterialModel_GOH06_Incompressible_Mixed::get_PK( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S)
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();
  trC = C.tr();
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  fE1 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a1Ca1) - 1.0;
  fE2 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a2Ca2) - 1.0;
  
  const double dfpsi1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1 );
  const double dfpsi2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2 );

  // P : I = I - 1/3 trC C^-1
  PxI.copy(Cinv); PxI.scale( (-1.0) * pt33 * trC );
  PxI.PY(I);

  // P : H = kd PxI + (1-3kd) a x a + (kd - 1/3) (a.Ca) C^-1
  PxH1.gen_zero(); PxH2.gen_zero();
  PxH1.AXPY( fkd, PxI );
  PxH2.AXPY( fkd, PxI );

  PxH1.AXPY( 1.0-3.0*fkd, a1xa1 );
  PxH2.AXPY( 1.0-3.0*fkd, a2xa2 );
  
  PxH1.AXPY( (fkd-pt33)*a1Ca1, Cinv );
  PxH2.AXPY( (fkd-pt33)*a2Ca2, Cinv );

  // S = mu J^-2/3 P:I + 2J^-2/3 dpsi_i P : H_i
  S.gen_zero();
  S.AXPY( mu * detFm0d67, PxI);
  
  S.AXPY(2.0 * detFm0d67 * dfpsi1, PxH1);
  S.AXPY(2.0 * detFm0d67 * dfpsi2, PxH2);

  // 1st PK
  P.MatMult(F,S);
}


void MaterialModel_GOH06_Incompressible_Mixed::get_PK_Stiffness(
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &CC )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();
  trC = C.tr();
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  fE1 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a1Ca1) - 1.0;
  fE2 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a2Ca2) - 1.0;
  
  const double dfpsi1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1 );
  const double dfpsi2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2 );

  // P : I = I - 1/3 trC C^-1
  PxI.copy(Cinv); PxI.scale( (-1.0) * pt33 * trC );
  PxI.PY(I);

  // P : H = kd PxI + (1-3kd) a x a + (kd - 1/3) (a.Ca) C^-1
  PxH1.gen_zero(); PxH2.gen_zero();
  PxH1.AXPY( fkd, PxI );
  PxH2.AXPY( fkd, PxI );

  PxH1.AXPY( 1.0-3.0*fkd, a1xa1 );
  PxH2.AXPY( 1.0-3.0*fkd, a2xa2 );
  
  PxH1.AXPY( (fkd-pt33)*a1Ca1, Cinv );
  PxH2.AXPY( (fkd-pt33)*a2Ca2, Cinv );

  // S = mu J^-2/3 P:I + 2J^-2/3 dpsi_i P : H_i
  S.gen_zero();
  S.AXPY( mu * detFm0d67, PxI);
  
  S.AXPY(2.0 * detFm0d67 * dfpsi1, PxH1);
  S.AXPY(2.0 * detFm0d67 * dfpsi2, PxH2);
  
  P.MatMult(F,S);

  // Assembly the stiffness matrix
  const double d2fpsi1 = fk1 * (1.0 + 2.0*fk2*fE1*fE1) * std::exp(fk2*fE1*fE1); 
  const double d2fpsi2 = fk1 * (1.0 + 2.0*fk2*fE2*fE2) * std::exp(fk2*fE2*fE2); 

  CC.gen_zero();

  const double val1 = 4.0 * detFm0d67 * detFm0d67;

  CC.add_OutProduct(val1 * d2fpsi1, PxH1, PxH1);
  CC.add_OutProduct(val1 * d2fpsi2, PxH2, PxH2);

  const double val2 = 2.0 * pt33 * detFm0d67 * ( mu * trC 
      + 2.0 * dfpsi1 * ( fkd * trC + (1.0-3.0*fkd)*a1Ca1 )
      + 2.0 * dfpsi2 * ( fkd * trC + (1.0-3.0*fkd)*a2Ca2 ) );

  CC.add_SymmProduct(val2, Cinv, Cinv);
  CC.add_OutProduct((-1.0)*pt33*val2, Cinv, Cinv);
  
  CC.add_OutProduct(mpt67, Cinv, S);
  CC.add_OutProduct(mpt67, S, Cinv);
}



double MaterialModel_GOH06_Incompressible_Mixed::get_strain_energy(
    const Matrix_3x3 &F )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();
  trC = C.tr();
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  fE1 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a1Ca1) - 1.0;
  fE2 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a2Ca2) - 1.0;

  const double PSI_iso = 0.5 * mu * (detFm0d67 * trC - 3.0);
  const double PSI_fi1 = 0.5 * (fk1 / fk2) * ( std::exp(fk2*fE1*fE1) - 1.0 );
  const double PSI_fi2 = 0.5 * (fk1 / fk2) * ( std::exp(fk2*fE2*fE2) - 1.0 );

  return PSI_iso + PSI_fi1 + PSI_fi2;
}


void MaterialModel_GOH06_Incompressible_Mixed::get_fibre_dir( const int &dir,
    double &fa1, double &fa2, double &fa3 ) const
{
  if(dir == 0)
  {
    fa1 = a1[0];
    fa2 = a1[1];
    fa3 = a1[2];
  }
  else if(dir == 1)
  {
    fa1 = a2[0];
    fa2 = a2[1];
    fa3 = a2[2];
  }
  else
  {
    SYS_T::print_fatal("Error: MaterialModel_GOH06_Incompressible_Mixed, wrong fibre direction. \n");
  }
}


// EOF
