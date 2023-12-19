#include "MaterialModel_GOH14_ST91_Mixed.hpp"

MaterialModel_GOH14_ST91_Mixed::MaterialModel_GOH14_ST91_Mixed(
    const double &in_rho, const double &in_E, const double &in_nu,
    const double &in_f1the, const double &in_f1phi,
    const double &in_f2the, const double &in_f2phi,
    const double &in_fk1, const double &in_fk2,
    const double &in_fkd )
: rho0( in_rho ), E(in_E), nu(in_nu), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 ),
  f1_the( in_f1the*MATH_T::PI/180.0 ), f1_phi( in_f1phi*MATH_T::PI/180.0 ),
  f2_the( in_f2the*MATH_T::PI/180.0 ), f2_phi( in_f2phi*MATH_T::PI/180.0 ),
  fk1(in_fk1), fk2(in_fk2), fkd(in_fkd)
{
  a1(0) = sin(f1_the) * cos(f1_phi);
  a1(1) = sin(f1_the) * sin(f1_phi);
  a1(2) = cos(f1_the);

  a2(0) = sin(f2_the) * cos(f2_phi);
  a2(1) = sin(f2_the) * sin(f2_phi);
  a2(2) = cos(f2_the);
}

MaterialModel_GOH14_ST91_Mixed::MaterialModel_GOH14_ST91_Mixed(
    const char * const &fname )
{
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_GOH06_ST91_Mixed constructor does not match h5 file.\n" );

  rho0   = h5r -> read_doubleScalar("/", "rho0");
  E      = h5r -> read_doubleScalar("/", "E");
  nu     = h5r -> read_doubleScalar("/", "nu");
  lambda = h5r -> read_doubleScalar("/", "lambda");
  mu     = h5r -> read_doubleScalar("/", "mu");
  kappa  = h5r -> read_doubleScalar("/", "kappa");
  f1_the = h5r -> read_doubleScalar("/", "f1_the");
  f1_phi = h5r -> read_doubleScalar("/", "f1_phi");
  f2_the = h5r -> read_doubleScalar("/", "f2_the");
  f2_phi = h5r -> read_doubleScalar("/", "f2_phi");
  fk1    = h5r -> read_doubleScalar("/", "fk1");
  fk2    = h5r -> read_doubleScalar("/", "fk2");
  fkd    = h5r -> read_doubleScalar("/", "fkd");

  delete h5r; H5Fclose(h5file);

  a1(0) = sin(f1_the) * cos(f1_phi);
  a1(1) = sin(f1_the) * sin(f1_phi);
  a1(2) = cos(f1_the);

  a2(0) = sin(f2_the) * cos(f2_phi);
  a2(1) = sin(f2_the) * sin(f2_phi);
  a2(2) = cos(f2_the);
}

void MaterialModel_GOH14_ST91_Mixed::print_info() const
{
  SYS_T::commPrint("\t  MaterialModel_GOH14_ST91_Mixed: \n");
  SYS_T::commPrint("\t  Ground Matrix Neo-Hookean: \n");
  SYS_T::commPrint("\t  Young's Modulus E  = %e \n", E);
  SYS_T::commPrint("\t  Possion's ratio nu = %e \n", nu);
  SYS_T::commPrint("\t  Shear modulus mu   = %e \n", mu);
  SYS_T::commPrint("\t  Lame coeff lambda  = %e \n", lambda);
  SYS_T::commPrint("\t  Bulk modulus kappa = %e \n", kappa);
  SYS_T::commPrint("\t  Fibre Fung: \n");
  SYS_T::commPrint("\t  Angle theta_1 (deg)= %e \n", f1_the*180/MATH_T::PI);
  SYS_T::commPrint("\t  Angle phi_1 (deg)  = %e \n", f1_phi*180/MATH_T::PI);
  SYS_T::commPrint("\t  Angle theta_1 (rad)= %e \n", f1_the);
  SYS_T::commPrint("\t  Angle phi_1 (rad)  = %e \n", f1_phi);
  SYS_T::commPrint("\t  Angle theta_2 (deg)= %e \n", f2_the*180/MATH_T::PI);
  SYS_T::commPrint("\t  Angle phi_2 (deg)  = %e \n", f2_phi*180/MATH_T::PI);
  SYS_T::commPrint("\t  Angle theta_2 (rad)= %e \n", f2_the);
  SYS_T::commPrint("\t  Angle phi_2 (rad)  = %e \n", f2_phi);
  SYS_T::commPrint("\t  a1: [ %e , %e , %e ] \n", a1(0), a1(1), a1(2));
  SYS_T::commPrint("\t  a2: [ %e , %e , %e ] \n", a2(0), a2(1), a2(2));
  SYS_T::commPrint("\t  Fibre k1   = %e \n", fk1);
  SYS_T::commPrint("\t  Fibre k2   = %e \n", fk2);
  SYS_T::commPrint("\t  Fibre k_dispersion = %e \n", fkd);
}

void MaterialModel_GOH14_ST91_Mixed::write_hdf5( const char * const &fname ) const
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
    h5w -> write_doubleScalar("fkd", fkd);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

void MaterialModel_GOH14_ST91_Mixed::get_PK( 
    const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  const double trC = C.tr();
  const double detFm0d67 = std::pow(F.det(), - 2.0 / 3.0);

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  const double fE1 = fkd*trC + (1.0-3.0*fkd)*a1Ca1 - 1.0;
  const double fE2 = fkd*trC + (1.0-3.0*fkd)*a2Ca2 - 1.0;
  
  const double dfpsi1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1 );
  const double dfpsi2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2 );

  // P : I = I - 1/3 trC C^-1
  Tensor2_3D PxI(Cinv); PxI.scale( (-1.0) * trC / 3.0 );
  PxI.AXPI( 1.0 );

  Tensor2_3D PxH1; Tensor2_3D PxH2;
  PxH1.gen_zero(); PxH2.gen_zero();
  PxH1.AXPI( fkd );
  PxH2.AXPI( fkd );

  Tensor2_3D a1xa1, a2xa2;
  a1xa1.gen_outprod(a1);
  a2xa2.gen_outprod(a2);

  PxH1.AXPY( 1.0-3.0*fkd, a1xa1 );
  PxH2.AXPY( 1.0-3.0*fkd, a2xa2 );

  S.gen_zero();
  S.AXPY( mu * detFm0d67, PxI);
  
  S.AXPY(2.0 * dfpsi1, PxH1);
  S.AXPY(2.0 * dfpsi2, PxH2);

  // 1st PK
  P.MatMult(F,S);
}

void MaterialModel_GOH14_ST91_Mixed::get_PK_Stiffness( 
    const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S, Tensor4_3D &CC ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  const double trC = C.tr();
  const double detFm0d67 = std::pow(F.det(), - 2.0 / 3.0);
  const Tensor2_3D I = Ten2::gen_id();

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  const double fE1 = fkd*trC + (1.0-3.0*fkd)*a1Ca1 - 1.0;
  const double fE2 = fkd*trC + (1.0-3.0*fkd)*a2Ca2 - 1.0;
  
  const double dfpsi1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1 );
  const double dfpsi2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2 );

  // P : I = I - 1/3 trC C^-1
  Tensor2_3D PxI(Cinv); PxI.scale( (-1.0) * trC / 3.0 );
  PxI.AXPI( 1.0 );

  Tensor2_3D PxH1; Tensor2_3D PxH2;
  PxH1.gen_zero(); PxH2.gen_zero();
  PxH1.AXPI( fkd );
  PxH2.AXPI( fkd );

  Tensor2_3D a1xa1, a2xa2;
  a1xa1.gen_outprod(a1);
  a2xa2.gen_outprod(a2);

  PxH1.AXPY( 1.0-3.0*fkd, a1xa1 );
  PxH2.AXPY( 1.0-3.0*fkd, a2xa2 );
  
  S.gen_zero();
  S.AXPY( mu * detFm0d67, PxI);
  
  S.AXPY(2.0 * dfpsi1, PxH1);
  S.AXPY(2.0 * dfpsi2, PxH2);
  
  P.MatMult(F,S);

  // Assembly the stiffness matrix
  const double d2fpsi1 = fk1 * (1.0 + 2.0*fk2*fE1*fE1) * std::exp(fk2*fE1*fE1); 
  const double d2fpsi2 = fk1 * (1.0 + 2.0*fk2*fE2*fE2) * std::exp(fk2*fE2*fE2); 

  CC.gen_zero();

  CC.add_OutProduct(4.0 * d2fpsi1, PxH1, PxH1);
  CC.add_OutProduct(4.0 * d2fpsi2, PxH2, PxH2);

  const double val = - 2.0 * mu * detFm0d67 / 3.0;

  CC.add_OutProduct(val, PxI, Cinv);
  CC.add_OutProduct(val, Cinv, I);
  CC.add_SymmProduct(-val * trC, Cinv, Cinv);
}

double MaterialModel_GOH14_ST91_Mixed::get_strain_energy(const Tensor2_3D &F ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  const double trC = C.tr();
  const double detFm0d67 = std::pow(F.det(), - 2.0 / 3.0);

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  const double fE1 = fkd*trC + (1.0-3.0*fkd)*a1Ca1 - 1.0;
  const double fE2 = fkd*trC + (1.0-3.0*fkd)*a2Ca2 - 1.0;

  const double PSI_iso = 0.5 * mu * (detFm0d67 * trC - 3.0);
  const double PSI_fi1 = 0.5 * (fk1 / fk2) * ( std::exp(fk2*fE1*fE1) - 1.0 );
  const double PSI_fi2 = 0.5 * (fk1 / fk2) * ( std::exp(fk2*fE2*fE2) - 1.0 );

  return PSI_iso + PSI_fi1 + PSI_fi2;
}

Vector_3 MaterialModel_GOH14_ST91_Mixed::get_fibre_dir( const int &dir ) const
{
  if(dir == 0)
  {
    return a1;
  }
  else if(dir == 1)
  {
    return a2;
  }
  else
  {
    SYS_T::print_fatal("Error: MaterialModel_GOH14_ST91_Mixed, wrong fibre direction. \n");
    return Vector_3();
  }
}

double MaterialModel_GOH14_ST91_Mixed::get_rho( const double &p ) const
{
  const double pk = p / kappa;
  return rho0 * (std::pow(pk*pk+1.0, 0.5) + pk);
}

double MaterialModel_GOH14_ST91_Mixed::get_drho_dp( const double &p ) const
{
  const double pk = p / kappa;
  return (rho0 / kappa) * (1 + pk * std::pow(pk*pk+1.0, -0.5));
}

double MaterialModel_GOH14_ST91_Mixed::get_beta( const double &p ) const
{
  const double pk = p / kappa;
  return 1.0 / (kappa * std::pow(pk*pk+1.0, 0.5));
}

double MaterialModel_GOH14_ST91_Mixed::get_dbeta_dp( const double &p ) const
{
  const double pk = p / kappa;
  return -pk / ( kappa*kappa*std::pow(pk*pk+1.0, 1.5) );
}

// EOF
