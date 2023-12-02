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
  fk1(in_fk1), fk2(in_fk2), fkd(in_fkd)
{
  a1(0) = sin(f1_the) * cos(f1_phi);
  a1(1) = sin(f1_the) * sin(f1_phi);
  a1(2) = cos(f1_the);

  a2(0) = sin(f2_the) * cos(f2_phi);
  a2(1) = sin(f2_the) * sin(f2_phi);
  a2(2) = cos(f2_the);

  a01(0) = sin(f1_the) * cos(f1_phi);
  a01(1) = sin(f1_the) * sin(f1_phi);
  a01(2) = cos(f1_the);

  a02(0) = sin(f2_the) * cos(f2_phi);
  a02(1) = sin(f2_the) * sin(f2_phi);
  a02(2) = cos(f2_the);
}

MaterialModel_GOH06_Incompressible_Mixed::MaterialModel_GOH06_Incompressible_Mixed(
    const char * const &fname)
: pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ), pi( MATH_T::PI )
{
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_GOH06_Incompressible_Mixed constructor does not match h5 file.\n" );

  rho0   = h5r -> read_doubleScalar("/", "rho0");
  E      = h5r -> read_doubleScalar("/", "E");
  nu     = h5r -> read_doubleScalar("/", "nu");
  mu     = h5r -> read_doubleScalar("/", "mu");
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
  
  a01(0) = sin(f1_the) * cos(f1_phi);
  a01(1) = sin(f1_the) * sin(f1_phi);
  a01(2) = cos(f1_the);

  a02(0) = sin(f2_the) * cos(f2_phi);
  a02(1) = sin(f2_the) * sin(f2_phi);
  a02(2) = cos(f2_the);
}

MaterialModel_GOH06_Incompressible_Mixed::~MaterialModel_GOH06_Incompressible_Mixed()
{}

void MaterialModel_GOH06_Incompressible_Mixed::print_info() const
{
  SYS_T::commPrint("\t  MaterialModel_GOH06_Incompressible_Mixed: \n");
  SYS_T::commPrint("\t  Density rho  = %e \n", rho0);
  SYS_T::commPrint("\t  Ground Matrix Neo-Hookean: \n");
  SYS_T::commPrint("\t  Young's Modulus E  = %e \n", E);
  SYS_T::commPrint("\t  Possion's ratio nu = %e \n", nu);
  SYS_T::commPrint("\t  Shear modulus mu   = %e \n", mu);
  SYS_T::commPrint("\t  Fibre Fung: \n");
  SYS_T::commPrint("\t  Angle theta_1 (deg)= %e \n", f1_the*180/pi);
  SYS_T::commPrint("\t  Angle phi_1 (deg)  = %e \n", f1_phi*180/pi);
  SYS_T::commPrint("\t  Angle theta_1 (rad)= %e \n", f1_the);
  SYS_T::commPrint("\t  Angle phi_1 (rad)  = %e \n", f1_phi);
  SYS_T::commPrint("\t  Angle theta_2 (deg)= %e \n", f2_the*180/pi);
  SYS_T::commPrint("\t  Angle phi_2 (deg)  = %e \n", f2_phi*180/pi);
  SYS_T::commPrint("\t  Angle theta_2 (rad)= %e \n", f2_the);
  SYS_T::commPrint("\t  Angle phi_2 (rad)  = %e \n", f2_phi);
  SYS_T::commPrint("\t  a1: [ %e , %e , %e ] \n", a1(0), a1(1), a1(2));
  SYS_T::commPrint("\t  a2: [ %e , %e , %e ] \n", a2(0), a2(1), a2(2));
  SYS_T::commPrint("\t  Fibre k1   = %e \n", fk1);
  SYS_T::commPrint("\t  Fibre k2   = %e \n", fk2);
  SYS_T::commPrint("\t  Fibre k_dispersion = %e \n", fkd);
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
    h5w -> write_doubleScalar("mu", mu);
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

void MaterialModel_GOH06_Incompressible_Mixed::get_PK( 
    const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  const double trC = C.tr();
  const double detFm0d67 = std::pow(F.det(), mpt67);

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  const double fE1 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a1Ca1) - 1.0;
  const double fE2 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a2Ca2) - 1.0;
  
  const double dfpsi1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1 );
  const double dfpsi2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2 );

  // PxI = P : I = I - 1/3 trC C^-1
  Tensor2_3D PxI(Cinv); PxI.scale( (-1.0) * pt33 * trC );
  PxI.AXPI( 1.0 );

  // PxH1(2) = P : H = kd PxI + (1-3kd) a x a + (kd - 1/3) (a.Ca) C^-1
  Tensor2_3D PxH1, PxH2;
  PxH1.gen_zero(); PxH2.gen_zero();
  PxH1.AXPY( fkd, PxI );
  PxH2.AXPY( fkd, PxI );

  Tensor2_3D a1xa1, a2xa2;
  a1xa1.gen_outprod(a1);
  a2xa2.gen_outprod(a2);

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
    const Tensor2_3D &F, Tensor2_3D &P, Tensor2_3D &S, Tensor4_3D &CC ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  const double trC = C.tr();
  const double detFm0d67 = std::pow(F.det(), mpt67);

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  const double fE1 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a1Ca1) - 1.0;
  const double fE2 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a2Ca2) - 1.0;
  
  const double dfpsi1 = fk1 * fE1 * std::exp( fk2 * fE1 * fE1 );
  const double dfpsi2 = fk1 * fE2 * std::exp( fk2 * fE2 * fE2 );

  // P : I = I - 1/3 trC C^-1
  Tensor2_3D PxI(Cinv); PxI.scale( (-1.0) * pt33 * trC );
  PxI.AXPI( 1.0 );

  // P : H = kd PxI + (1-3kd) a x a + (kd - 1/3) (a.Ca) C^-1
  Tensor2_3D PxH1, PxH2;
  PxH1.gen_zero(); PxH2.gen_zero();
  PxH1.AXPY( fkd, PxI );
  PxH2.AXPY( fkd, PxI );

  Tensor2_3D a1xa1, a2xa2;
  a1xa1.gen_outprod(a1);
  a2xa2.gen_outprod(a2);

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
    const Tensor2_3D &F ) const
{
  Tensor2_3D C; C.MatMultTransposeLeft(F);
  Tensor2_3D Cinv = Ten2::inverse(C);
  const double trC = C.tr();
  const double detFm0d67 = std::pow(F.det(), mpt67);

  const double a1Ca1 = C.VecMatVec(a1, a1);
  const double a2Ca2 = C.VecMatVec(a2, a2);

  const double fE1 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a1Ca1) - 1.0;
  const double fE2 = detFm0d67 * (fkd*trC + (1.0-3.0*fkd)*a2Ca2) - 1.0;

  const double PSI_iso = 0.5 * mu * (detFm0d67 * trC - 3.0);
  const double PSI_fi1 = 0.5 * (fk1 / fk2) * ( std::exp(fk2*fE1*fE1) - 1.0 );
  const double PSI_fi2 = 0.5 * (fk1 / fk2) * ( std::exp(fk2*fE2*fE2) - 1.0 );

  return PSI_iso + PSI_fi1 + PSI_fi2;
}

Vector_3 MaterialModel_GOH06_Incompressible_Mixed::get_fibre_dir( const int &dir ) const
{
  if(dir == 0)
    return a1;
  else if(dir == 1)
    return a2;
  else
  {
    SYS_T::print_fatal("Error: MaterialModel_GOH06_Incompressible_Mixed, wrong fibre direction. \n");
    return Vector_3();
  }
}

void MaterialModel_GOH06_Incompressible_Mixed::update_fibre_dir( const Vector_3 &basis_r,
    const Vector_3 &basis_c, const Vector_3 &basis_l )
{
  a1 = a01(0) * basis_r + a01(1) * basis_c + a01(2) * basis_l;
  a2 = a02(0) * basis_r + a02(1) * basis_c + a02(2) * basis_l;
}

// EOF
