#include "MaterialModel_Guccione_Incompressible_Mixed.hpp"

MaterialModel_Guccione_Incompressible_Mixed::MaterialModel_Guccione_Incompressible_Mixed(
    const double &in_rho, const double &in_C,
    const double &in_bf, const double &in_bt, const double &in_bft,
    const double &fx, const double &fy, const double &fz,
    const double &sx, const double &sy, const double &sz )
: pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  rho0( in_rho ), Cq(in_C), b_f(in_bf), b_t(in_bt),
  b_ft(in_bft),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{
  f[0] = fx; f[1] = fy; f[2] = fz;

  // Check to make sure f is a unit vector
  if( !MATH_T::equals( MATH_T::norm2(fx,fy,fz), 1.0, 1.0e-12) )
  {
    SYS_T::commPrint("Guccione model, input f vector is not unit.\n");
    MATH_T::normalize3d(f[0], f[1], f[2]);
  }

  // Check to make sure s is a unit vector
  s[0] = sx; s[1] = sy; s[2] = sz;
  if( !MATH_T::equals( MATH_T::norm2(sx, sy, sz), 1.0, 1.0e-12) )
  {
    SYS_T::commPrint("Guccione model, input s vector is not unit.\n");
    MATH_T::normalize3d(s[0], s[1], s[2]);
  }
 
  // f x s / || f x s || = n 
  MATH_T::cross3d( f[0], f[1], f[2], s[0], s[1], s[2], n[0], n[1], n[2] );
  MATH_T::normalize3d( n[0], n[1], n[2] ); 

  // Define the roatation matrix R
  R(0,0) = f[0]; R(0,1) = f[1]; R(0,2) = f[2];
  R(1,0) = s[0]; R(1,1) = s[1]; R(1,2) = s[2];
  R(2,0) = n[0]; R(2,1) = n[1]; R(2,2) = n[2];

  Rt.copy(R); Rt.transpose();
}

MaterialModel_Guccione_Incompressible_Mixed::MaterialModel_Guccione_Incompressible_Mixed(
        const char * const &fname )
: pt33( 1.0 / 3.0 ), mpt67( -2.0 * pt33 ),
  I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
{
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );
  
  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_Guccione_Incompressible_Mixed constructor does not match h5 file.\n" );

  rho0 = h5r -> read_doubleScalar("/", "rho0");
  Cq   = h5r -> read_doubleScalar("/", "Cq");
  b_f  = h5r -> read_doubleScalar("/", "b_f");
  b_t  = h5r -> read_doubleScalar("/", "b_t");
  b_ft = h5r -> read_doubleScalar("/", "b_ft");

  const std::vector<double> temp_f = h5r -> read_doubleVector( "/", "f" );
  f[0] = temp_f[0]; f[1] = temp_f[1]; f[2] = temp_f[2];

  const std::vector<double> temp_s = h5r -> read_doubleVector( "/", "s" );
  s[0] = temp_s[0]; s[1] = temp_s[1]; s[2] = temp_s[2];

  delete h5r; H5Fclose(h5file);
  
  // Check to make sure f is a unit vector
  if( !MATH_T::equals( MATH_T::norm2(f[0],f[1],f[2]), 1.0, 1.0e-12) )
  {
    SYS_T::commPrint("Guccione model, input f vector is not unit.\n");
    MATH_T::normalize3d(f[0], f[1], f[2]);
  }

  // Check to make sure s is a unit vector
  if( !MATH_T::equals( MATH_T::norm2(s[0], s[1], s[2]), 1.0, 1.0e-12) )
  {
    SYS_T::commPrint("Guccione model, input s vector is not unit.\n");
    MATH_T::normalize3d(s[0], s[1], s[2]);
  }

  // f x s / || f x s || = n
  MATH_T::cross3d( f[0], f[1], f[2], s[0], s[1], s[2], n[0], n[1], n[2] );
  MATH_T::normalize3d( n[0], n[1], n[2] );

  // Define the roatation matrix R
  R(0,0) = f[0]; R(0,1) = f[1]; R(0,2) = f[2];
  R(1,0) = s[0]; R(1,1) = s[1]; R(1,2) = s[2];
  R(2,0) = n[0]; R(2,1) = n[1]; R(2,2) = n[2];

  Rt.copy(R); Rt.transpose();
}


MaterialModel_Guccione_Incompressible_Mixed::~MaterialModel_Guccione_Incompressible_Mixed()
{}


void MaterialModel_Guccione_Incompressible_Mixed::print_info() const
{
  SYS_T::commPrint("\t  MaterialModel_Guccione_Incompressible_Mixed: \n");
  SYS_T::commPrint("\t  Density rho  = %e \n", rho0);
  SYS_T::commPrint("\t  Para C  = %e \n", Cq);
  SYS_T::commPrint("\t  Para b_f  = %e \n", b_f);
  SYS_T::commPrint("\t  Para b_t  = %e \n", b_t);
  SYS_T::commPrint("\t  Para b_ft = %e \n", b_ft);
  SYS_T::commPrint("\t  Fibre dir = [%e %e %e] \n", f[0], f[1], f[2]); 
  SYS_T::commPrint("\t  Sheet normal dir = [%e %e %e] \n", s[0], s[1], s[2]); 
  SYS_T::commPrint("\t  Third n dir = [%e %e %e] \n", n[0], n[1], n[2]); 
}

void MaterialModel_Guccione_Incompressible_Mixed::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar( "rho0", rho0 );
    h5w -> write_doubleScalar("Cq", Cq);
    h5w -> write_doubleScalar("b_f", b_f);
    h5w -> write_doubleScalar("b_t", b_t);
    h5w -> write_doubleScalar("b_ft", b_ft);
    h5w -> write_doubleVector("f",f,3);
    h5w -> write_doubleVector("s",s,3);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}


void MaterialModel_Guccione_Incompressible_Mixed::get_PK( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();
  trC = C.tr();
  trC2 = C.MatContraction( C );
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  // E_bar = 0.5 * (J^-2/3 C - I )
  E_bar.copy(C); E_bar.scale(0.5 * detFm0d67); E_bar.AXPY(-0.5, I);

  // E* = R^T E_bar R
  E_star.MatMult(Rt, E_bar); E_star.MatMult(E_star, R);

  // PxE_bar = E_bar - 1/6 (J^-2/3 C:C  - trC ) C^-1.
  PxE_bar.copy(E_bar);
  PxE_bar.AXPY(-0.5*pt33*(detFm0d67 * trC2 - trC), Cinv);

  // Calculate Q
  const double Q = b_f * E_star(0) * E_star(0) + b_t * ( E_star(4) * E_star(4)
      + E_star(8) * E_star(8) + E_star(5) * E_star(5) + E_star(7) * E_star(7) ) 
    + b_ft * ( E_star(1) * E_star(1) + E_star(3) * E_star(3) + E_star(2) * E_star(2)
       + E_star(6) * E_star(6) ); 

  const double coeff = detFm0d67 * Cq * std::exp(Q);

  // 2nd PK
  S(0) = coeff * b_f  * PxE_bar(0);
  S(1) = coeff * b_ft * PxE_bar(1);
  S(2) = coeff * b_ft * PxE_bar(2);
  S(3) = coeff * b_ft * PxE_bar(3);
  S(4) = coeff * b_t  * PxE_bar(4);
  S(5) = coeff * b_t  * PxE_bar(5);
  S(6) = coeff * b_ft * PxE_bar(6);
  S(7) = coeff * b_t  * PxE_bar(7);
  S(8) = coeff * b_t  * PxE_bar(8);

  // 1st PK
  P.MatMult(F,S);
}


void MaterialModel_Guccione_Incompressible_Mixed::get_PK_Stiffness( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &CC )
{
  C.MatMultTransposeLeft(F);
  Cinv.copy(C); Cinv.inverse();
  trC = C.tr();
  trC2 = C.MatContraction( C );
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  // E_bar = 0.5 * (J^-2/3 C - I )
  E_bar.copy(C); E_bar.scale(0.5 * detFm0d67); E_bar.AXPY(-0.5, I);

  // E* = R^T E_bar R
  E_star.MatMult(Rt, E_bar); E_star.MatMult(E_star, R);

  // PxE_bar = E_bar - 1/6 (J^-2/3 C:C  - trC ) C^-1.
  PxE_bar.copy(E_bar);
  PxE_bar.AXPY(-0.5*pt33*(detFm0d67 * trC2 - trC), Cinv);

  // Calculate Q
  const double Q = b_f * E_star(0) * E_star(0) + b_t * ( E_star(4) * E_star(4)
      + E_star(8) * E_star(8) + E_star(5) * E_star(5) + E_star(7) * E_star(7) ) 
    + b_ft * ( E_star(1) * E_star(1) + E_star(3) * E_star(3) + E_star(2) * E_star(2)
       + E_star(6) * E_star(6) ); 

  const double coeff = detFm0d67 * Cq * std::exp(Q);

  // 2nd PK
  S(0) = coeff * b_f  * PxE_bar(0);
  S(1) = coeff * b_ft * PxE_bar(1);
  S(2) = coeff * b_ft * PxE_bar(2);
  S(3) = coeff * b_ft * PxE_bar(3);
  S(4) = coeff * b_t  * PxE_bar(4);
  S(5) = coeff * b_t  * PxE_bar(5);
  S(6) = coeff * b_ft * PxE_bar(6);
  S(7) = coeff * b_t  * PxE_bar(7);
  S(8) = coeff * b_t  * PxE_bar(8);

  // 1st PK
  P.MatMult(F,S);
  
  // C_tilde = J^-4/3 C_q exp(Q) * b * symmid
  CC.gen_symm_id(); CC.scale(coeff * detFm0d67);
  for(int kk=0; kk<3; ++kk)
  {
    for(int ll=0; ll<3; ++ll)
    {
      CC(0,0,kk,ll) *= b_f;
      CC(0,1,kk,ll) *= b_ft;
      CC(0,2,kk,ll) *= b_ft;
      CC(1,0,kk,ll) *= b_ft;
      CC(1,1,kk,ll) *= b_t;
      CC(1,2,kk,ll) *= b_t;
      CC(2,0,kk,ll) *= b_ft;
      CC(2,1,kk,ll) *= b_t;
      CC(2,2,kk,ll) *= b_t;
    }
  }

  // P^T : C_tilde : P
  Tensor4_3D PP; PP.gen_P( C, Cinv ); CC.TenPMult( PP );

  const double val = coeff * ( b_f * E_bar(0) * C(0)
      + b_ft *(E_bar(1)*C(1) + E_bar(2)*C(2) + E_bar(3)*C(3) + E_bar(6)*C(6))
      + b_t * (E_bar(4)*C(4) + E_bar(5)*C(5) + E_bar(7)*C(7) + E_bar(8)*C(8)) );

  CC.add_SymmProduct(val, Cinv, Cinv);
  CC.add_OutProduct((-1.0)*pt33*val, Cinv, Cinv);

  CC.add_OutProduct(mpt67, Cinv, S);
  CC.add_OutProduct(mpt67, S, Cinv);
}


double MaterialModel_Guccione_Incompressible_Mixed::get_strain_energy( 
    const Matrix_3x3 &F )
{
  C.MatMultTransposeLeft(F);
  detF = F.det();
  detFm0d67 = std::pow(detF, mpt67);

  // E_bar = 0.5 * (J^-2/3 C - I )
  E_bar.copy(C); E_bar.scale(0.5 * detFm0d67); E_bar.AXPY(-0.5, I);

  // E* = R^T E_bar R
  E_star.MatMult(Rt, E_bar); E_star.MatMult(E_star, R);
  
  const double Q = b_f * E_star(0) * E_star(0) + b_t * ( E_star(4) * E_star(4)
      + E_star(8) * E_star(8) + E_star(5) * E_star(5) + E_star(7) * E_star(7) ) 
    + b_ft * ( E_star(1) * E_star(1) + E_star(3) * E_star(3) + E_star(2) * E_star(2)
       + E_star(6) * E_star(6) ); 

  return 0.5 * Cq * ( std::exp(Q) - 1.0 );
}

// EOF
