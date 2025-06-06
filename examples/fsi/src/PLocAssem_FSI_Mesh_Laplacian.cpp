#include "PLocAssem_FSI_Mesh_Laplacian.hpp"

PLocAssem_FSI_Mesh_Laplacian::PLocAssem_FSI_Mesh_Laplacian( 
  const FEType &in_type, const int &in_nqp_v, const int &in_nqp_s )
: elemType(in_type), nqpv(in_nqp_v), nqps(in_nqp_s),
  elementv( ElementFactory::createVolElement(elemType, nqpv) ),
  elements( ElementFactory::createSurElement(elemType, nqps) ),
  quadv( QuadPtsFactory::createVolQuadrature(elemType, nqpv) ),
  quads( QuadPtsFactory::createSurQuadrature(elemType, nqps) ),
  nLocBas( elementv->get_nLocBas() ), snLocBas( elements->get_nLocBas() ), 
  num_ebc_fun(0), vec_size(nLocBas*3)
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();

  if( num_ebc_fun == 0 ) flist = nullptr;
  else flist = new locassem_fsi_mesh_lap_funs [num_ebc_fun];

  print_info();
}

PLocAssem_FSI_Mesh_Laplacian::~PLocAssem_FSI_Mesh_Laplacian()
{
  delete [] Tangent; delete [] Residual; Tangent = nullptr; Residual = nullptr;
  if(num_ebc_fun > 0) delete [] flist;
}

void PLocAssem_FSI_Mesh_Laplacian::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional Laplacian equation: \n");
  elementv->print_info();
  SYS_T::commPrint("  Spatial: Galerkin Finite element \n");
  SYS_T::commPrint("  This solver is for the fluid sub-domain mesh motion in FSI problems.\n");
  SYS_T::print_sep_line();
}

void PLocAssem_FSI_Mesh_Laplacian::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 0.0;
}

void PLocAssem_FSI_Mesh_Laplacian::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}

void PLocAssem_FSI_Mesh_Laplacian::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii) Tangent[ii] = 1.0;
}

void PLocAssem_FSI_Mesh_Laplacian::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z )
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Tangent_Residual();
  
  const double curr = time;
  
  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqpv; ++qua)
  {
    elementv->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    double ux = 0.0, uy = 0.0, uz = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double wx = 0.0, wy = 0.0, wz = 0.0;
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
    {
      ux += vec_b[ii*3+0] * dR_dx[ii];
      uy += vec_b[ii*3+0] * dR_dy[ii];
      uz += vec_b[ii*3+0] * dR_dz[ii];

      vx += vec_b[ii*3+1] * dR_dx[ii];
      vy += vec_b[ii*3+1] * dR_dy[ii];
      vz += vec_b[ii*3+1] * dR_dz[ii];

      wx += vec_b[ii*3+2] * dR_dx[ii];
      wy += vec_b[ii*3+2] * dR_dy[ii];
      wz += vec_b[ii*3+2] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    
    const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);
    
    const Vector_3 f_body = get_f(coor_x, coor_y, coor_z, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[3*A]   += gwts * ( NA_x * ux + NA_y * uy + NA_z * uz - R[A] * f_body.x() );
      
      Residual[3*A+1] += gwts * ( NA_x * vx + NA_y * vy + NA_z * vz - R[A] * f_body.y() );

      Residual[3*A+2] += gwts * ( NA_x * wx + NA_y * wy + NA_z * wz - R[A] * f_body.z() );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];

        Tangent[ 3*nLocBas*(3*A+0) + 3*B + 0 ] += gwts * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);

        Tangent[ 3*nLocBas*(3*A+1) + 3*B + 1 ] += gwts * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);

        Tangent[ 3*nLocBas*(3*A+2) + 3*B + 2 ] += gwts * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z);
      }
    }
  }
}

void PLocAssem_FSI_Mesh_Laplacian::Assem_Residual(
    const double &time, const double &dt,
    const double * const &vec_a,
    const double * const &vec_b,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z )
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Residual();

  const double curr = time;
  
  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqpv; ++qua)
  {
    elementv->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    double ux = 0.0, uy = 0.0, uz = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double wx = 0.0, wy = 0.0, wz = 0.0;

    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    for(int ii=0; ii<nLocBas; ++ii)
    {
      ux += vec_b[ii*3+0] * dR_dx[ii];
      uy += vec_b[ii*3+0] * dR_dy[ii];
      uz += vec_b[ii*3+0] * dR_dz[ii];

      vx += vec_b[ii*3+1] * dR_dx[ii];
      vy += vec_b[ii*3+1] * dR_dy[ii];
      vz += vec_b[ii*3+1] * dR_dz[ii];

      wx += vec_b[ii*3+2] * dR_dx[ii];
      wy += vec_b[ii*3+2] * dR_dy[ii];
      wz += vec_b[ii*3+2] * dR_dz[ii];

      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
      coor_z += eleCtrlPts_z[ii] * R[ii];
    }
    
    const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);

    const Vector_3 f_body = get_f(coor_x, coor_y, coor_z, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[3*A  ] += gwts * ( NA_x * ux + NA_y * uy + NA_z * uz - R[A] * f_body.x() );
      
      Residual[3*A+1] += gwts * ( NA_x * vx + NA_y * vy + NA_z * vz - R[A] * f_body.y() );

      Residual[3*A+2] += gwts * ( NA_x * wx + NA_y * wy + NA_z * wz - R[A] * f_body.z() );
    }
  }
}

// EOF
