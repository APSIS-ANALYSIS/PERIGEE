#include "PLocAssem_Elastodynamics_GenAlpha.hpp"

PLocAssem_Elastodynamics_GenAlpha::PLocAssem_Elastodynamics_GenAlpha(
    const double &in_rho, const double &in_module_E, const double &in_nu,
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_snlocbas, 
    const int &in_num_ebc_fun )
: rho( in_rho ), module_E( in_module_E ), nu( in_nu ),
  lambda( in_nu * in_module_E / ((1.0 + in_nu) * (1.0 - 2.0 * in_nu)) ),
  mu( 0.5 * in_module_E / (1.0 + in_nu) ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), num_ebc_fun( in_num_ebc_fun ),
  nLocBas( in_nlocbas ), snLocBas( in_snlocbas ),
  vec_size( 3 * in_nlocbas ), sur_size( 3 * in_snlocbas )
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];
  
  sur_Residual = new PetscScalar[sur_size];

  Zero_Tangent_Residual();

  Zero_sur_Residual();

  if( num_ebc_fun == 0 ) flist = nullptr;
  else flist = new locassem_transport_funs [num_ebc_fun];

  //flist[0] = &PLocAssem_Elastodynamics_GenAlpha::get_g_0;
  //flist[1] = &PLocAssem_Elastodynamics_GenAlpha::get_g_1;
 
  print_info();
}

PLocAssem_Elastodynamics_GenAlpha::~PLocAssem_Elastodynamics_GenAlpha()
{
  delete [] Tangent; Tangent = nullptr;
  delete [] Residual; Residual = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;

  if(num_ebc_fun > 0) delete [] flist;
}

void PLocAssem_Elastodynamics_GenAlpha::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional elastodynamics equation: \n");
  if(nLocBas == 4)
    SYS_T::commPrint("  FEM: 4-node Tetrahedral element \n");
  else if(nLocBas == 10)
    SYS_T::commPrint("  FEM: 10-node Tetrahedral element \n");
  else if(nLocBas == 8)
    SYS_T::commPrint("  FEM: 8-node Hexahedral element \n");
  else if(nLocBas == 27)
    SYS_T::commPrint("  FEM: 27-node Hexahedral element \n");
  else SYS_T::print_fatal("Error: unknown elem type.\n");
  SYS_T::commPrint("  Spatial: finite element \n");
  SYS_T::commPrint("  Temporal: Generalized-alpha Method \n");
  SYS_T::commPrint("  Consistent tangent matrix used. \n");
  SYS_T::print_sep_line();
}

void PLocAssem_Elastodynamics_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol_velo,
    const double * const &sol_disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );
  
  const double curr = time + alpha_f * dt;

  Zero_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double ux_x = 0.0, ux_y = 0.0, ux_z = 0.0;
    double uy_x = 0.0, uy_y = 0.0, uy_z = 0.0;
    double uz_x = 0.0, uz_y = 0.0, uz_z = 0.0;

    double vx_t = 0.0, vy_t = 0.0, vz_t = 0.0;

    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii3 = 3 * ii;

      ux_x += sol_disp[ii3  ] * dR_dx[ii];
      uy_x += sol_disp[ii3+1] * dR_dx[ii];
      uz_x += sol_disp[ii3+2] * dR_dx[ii];

      ux_y += sol_disp[ii3  ] * dR_dy[ii];
      uy_y += sol_disp[ii3+1] * dR_dy[ii];
      uz_y += sol_disp[ii3+2] * dR_dy[ii];

      ux_z += sol_disp[ii3  ] * dR_dz[ii];
      uy_z += sol_disp[ii3+1] * dR_dz[ii];
      uz_z += sol_disp[ii3+2] * dR_dz[ii];

      vx_t += dot_sol_velo[ii3  ] * R[ii];
      vy_t += dot_sol_velo[ii3+1] * R[ii];
      vz_t += dot_sol_velo[ii3+2] * R[ii];
      
      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }
    
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    const Vector_3 f_body = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[3*A  ] += gwts * ( NA * rho * vx_t
          + NA_x * ((2.0 * mu + lambda) * ux_x + lambda * (uy_y + uz_z))
          + NA_y * lambda * (ux_y + uy_x)
          + NA_z * lambda * (ux_z + uz_x)
          - NA * rho * f_body.x() );

      Residual[3*A+1] += gwts * ( NA * rho * vy_t
          + NA_x * lambda * (ux_y + uy_x)
          + NA_y * ((2.0 * mu + lambda) * uy_y + lambda * (ux_x + uz_z))
          + NA_z * lambda * (uy_z + uz_y)
          - NA * rho * f_body.y() );

      Residual[3*A+2] += gwts * ( NA * rho * vz_t
          + NA_x * lambda * (ux_z + uz_x)
          + NA_y * lambda * (uy_z + uz_y)
          + NA_z * ((2.0 * mu + lambda) * uz_z + lambda * (ux_x + uy_y))
          - NA * rho * f_body.z() );
    }
  }
}

void PLocAssem_Elastodynamics_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol_velo,
    const double * const &sol_disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();
  
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );
  
  const double curr = time + alpha_f * dt;

  const double dd_dv = alpha_f * alpha_f * gamma * gamma * dt * dt / alpha_m;

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double ux_x = 0.0, ux_y = 0.0, ux_z = 0.0;
    double uy_x = 0.0, uy_y = 0.0, uy_z = 0.0;
    double uz_x = 0.0, uz_y = 0.0, uz_z = 0.0;

    double vx_t = 0.0, vy_t = 0.0, vz_t = 0.0;
    
    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int ii3 = 3 * ii;

      ux_x += sol_disp[ii3  ] * dR_dx[ii];
      uy_x += sol_disp[ii3+1] * dR_dx[ii];
      uz_x += sol_disp[ii3+2] * dR_dx[ii];

      ux_y += sol_disp[ii3  ] * dR_dy[ii];
      uy_y += sol_disp[ii3+1] * dR_dy[ii];
      uz_y += sol_disp[ii3+2] * dR_dy[ii];

      ux_z += sol_disp[ii3  ] * dR_dz[ii];
      uy_z += sol_disp[ii3+1] * dR_dz[ii];
      uz_z += sol_disp[ii3+2] * dR_dz[ii];

      vx_t += dot_sol_velo[ii3  ] * R[ii];
      vy_t += dot_sol_velo[ii3+1] * R[ii];
      vz_t += dot_sol_velo[ii3+2] * R[ii];
      
      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }
    
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    const Vector_3 f_body = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[3*A  ] += gwts * ( NA * rho * vx_t
          + NA_x * ((2.0 * mu + lambda) * ux_x + lambda * (uy_y + uz_z))
          + NA_y * lambda * (ux_y + uy_x)
          + NA_z * lambda * (ux_z + uz_x)
          - NA * rho * f_body.x() );

      Residual[3*A+1] += gwts * ( NA * rho * vy_t
          + NA_x * lambda * (ux_y + uy_x)
          + NA_y * ((2.0 * mu + lambda) * uy_y + lambda * (ux_x + uz_z))
          + NA_z * lambda * (uy_z + uz_y)
          - NA * rho * f_body.y() );

      Residual[3*A+2] += gwts * ( NA * rho * vz_t
          + NA_x * lambda * (ux_z + uz_x)
          + NA_y * lambda * (uy_z + uz_y)
          + NA_z * ((2.0 * mu + lambda) * uz_z + lambda * (ux_x + uy_y))
          - NA * rho * f_body.z() );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];

        Tangent[9*nLocBas*A+3*B  ] += gwts * (NA * rho * alpha_m * NB 
            + dd_dv * ((2.0 * mu + lambda) * NA_x * NB_x
            + lambda * (NA_y * NB_y + NA_z * NB_z)));

        Tangent[9*nLocBas*A+3*B+1] += gwts * dd_dv * lambda * (NA_x * NB_y + NA_y * NB_x);

        Tangent[9*nLocBas*A+3*B+2] += gwts * dd_dv * lambda * (NA_x * NB_z + NA_z * NB_x);

        Tangent[3*nLocBas*(3*A+1)+3*B  ] += gwts * dd_dv * lambda * (NA_x * NB_y + NA_y * NB_x);

        Tangent[3*nLocBas*(3*A+1)+3*B+1] += gwts * ( NA * rho * alpha_m * NB 
            + dd_dv * ((2.0 * mu + lambda) * NA_y * NB_y
            + lambda * (NA_x * NB_x + NA_z * NB_z)));

        Tangent[3*nLocBas*(3*A+1)+3*B+2] += gwts * dd_dv * lambda * (NA_y * NB_z + NA_z * NB_y);

        Tangent[3*nLocBas*(3*A+2)+3*B  ] += gwts * dd_dv * lambda * (NA_x * NB_z + NA_z * NB_x);

        Tangent[3*nLocBas*(3*A+2)+3*B+1] += gwts * dd_dv * lambda * (NA_y * NB_z + NA_z * NB_y);

        Tangent[3*nLocBas*(3*A+2)+3*B+2] += gwts * (NA * rho * alpha_m * NB 
            + dd_dv * ((2.0 * mu + lambda) * NA_z * NB_z
            + lambda * (NA_x * NB_x + NA_y * NB_y)));     
      } 
    }
  } // End-of-quadrature-loop
}

void PLocAssem_Elastodynamics_GenAlpha::Assem_Mass_Residual(
    const double * const &sol_disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();
  
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const double curr = 0.0;

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double ux_x = 0.0, ux_y = 0.0, ux_z = 0.0;
    double uy_x = 0.0, uy_y = 0.0, uy_z = 0.0;
    double uz_x = 0.0, uz_y = 0.0, uz_z = 0.0;

    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      int ii3 = ii * 3;

      ux_x += sol_disp[ii3  ] * dR_dx[ii];
      uy_x += sol_disp[ii3+1] * dR_dx[ii];
      uz_x += sol_disp[ii3+2] * dR_dx[ii];

      ux_y += sol_disp[ii3  ] * dR_dy[ii];
      uy_y += sol_disp[ii3+1] * dR_dy[ii];
      uz_y += sol_disp[ii3+2] * dR_dy[ii];

      ux_z += sol_disp[ii3  ] * dR_dz[ii];
      uy_z += sol_disp[ii3+1] * dR_dz[ii];
      uz_z += sol_disp[ii3+2] * dR_dz[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const Vector_3 f_body = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[3*A  ] += gwts * (
          + NA_x * ((2.0 * mu + lambda) * ux_x + lambda * (uy_y + uz_z))
          + NA_y * lambda * (ux_y + uy_x)
          + NA_z * lambda * (ux_z + uz_x)
          - NA * rho * f_body.x() );

      Residual[3*A+1] += gwts * (
          + NA_x * lambda * (ux_y + uy_x)
          + NA_y * ((2.0 * mu + lambda) * uy_y + lambda * (ux_x + uz_z))
          + NA_z * lambda * (uy_z + uz_y)
          - NA * rho * f_body.y() );

      Residual[3*A+2] += gwts * (
          + NA_x * lambda * (ux_z + uz_x)
          + NA_y * lambda * (uy_z + uz_y)
          + NA_z * ((2.0 * mu + lambda) * uz_z + lambda * (ux_x + uy_y))
          - NA * rho * f_body.z() );

      for(int B=0; B<nLocBas; ++B)
      {
        Tangent[3*nLocBas*(3*A) + 3*B    ] += gwts * rho * NA * R[B];
        Tangent[3*nLocBas*(3*A+1) + 3*B+1] += gwts * rho * NA * R[B];
        Tangent[3*nLocBas*(3*A+2) + 3*B+2] += gwts * rho * NA * R[B];
      }
    }
  } // End-of-quadrature-loop
}

void PLocAssem_Elastodynamics_GenAlpha::Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  const int face_nqp = quad -> get_num_quadPts();

  Zero_sur_Residual();

  const double curr = time + alpha_f * dt;

  for(int qua=0; qua < face_nqp; ++qua)
  {
    const std::vector<double> R = element->get_R(qua);

    double surface_area;
    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    Vector_3 coor(0.0, 0.0, 0.0);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const Vector_3 traction = get_ebc_fun( ebc_id, coor, curr, n_out );

    const double gwts = surface_area * quad -> get_qw( qua );

    for(int A=0; A<snLocBas; ++A)
    {
      sur_Residual[3*A  ] -= gwts * R[A] * traction.x();
      sur_Residual[3*A+1] -= gwts * R[A] * traction.y();
      sur_Residual[3*A+2] -= gwts * R[A] * traction.z();
    }
  }
}

// EOF
