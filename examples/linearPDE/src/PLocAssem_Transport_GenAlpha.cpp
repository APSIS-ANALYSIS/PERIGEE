#include "PLocAssem_Transport_GenAlpha.hpp"

PLocAssem_Transport_GenAlpha::PLocAssem_Transport_GenAlpha(
    const double &in_rho, const double &in_cap, const double &in_kappa,
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_snlocbas, 
    const int &in_num_ebc_fun )
: rho( in_rho ), cap( in_cap ), kappa( in_kappa ), 
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), num_ebc_fun( in_num_ebc_fun ),
  nLocBas( in_nlocbas ), snLocBas( in_snlocbas ),
  vec_size( in_nlocbas ), sur_size( in_snlocbas )
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];
  
  sur_Residual = new PetscScalar[sur_size];

  Zero_Tangent_Residual();

  Zero_sur_Residual();

  if( num_ebc_fun == 0 ) flist = nullptr;
  else flist = new locassem_transport_funs [num_ebc_fun];

  //flist[0] = &PLocAssem_Transport_GenAlpha::get_g_0;
  //flist[1] = &PLocAssem_Transport_GenAlpha::get_g_1;
 
  print_info();
}

PLocAssem_Transport_GenAlpha::~PLocAssem_Transport_GenAlpha()
{
  delete [] Tangent; Tangent = nullptr;
  delete [] Residual; Residual = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;

  if(num_ebc_fun > 0) delete [] flist;
}

void PLocAssem_Transport_GenAlpha::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional transport equation: \n");
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

void PLocAssem_Transport_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
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
    double u_t = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;

    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for(int ii=0; ii<nLocBas; ++ii)
    {
      u_t += dot_sol[ii] * R[ii];
      u_x += sol[ii]     * dR_dx[ii];
      u_y += sol[ii]     * dR_dy[ii];
      u_z += sol[ii]     * dR_dz[ii];
      
      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }
    
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    const double ff = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      Residual[A] += gwts * ( R[A] * (rho * cap * u_t - ff) + kappa * dR_dx[A] * u_x
         + kappa * dR_dy[A] * u_y + kappa * dR_dz[A] * u_z ); 
    }
  }
}

void PLocAssem_Transport_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();
  
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );
  
  const double curr = time + alpha_f * dt;

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double u_t = 0.0, u_x = 0.0, u_y = 0.0, u_z = 0.0;
    
    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for(int ii=0; ii<nLocBas; ++ii)
    {
      u_t += dot_sol[ii] * R[ii];
      u_x += sol[ii]     * dR_dx[ii];
      u_y += sol[ii]     * dR_dy[ii];
      u_z += sol[ii]     * dR_dz[ii];
      
      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }
    
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    const double ff = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[A] += gwts * ( NA * (rho * cap * u_t - ff) + kappa * ( NA_x * u_x
         + NA_y * u_y + NA_z * u_z ) );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B], NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];

        Tangent[nLocBas * A + B] += gwts * ( alpha_m * rho * cap * NA * NB
           + alpha_f * gamma * dt * kappa * (NA_x * NB_x + NA_y * NB_y + NA_z * NB_z) ); 
      } 
    }
  } // End-of-quadrature-loop
}

void PLocAssem_Transport_GenAlpha::Assem_Mass_Residual(
    const double * const &sol,
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
    double u_x = 0.0, u_y = 0.0, u_z = 0.0;
    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      u_x += sol[ii]     * dR_dx[ii];
      u_y += sol[ii]     * dR_dy[ii];
      u_z += sol[ii]     * dR_dz[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const double ff = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A], NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      Residual[A] += gwts * ( kappa * ( NA_x * u_x + NA_y * u_y + NA_z * u_z ) - NA * ff );

      for(int B=0; B<nLocBas; ++B)
        Tangent[nLocBas * A + B] += gwts * rho * cap * R[A] * R[B];
    }
  } // End-of-quadrature-loop
}

void PLocAssem_Transport_GenAlpha::Assem_Residual_EBC(
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

    const double gg = get_ebc_fun( ebc_id, coor, curr );

    const double gwts = surface_area * quad -> get_qw( qua );

    for(int A=0; A<snLocBas; ++A)
      sur_Residual[A] -= gwts * R[A] * gg;
  }
}

// EOF