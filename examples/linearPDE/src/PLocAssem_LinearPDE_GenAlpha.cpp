#include "PLocAssem_LinearPDE_GenAlpha.hpp"

PLocAssem_LinearPDE_GenAlpha::PLocAssem_LinearPDE_GenAlpha(
    const double &in_Young_modulus, const double &in_Poisson_ratio,
    const double &in_rho,
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_nlocbas, const int &in_snlocbas, 
    const int &in_num_ebc_fun, const int &in_dof,
    const int &in_dof_mat, const int &elemtype )
: Young_modulus( in_Young_modulus ), Possion_ratio( in_Poisson_ratio ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), num_ebc_fun( in_num_ebc_fun ),
  dof( in_dof ), dof_mat( in_dof_mat )
{
  lambda = Possion_ratio * Young_modulus / (1.0 + Possion_ratio) * (1.0 - 2.0 * Possion_ratio);
  mu = 0.5 * Young_modulus / (1.0 + Possion_ratio);

  if(elemtype == 501)
  {
    // 501 is linear tet element
    nLocBas = 4; snLocBas = 3;
  }
  else if(elemtype == 502)
  {
    // 502 is quadratic tet element
    nLocBas = 10; snLocBas = 6;
  }
  else if(elemtype == 601)
  {
    // 601 is tri-linear hex element
    nLocBas = 8; snLocBas = 4;
  }
  else if(elemtype == 602)
  {
    // 602 is tri-quadratic hex element
    nLocBas = 27; snLocBas = 9;
  }
  else SYS_T::print_fatal("Error: unknown elem type.\n");

  vec_size = nLocBas * dof_mat;
  sur_size = snLocBas * dof_mat;

  Mass = new PetscScalar[vec_size * vec_size];
  Stiffness = new PetscScalar[vec_size * vec_size];
  Load = new PetscScalar[vec_size];

  sur_Mass = new PetscScalar[sur_size * sur_size];
  sur_Load = new PetscScalar[sur_size];

  Zero_Mass_Stiffness_Load();

  Zero_sur_Mass_Load();

  if( num_ebc_fun == 0 ) flist = nullptr;
  else flist = new locassem_load_funs [num_ebc_fun];

  if( num_ebc_fun == 0 ) colist = nullptr;
  else colist = new locassem_robin_coefficient [num_ebc_fun];

  //flist[0] = &PLocAssem_LinearPDE_GenAlpha::get_g_0;
  //flist[1] = &PLocAssem_LinearPDE_GenAlpha::get_g_1;
 
  print_info();
}

PLocAssem_LinearPDE_GenAlpha::~PLocAssem_LinearPDE_GenAlpha()
{
  delete [] Mass; Mass = nullptr;
  delete [] Stiffness; Stiffness = nullptr;
  delete [] Load; Load = nullptr;
  delete [] sur_Mass; sur_Mass = nullptr;
  delete [] sur_Load; sur_Load = nullptr;

  if(num_ebc_fun > 0) delete [] flist;
}

void PLocAssem_LinearPDE_GenAlpha::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Three-dimensional linear problem: \n");
  if(nLocBas == 4)
    SYS_T::commPrint("  FEM: 4-node Tetrahedral element \n");
  else if(nLocBas == 10)
    SYS_T::commPrint("  FEM: 10-node Tetrahedral element \n");
  else if(nLocBas == 8)
    SYS_T::commPrint("  FEM: 8-node Hexagonal element \n");
  else if(nLocBas == 27)
    SYS_T::commPrint("  FEM: 27-node Hexagonal element \n");
  else SYS_T::print_fatal("Error: unknown elem type.\n");
  SYS_T::commPrint("  Spatial: finite element \n");
  SYS_T::commPrint("  Temporal: Generalized-alpha Method \n");
  SYS_T::print_sep_line();
}

void PLocAssem_LinearPDE_GenAlpha::Assem_Load(
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

  Zero_Load();

  std::vector<double> R(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    Vector_3 coor(0.0, 0.0, 0.0);

    element->get_R( qua, &R[0]);
    
    for(int ii=0; ii<nLocBas; ++ii)
    { 
      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }
    
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    const Vector_3 f_body = get_f(coor, curr);

    for(int A=0; A<nLocBas; ++A)
    {
      Load[dof_mat*A  ] += gwts * R[A] * f_body.x();
      Load[dof_mat*A+1] += gwts * R[A] * f_body.y();
      Load[dof_mat*A+2] += gwts * R[A] * f_body.z();
    }
  }
}


void PLocAssem_LinearPDE_GenAlpha::Assem_Stiffness(
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

  Zero_Stiffness();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    element->get_gradR( qua, &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    const double row = dof_mat * nLocBas;
    const double row_mat = dof_mat * row;

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB_x = dR_dx[B], NB_y = dR_dy[B], NB_z = dR_dz[B];

        // 11
        Stiffness[row_mat*A+dof_mat*B        ] += gwts * ( (2.0 * mu + lambda) * NA_x * NB_x
            + mu * (NA_y * NB_y + NA_z * NB_z) );
        
        // 12
        Stiffness[row_mat*A+dof_mat*B+1      ] += gwts * ( lambda * NA_x * NB_y + mu * NA_y * NB_x );

        // 13
        Stiffness[row_mat*A+dof_mat*B+2      ] += gwts * ( lambda * NA_x * NB_z + mu * NA_z * NB_x );

        // 21
        Stiffness[row_mat*A+row+dof_mat*B    ] += gwts * ( lambda * NA_y * NB_x + mu * NA_x * NB_y );

        // 22
        Stiffness[row_mat*A+row+dof_mat*B+1  ] += gwts * ( (2.0 * mu + lambda) * NA_y * NB_y
            + mu * (NA_x * NB_x + NA_z * NB_z) );
        
        // 23
        Stiffness[row_mat*A+row+dof_mat*B+2  ] += gwts * ( lambda * NA_y * NB_z + mu * NA_z * NB_y );

        // 31
        Stiffness[row_mat*A+2*row+dof_mat*B  ] += gwts * ( lambda * NA_z * NB_x + mu * NA_x * NB_y );

        // 32
        Stiffness[row_mat*A+2*row+dof_mat*B+1] += gwts * ( lambda * NA_z * NB_y + mu * NA_y * NB_z );

        // 33
        Stiffness[row_mat*A+2*row+dof_mat*B+2] += gwts * ( (2.0 * mu + lambda) * NA_z * NB_z
            + mu * (NA_x * NB_x + NA_y * NB_y) );
      } 
    }
  } // End-of-quadrature-loop
}


void PLocAssem_LinearPDE_GenAlpha::Assem_Mass(
    const double * const &sol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();
  
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Mass();

  std::vector<double> R(nLocBas, 0.0);

  const double row = dof_mat * nLocBas;
  const double row_mat = dof_mat * row;

  for(int qua=0; qua<nqp; ++qua)
  {
    element->get_R( qua, &R[0]);

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A];
      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B];
        // 11
        Mass[row_mat*A+dof_mat*B        ] += gwts * rho * NA * NB;

        // 22
        Mass[row_mat*A+row+dof_mat*B+1  ] += gwts * rho * NA * NB;

        // 33
        Mass[row_mat*A+2*row+dof_mat*B+2] += gwts * rho * NA * NB;
      }
    }
  } // End-of-quadrature-loop
}

void PLocAssem_LinearPDE_GenAlpha::Assem_Load_EBC(
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

  Zero_sur_Load();

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

    const Vector_3 gg = get_ebc_fun( ebc_id, coor, curr );

    const double gwts = surface_area * quad -> get_qw( qua );

    for(int A=0; A<snLocBas; ++A)
    {
      const double NA = R[A];
      sur_Load[dof_mat*A  ] += gwts * NA * gg.x();
      sur_Load[dof_mat*A+1] += gwts * NA * gg.y();
      sur_Load[dof_mat*A+2] += gwts * NA * gg.z();
    }
  }
}

// EOF
