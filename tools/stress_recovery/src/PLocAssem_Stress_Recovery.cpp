#include "PLocAssem_Stress_Recovery.hpp"

PLocAssem_Stress_Recovery::PLocAssem_Stress_Recovery( const int &in_nlocbas )
: nLocBas( in_nlocbas ), vec_size( 9 * nLocBas )
{
  Tangent = new PetscScalar[ vec_size * vec_size ];
  Residual = new PetscScalar[ vec_size ];

  Zero_Tangent_Residual();

  print_info();
}

PLocAssem_Stress_Recovery::~PLocAssem_Stress_Recovery()
{
  delete [] Tangent; Tangent = nullptr;
  delete [] Residual; Residual = nullptr;
}

void PLocAssem_Stress_Recovery::print_info() const
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
  SYS_T::print_sep_line();
}

void PLocAssem_Stress_Recovery::Assem_Residual(
    const double * const &isol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double ux_x = 0.0, ux_y = 0.0, ux_z = 0.0;
    double uy_x = 0.0, uy_y = 0.0, uy_z = 0.0;
    double uz_x = 0.0, uz_y = 0.0, uz_z = 0.0;

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int offset = 3 * ii;

      ux_x += isol[offset  ] * dR_dx[ii];
      uy_x += isol[offset+1] * dR_dx[ii];
      uz_x += isol[offset+2] * dR_dx[ii];

      ux_y += isol[offset  ] * dR_dy[ii];
      uy_y += isol[offset+1] * dR_dy[ii];
      uz_y += isol[offset+2] * dR_dy[ii];

      ux_z += isol[offset  ] * dR_dz[ii];
      uy_z += isol[offset+1] * dR_dz[ii];
      uz_z += isol[offset+2] * dR_dz[ii];
    }

    const Tensor2_3D Grad_disp(ux_x, ux_y, ux_z,
                               uy_x, uy_y, uy_z,
                               uz_x, uz_y, uz_z);

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A];

      Residual[9*A  ] += gwts * NA * Grad_disp(0);
      Residual[9*A+1] += gwts * NA * Grad_disp(1);
      Residual[9*A+2] += gwts * NA * Grad_disp(2);
      Residual[9*A+3] += gwts * NA * Grad_disp(3);
      Residual[9*A+4] += gwts * NA * Grad_disp(4);
      Residual[9*A+5] += gwts * NA * Grad_disp(5);
      Residual[9*A+6] += gwts * NA * Grad_disp(6);
      Residual[9*A+7] += gwts * NA * Grad_disp(7);
      Residual[9*A+8] += gwts * NA * Grad_disp(8);
    }
  }
}

void PLocAssem_Stress_Recovery::Assem_Mass_Residual(
    const double * const &isol,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const IQuadPts * const &quad )
{
  const int nqp = quad -> get_num_quadPts();

  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );

  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqp; ++qua)
  {
    double ux_x = 0.0, ux_y = 0.0, ux_z = 0.0;
    double uy_x = 0.0, uy_y = 0.0, uy_z = 0.0;
    double uz_x = 0.0, uz_y = 0.0, uz_z = 0.0;

    element->get_R_gradR( qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0] );
    
    for(int ii=0; ii<nLocBas; ++ii)
    {
      const int offset = 3 * ii;

      ux_x += isol[offset  ] * dR_dx[ii];
      uy_x += isol[offset+1] * dR_dx[ii];
      uz_x += isol[offset+2] * dR_dx[ii];

      ux_y += isol[offset  ] * dR_dy[ii];
      uy_y += isol[offset+1] * dR_dy[ii];
      uz_y += isol[offset+2] * dR_dy[ii];

      ux_z += isol[offset  ] * dR_dz[ii];
      uy_z += isol[offset+1] * dR_dz[ii];
      uz_z += isol[offset+2] * dR_dz[ii];
    }

    const Tensor2_3D Grad_disp(ux_x, ux_y, ux_z,
                               uy_x, uy_y, uy_z,
                               uz_x, uz_y, uz_z);

    const double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    
    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A];

      Residual[9*A  ] += gwts * NA * Grad_disp(0);
      Residual[9*A+1] += gwts * NA * Grad_disp(1);
      Residual[9*A+2] += gwts * NA * Grad_disp(2);
      Residual[9*A+3] += gwts * NA * Grad_disp(3);
      Residual[9*A+4] += gwts * NA * Grad_disp(4);
      Residual[9*A+5] += gwts * NA * Grad_disp(5);
      Residual[9*A+6] += gwts * NA * Grad_disp(6);
      Residual[9*A+7] += gwts * NA * Grad_disp(7);
      Residual[9*A+8] += gwts * NA * Grad_disp(8);

      for(int B=0; B<nLocBas; ++B)
      {
        // Non-zero diagonal of each A-B block
        for(int C=0; C<9; ++C)
        {
          Tangent[9*nLocBas*(9*A + C) + 9*B + C] += gwts * NA * R[B];
        }
      }
    }
  }
}

//EOF