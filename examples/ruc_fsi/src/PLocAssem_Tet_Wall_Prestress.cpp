#include "PLocAssem_Tet_Wall_Prestress.hpp"

PLocAssem_Tet_Wall_Prestress::PLocAssem_Tet_Wall_Prestress(
    const TimeMethod_GenAlpha * const &tm_gAlpha,
    const int &in_face_nqp, const double &in_rho, 
    const double &in_vis_mu, const double &in_wall_rho,
    const double &in_nu, const double &in_kappa,
    const int &elemtype )
: rho0( in_rho ), vis_mu( in_vis_mu ),
  alpha_f(tm_gAlpha->get_alpha_f()), alpha_m(tm_gAlpha->get_alpha_m()),
  gamma(tm_gAlpha->get_gamma()), rho_w(in_wall_rho),
  nu_w(in_nu), kappa_w(in_kappa), face_nqp(in_face_nqp)
{
  if(elemtype == 501)
  {
    // 501 is linear element
    snLocBas = 3;
  }
  else if(elemtype == 502)
  {
    // 502 is quadratic element
    snLocBas = 6;
  }
  else SYS_T::print_fatal("Error: unknown elem type.\n");

  sur_size = snLocBas * 4; // dof_per_node = 4

  // We do not need volumetric element matrix or vector
  Tangent = nullptr;
  Residual = nullptr;

  sur_Tangent = new PetscScalar[sur_size * sur_size];
  sur_Residual = new PetscScalar[sur_size];

  Zero_sur_Tangent_Residual();

  print_info();
}

PLocAssem_Tet_Wall_Prestress::~PLocAssem_Tet_Wall_Prestress()
{
  delete [] sur_Tangent; sur_Tangent = nullptr;
  delete [] sur_Residual; sur_Residual = nullptr;
}

void PLocAssem_Tet_Wall_Prestress::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("  Prestress generation for CMM type wall surface: \n");
  if(snLocBas == 3)
    SYS_T::commPrint("  FEM: linear tetrahedral/triangle element \n");
  else if(snLocBas == 6)
    SYS_T::commPrint("  FEM: quadratic tetrahedral/triangle element \n");
  else SYS_T::print_fatal("Error: unknown elem type.\n");
  SYS_T::commPrint("  Density rho = %e \n", rho0);
  SYS_T::commPrint("  Dynamic Viscosity mu = %e \n", vis_mu);
  SYS_T::commPrint("  Kienmatic Viscosity nu = %e \n", vis_mu / rho0);
  SYS_T::commPrint("  Note: \n");
  SYS_T::commPrint("  1. Only pressure is applied to generate the wall prestress.\n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
}

void PLocAssem_Tet_Wall_Prestress::Assem_Tangent_Residual_EBC_Wall(
    const double &time, const double &dt,
    const double * const &dot_sol,
    const double * const &sol,
    const double * const &sol_wall_disp,
    FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const double * const &ele_thickness,
    const double * const &ele_youngsmod,
    const double * const &ele_springconst,
    const double * const &ele_dampingconst,
    const double * const &qua_prestress,
    const IQuadPts * const &quad )
{
  element->buildBasis( quad, eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z ); 

  const int dim = 3;

  const double dd_dv = alpha_f * gamma * dt;
  const double dd_du = dd_dv * dd_dv / alpha_m;

  const double curr = time + alpha_f * dt;

  // Global Cauchy stress at all quadrature points
  std::vector<Tensor2_3D> sigma( face_nqp, Tensor2_3D() );
  get_Wall_CauchyStress(sol_wall_disp, element, ele_youngsmod, sigma );

  Zero_sur_Tangent_Residual();

  for(int qua=0; qua<face_nqp; ++qua)
  {
    const std::vector<double> R = element -> get_R( qua );
    
    // For membrane elements, basis function gradients are computed
    // with respect to lamina coords
    const std::vector<double> dR_dxl = element -> get_dR_dx( qua );
    const std::vector<double> dR_dyl = element -> get_dR_dy( qua );

    // Lamina and global stiffness matrices
    double * Kl = new double [ (snLocBas*dim) * (snLocBas*dim) ] {};
    double * Kg = new double [ (snLocBas*dim) * (snLocBas*dim) ] {};

    // Global-to-local rotation matrix Q
    const Tensor2_3D Q = element->get_rotationMatrix(qua);

    double u_t = 0.0, v_t = 0.0, w_t = 0.0;
    double u = 0.0, v = 0.0, w = 0.0, pp = 0.0; 
    double disp_x = 0.0, disp_y = 0.0, disp_z = 0.0;
    double h_w = 0.0, E_w = 0.0, ks_w = 0.0, cs_w = 0.0;
    Vector_3 coor(0.0, 0.0, 0.0);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      u_t += dot_sol[ii*4+1] * R[ii];
      v_t += dot_sol[ii*4+2] * R[ii];
      w_t += dot_sol[ii*4+3] * R[ii];

      pp += sol[ii*4]   * R[ii];
      u  += sol[ii*4+1] * R[ii];
      v  += sol[ii*4+2] * R[ii];
      w  += sol[ii*4+3] * R[ii];

      disp_x += sol_wall_disp[ii*3+0] * R[ii];
      disp_y += sol_wall_disp[ii*3+1] * R[ii];
      disp_z += sol_wall_disp[ii*3+2] * R[ii];

      h_w += ele_thickness[ii] * R[ii];
      E_w += ele_youngsmod[ii] * R[ii];

      ks_w += ele_springconst[ii]  * R[ii];
      cs_w += ele_dampingconst[ii] * R[ii];

      coor.x() += eleCtrlPts_x[ii] * R[ii];
      coor.y() += eleCtrlPts_y[ii] * R[ii];
      coor.z() += eleCtrlPts_z[ii] * R[ii];
    }

    // Body force acting on the wall
    const Vector_3 fw = get_fw(coor, curr);

    double surface_area;
    const Vector_3 n_out = element->get_2d_normal_out(qua, surface_area);

    const double gwts = surface_area * quad->get_qw(qua);

    const double coef = E_w / (1.0 - nu_w * nu_w);

    // Add prestress: convert from Voigt notation (comps 11, 22, 33, 23, 13, 12)
    sigma[qua].xx() += qua_prestress[qua*6];
    sigma[qua].xy() += qua_prestress[qua*6+5];
    sigma[qua].xz() += qua_prestress[qua*6+4];
    sigma[qua].yx() += qua_prestress[qua*6+5];
    sigma[qua].yy() += qua_prestress[qua*6+1];
    sigma[qua].yz() += qua_prestress[qua*6+3];
    sigma[qua].zx() += qua_prestress[qua*6+4];
    sigma[qua].zy() += qua_prestress[qua*6+3];
    sigma[qua].zz() += qua_prestress[qua*6+2];

    // Basis function gradients with respect to global coords
    // dR/dx_{i} = Q_{ji} * dR/dxl_{j}. Note that dR/dzl = 0.0
    std::vector<double> dR_dx(snLocBas, 0.0), dR_dy(snLocBas, 0.0), dR_dz(snLocBas, 0.0);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      dR_dx[ii] = Q.xx() * dR_dxl[ii] + Q.yx() * dR_dyl[ii];
      dR_dy[ii] = Q.xy() * dR_dxl[ii] + Q.yy() * dR_dyl[ii];
      dR_dz[ii] = Q.xz() * dR_dxl[ii] + Q.yz() * dR_dyl[ii];
    }

    // Stiffness tensor in lamina coords
    // Bl^T * D * Bl = Bl_{ki} * D_{kl} * Bl_{lj}
    for(int A=0; A<snLocBas; ++A)
    {
      const double NA_xl = dR_dxl[A], NA_yl = dR_dyl[A];

      for(int B=0; B<snLocBas; ++B)
      {
        const double NB_xl = dR_dxl[B], NB_yl = dR_dyl[B];

        // Momentum-x with respect to u1, u2 
        Kl[(snLocBas*dim)*(A*dim) + (B*dim)]     += coef * ( NA_xl * NB_xl
            + 0.5*(1.0-nu_w) * NA_yl * NB_yl );
        Kl[(snLocBas*dim)*(A*dim) + (B*dim+1)]   += coef * ( nu_w * NA_xl * NB_yl
            + 0.5*(1.0-nu_w) * NA_yl * NB_xl );

        // Momentum-y with respect to u1, u2 
        Kl[(snLocBas*dim)*(A*dim+1) + (B*dim)]   += coef * ( nu_w * NA_yl * NB_xl
            + 0.5*(1.0-nu_w) * NA_xl * NB_yl );
        Kl[(snLocBas*dim)*(A*dim+1) + (B*dim+1)] += coef * ( NA_yl * NB_yl
            + 0.5*(1.0-nu_w) * NA_xl * NB_xl );

        // Momentum-z with respect to u3 
        Kl[(snLocBas*dim)*(A*dim+2) + (B*dim+2)] += coef * 0.5*kappa_w*(1.0-nu_w) * (
            NA_xl * NB_xl + NA_yl * NB_yl );
      }
    }

    // Stiffness tensor in global coords
    // theta^T * Kl * theta, where theta = [Q, 0, 0; 0, Q, 0; 0, 0, Q]
    // or Q^T * Kl_[AB] * Q = Q_{ki} * Kl_[AB]{kl} * Q_{lj}
    for(int A=0; A<snLocBas; ++A)
    {
      for(int B=0; B<snLocBas; ++B)
      {
        for(int ii=0; ii<dim; ++ii)
        {
          for(int jj=0; jj<dim; ++jj)
          {
            // Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(kk,ii) * Kl[ (A*dim+kk)*(snLocBas*dim) + (B*dim+ll) ] * Q(ll,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(snLocBas*dim) + (B*dim+0) ] * Q(0,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(snLocBas*dim) + (B*dim+1) ] * Q(1,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(0,ii) * Kl[ (A*dim+0)*(snLocBas*dim) + (B*dim+2) ] * Q(2,jj);

            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(snLocBas*dim) + (B*dim+0) ] * Q(0,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(snLocBas*dim) + (B*dim+1) ] * Q(1,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(1,ii) * Kl[ (A*dim+1)*(snLocBas*dim) + (B*dim+2) ] * Q(2,jj);

            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(snLocBas*dim) + (B*dim+0) ] * Q(0,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(snLocBas*dim) + (B*dim+1) ] * Q(1,jj);
            Kg[ (snLocBas*dim)*(A*dim+ii) + (B*dim+jj) ] += Q(2,ii) * Kl[ (A*dim+2)*(snLocBas*dim) + (B*dim+2) ] * Q(2,jj);
          }
        }
      }
    }

    for(int A=0; A<snLocBas; ++A)
    {
      const double NA_x = dR_dx[A], NA_y = dR_dy[A], NA_z = dR_dz[A];

      sur_Residual[4*A+1] += gwts * h_w * ( R[A] * rho_w * ( u_t - fw.x() )
          + NA_x * sigma[qua].xx() + NA_y * sigma[qua].xy() + NA_z * sigma[qua].xz() )
          + gwts * R[A] * ( ks_w * disp_x + cs_w * u )
          - gwts * R[A] * pp * n_out.x(); 

      sur_Residual[4*A+2] += gwts * h_w * ( R[A] * rho_w * ( v_t - fw.y() )
          + NA_x * sigma[qua].yx() + NA_y * sigma[qua].yy() + NA_z * sigma[qua].yz() ) 
          + gwts * R[A] * ( ks_w * disp_y + cs_w * v )
          - gwts * R[A] * pp * n_out.y(); 

      sur_Residual[4*A+3] += gwts * h_w * ( R[A] * rho_w * ( w_t - fw.z() )
          + NA_x * sigma[qua].zx() + NA_y * sigma[qua].zy() + NA_z * sigma[qua].zz() ) 
          + gwts * R[A] * ( ks_w * disp_z + cs_w * w )
          - gwts * R[A] * pp * n_out.z(); 

      for(int B=0; B<snLocBas; ++B)
      {
        // Momentum-x with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+1 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim) ] )
            + gwts * R[A] * R[B] * ( dd_du * ks_w + dd_dv * cs_w );

        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+2 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim+1) ] );

        sur_Tangent[ 4*snLocBas*(4*A+1) + 4*B+3 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim) + (B*dim+2) ] );

        // Momentum-y with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+1 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim) ] );

        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+2 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim+1) ] )
            + gwts * R[A] * R[B] * ( dd_du * ks_w + dd_dv * cs_w );

        sur_Tangent[ 4*snLocBas*(4*A+2) + 4*B+3 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+1) + (B*dim+2) ] );

        // Momentum-z with respect to u, v, w
        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+1 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim) ] );

        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+2 ] += gwts * h_w * (
            dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim+1) ] );

        sur_Tangent[ 4*snLocBas*(4*A+3) + 4*B+3 ] += gwts * h_w * (
            alpha_m * rho_w * R[A] * R[B]
            + dd_du * Kg[ (snLocBas*dim)*(A*dim+2) + (B*dim+2) ] )
            + gwts * R[A] * R[B] * ( dd_du * ks_w + dd_dv * cs_w );

      } // end B loop
    } // end A loop

    delete [] Kl; delete [] Kg; Kl = nullptr; Kg = nullptr;

  } // end qua loop

}

void PLocAssem_Tet_Wall_Prestress::get_Wall_CauchyStress(
    const double * const &sol_wall_disp,
    const FEAElement * const &element,
    const double * const &ele_youngsmod,
    std::vector<Tensor2_3D> &sigma )
{
  SYS_T::print_fatal_if( element -> get_numQuapts() != face_nqp, "Error: The element's data structure is incompatible with the face_nqp in the local assembly.\n" );

  const int dim = 3;

  // Lamina displacements
  double * sol_wall_disp_l = new double [ snLocBas * dim ];

  for(int qua=0; qua<face_nqp; ++qua)
  {
    const std::vector<double> R = element -> get_R( qua );
    
    // For membrane elements, basis function gradients are computed
    // with respect to lamina coords
    const std::vector<double> dR_dxl = element -> get_dR_dx( qua );
    const std::vector<double> dR_dyl = element -> get_dR_dy( qua );

    // Global-to-local rotation matrix Q
    const Tensor2_3D Q = element->get_rotationMatrix(qua);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      sol_wall_disp_l[dim*ii]   = sol_wall_disp[dim*ii] * Q.xx() 
        + sol_wall_disp[dim*ii+1] * Q.xy() + sol_wall_disp[dim*ii+2] * Q.xz();

      sol_wall_disp_l[dim*ii+1] = sol_wall_disp[dim*ii] * Q.yx()
        + sol_wall_disp[dim*ii+1] * Q.yy() + sol_wall_disp[dim*ii+2] * Q.yz();

      sol_wall_disp_l[dim*ii+2] = sol_wall_disp[dim*ii] * Q.zx()
        + sol_wall_disp[dim*ii+1] * Q.zy() + sol_wall_disp[dim*ii+2] * Q.zz();
    }

    double E_w = 0.0;
    double u1l_xl = 0.0, u2l_xl = 0.0, u3l_xl = 0.0;
    double u1l_yl = 0.0, u2l_yl = 0.0, u3l_yl = 0.0;

    for(int ii=0; ii<snLocBas; ++ii)
    {
      E_w += ele_youngsmod[ii] * R[ii];

      u1l_xl += sol_wall_disp_l[dim*ii]   * dR_dxl[ii];
      u1l_yl += sol_wall_disp_l[dim*ii]   * dR_dyl[ii];
      u2l_xl += sol_wall_disp_l[dim*ii+1] * dR_dxl[ii];
      u2l_yl += sol_wall_disp_l[dim*ii+1] * dR_dyl[ii];
      u3l_xl += sol_wall_disp_l[dim*ii+2] * dR_dxl[ii];
      u3l_yl += sol_wall_disp_l[dim*ii+2] * dR_dyl[ii];
    }

    const double coef = E_w / (1.0 - nu_w * nu_w);

    // Lamina Cauchy stress
    sigma[qua] = Tensor2_3D(
        u1l_xl + nu_w*u2l_yl,               0.5*(1.0-nu_w) * (u1l_yl + u2l_xl), 0.5*kappa_w*(1.0-nu_w) * u3l_xl,
        0.5*(1.0-nu_w) * (u1l_yl + u2l_xl), nu_w*u1l_xl + u2l_yl,               0.5*kappa_w*(1.0-nu_w) * u3l_yl,
        0.5*kappa_w*(1.0-nu_w) * u3l_xl,    0.5*kappa_w*(1.0-nu_w) * u3l_yl,    0.0 );
    sigma[qua] *= coef;

    // Global Cauchy stress: Q^T * lamina_Cauchy * Q
    sigma[qua].MatRot(Q);
  }

  delete [] sol_wall_disp_l; sol_wall_disp_l = nullptr;
}

// EOF
