#include "PLocAssem_Transport_VMS_NS_GenAlpha.hpp"

PLocAssem_Transport_VMS_NS_GenAlpha::PLocAssem_Transport_VMS_NS_GenAlpha(
    const FEType &in_type, const int &in_nqp_v, const double &in_ct,
    const double &in_vis_mu, const TimeMethod_GenAlpha * const &tm_gAlpha,
    std::unique_ptr<IHIModel> in_hi_model )
: elemType(in_type),
  nqpv(in_nqp_v),
  elementv( ElementFactory::createVolElement(elemType, nqpv) ),
  quadv( QuadPtsFactory::createVolQuadrature(elemType, nqpv) ),
  hi_model(std::move(in_hi_model)),
  alpha_f( tm_gAlpha->get_alpha_f() ),
  alpha_m( tm_gAlpha->get_alpha_m() ),
  gamma( tm_gAlpha->get_gamma() ),
  Ct( in_ct ),
  vis_mu( in_vis_mu ),
  nLocBas( elementv->get_nLocBas() ),
  vec_size( nLocBas )
{
  Tangent = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  Zero_Tangent_Residual();
  print_info();
}

PLocAssem_Transport_VMS_NS_GenAlpha::~PLocAssem_Transport_VMS_NS_GenAlpha()
{
  delete [] Tangent; Tangent = nullptr;
  delete [] Residual; Residual = nullptr;
}

void PLocAssem_Transport_VMS_NS_GenAlpha::print_info() const
{
  SYS_T::print_sep_line();
  SYS_T::commPrint("  Weakly coupled scalar transport equation:\n");
  elementv->print_info();
  SYS_T::commPrint("  Spatial: finite element with SUPG/RBVMS stabilization\n");
  SYS_T::commPrint("  Temporal: Generalized-alpha Method\n");
  SYS_T::commPrint("  Transport stress viscosity mu = %e \n", vis_mu);
  if(hi_model == nullptr)
  {
    SYS_T::commPrint("  Source model: none (pure transport / zero source)\n");
  }
  else
  {
    SYS_T::commPrint("  Source model: %s \n", hi_model->get_model_name().c_str());
    hi_model->print_info();
  }
  SYS_T::print_sep_line();
}

void PLocAssem_Transport_VMS_NS_GenAlpha::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 0.0;
}

void PLocAssem_Transport_VMS_NS_GenAlpha::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii) Residual[ii] = 0.0;
}

void PLocAssem_Transport_VMS_NS_GenAlpha::Assem_Estimate()
{
  for(int ii=0; ii<vec_size*vec_size; ++ii) Tangent[ii] = 1.0;
}

double PLocAssem_Transport_VMS_NS_GenAlpha::get_tau(
    const double &dt, const std::array<double, 9> &dxi_dx,
    const double &ux, const double &uy, const double &uz ) const
{
  const double g11 = dxi_dx[0]*dxi_dx[0] + dxi_dx[3]*dxi_dx[3] + dxi_dx[6]*dxi_dx[6];
  const double g12 = dxi_dx[0]*dxi_dx[1] + dxi_dx[3]*dxi_dx[4] + dxi_dx[6]*dxi_dx[7];
  const double g13 = dxi_dx[0]*dxi_dx[2] + dxi_dx[3]*dxi_dx[5] + dxi_dx[6]*dxi_dx[8];
  const double g22 = dxi_dx[1]*dxi_dx[1] + dxi_dx[4]*dxi_dx[4] + dxi_dx[7]*dxi_dx[7];
  const double g23 = dxi_dx[1]*dxi_dx[2] + dxi_dx[4]*dxi_dx[5] + dxi_dx[7]*dxi_dx[8];
  const double g33 = dxi_dx[2]*dxi_dx[2] + dxi_dx[5]*dxi_dx[5] + dxi_dx[8]*dxi_dx[8];

  const double ugu = ux * (g11*ux + g12*uy + g13*uz)
                   + uy * (g12*ux + g22*uy + g23*uz)
                   + uz * (g13*ux + g23*uy + g33*uz);

  return 1.0 / std::sqrt( std::pow(2.0 / dt, 2.0) + Ct * ugu );
}

double PLocAssem_Transport_VMS_NS_GenAlpha::get_scalar_shear_stress(
    const double &u_x, const double &u_y, const double &u_z,
    const double &v_x, const double &v_y, const double &v_z,
    const double &w_x, const double &w_y, const double &w_z ) const
{
  const double two_mu = 2.0 * vis_mu;
  const SymmTensor2_3D tau(
      two_mu * u_x,
      two_mu * v_y,
      two_mu * w_z,
      vis_mu * (v_z + w_y),
      vis_mu * (u_z + w_x),
      vis_mu * (u_y + v_x) );

  return std::sqrt(0.5 * tau.MatContraction());
}

void PLocAssem_Transport_VMS_NS_GenAlpha::Assem_Residual(
    const double &time, const double &dt,
    const double * const &dot_phi,
    const double * const &phi,
    const double * const &flow_sol,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z )
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );
  Zero_Residual();
  UNUSED(time);

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqpv; ++qua)
  {
    double phi_t = 0.0, phi_x = 0.0, phi_y = 0.0, phi_z = 0.0;
    double ux = 0.0, uy = 0.0, uz = 0.0;
    double u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w_x = 0.0, w_y = 0.0, w_z = 0.0;

    elementv->get_R_gradR(qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0]);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      phi_t += dot_phi[ii] * R[ii];
      phi_x += phi[ii] * dR_dx[ii];
      phi_y += phi[ii] * dR_dy[ii];
      phi_z += phi[ii] * dR_dz[ii];

      ux += flow_sol[4*ii+1] * R[ii];
      uy += flow_sol[4*ii+2] * R[ii];
      uz += flow_sol[4*ii+3] * R[ii];
      u_x += flow_sol[4*ii+1] * dR_dx[ii];
      u_y += flow_sol[4*ii+1] * dR_dy[ii];
      u_z += flow_sol[4*ii+1] * dR_dz[ii];
      v_x += flow_sol[4*ii+2] * dR_dx[ii];
      v_y += flow_sol[4*ii+2] * dR_dy[ii];
      v_z += flow_sol[4*ii+2] * dR_dz[ii];
      w_x += flow_sol[4*ii+3] * dR_dx[ii];
      w_y += flow_sol[4*ii+3] * dR_dy[ii];
      w_z += flow_sol[4*ii+3] * dR_dz[ii];
    }

    const double adv = ux * phi_x + uy * phi_y + uz * phi_z;
    const double tau_s = get_scalar_shear_stress(u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z);
    const double source = get_source(tau_s);
    const double strong = phi_t + adv - source;
    const double tau = get_tau(dt, elementv->get_invJacobian(qua), ux, uy, uz);
    const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);

    for(int A=0; A<nLocBas; ++A)
    {
      const double u_grad_NA = ux * dR_dx[A] + uy * dR_dy[A] + uz * dR_dz[A];
      Residual[A] += gwts * ( R[A] * strong + tau * u_grad_NA * strong );
    }
  }
}

void PLocAssem_Transport_VMS_NS_GenAlpha::Assem_Tangent_Residual(
    const double &time, const double &dt,
    const double * const &dot_phi,
    const double * const &phi,
    const double * const &flow_sol,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z )
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );
  Zero_Tangent_Residual();
  UNUSED(time);

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqpv; ++qua)
  {
    double phi_t = 0.0, phi_x = 0.0, phi_y = 0.0, phi_z = 0.0;
    double ux = 0.0, uy = 0.0, uz = 0.0;
    double u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w_x = 0.0, w_y = 0.0, w_z = 0.0;

    elementv->get_R_gradR(qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0]);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      phi_t += dot_phi[ii] * R[ii];
      phi_x += phi[ii] * dR_dx[ii];
      phi_y += phi[ii] * dR_dy[ii];
      phi_z += phi[ii] * dR_dz[ii];

      ux += flow_sol[4*ii+1] * R[ii];
      uy += flow_sol[4*ii+2] * R[ii];
      uz += flow_sol[4*ii+3] * R[ii];
      u_x += flow_sol[4*ii+1] * dR_dx[ii];
      u_y += flow_sol[4*ii+1] * dR_dy[ii];
      u_z += flow_sol[4*ii+1] * dR_dz[ii];
      v_x += flow_sol[4*ii+2] * dR_dx[ii];
      v_y += flow_sol[4*ii+2] * dR_dy[ii];
      v_z += flow_sol[4*ii+2] * dR_dz[ii];
      w_x += flow_sol[4*ii+3] * dR_dx[ii];
      w_y += flow_sol[4*ii+3] * dR_dy[ii];
      w_z += flow_sol[4*ii+3] * dR_dz[ii];
    }

    const double adv = ux * phi_x + uy * phi_y + uz * phi_z;
    const double tau_s = get_scalar_shear_stress(u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z);
    const double source = get_source(tau_s);
    const double strong = phi_t + adv - source;
    const double tau = get_tau(dt, elementv->get_invJacobian(qua), ux, uy, uz);
    const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);

    for(int A=0; A<nLocBas; ++A)
    {
      const double NA = R[A];
      const double u_grad_NA = ux * dR_dx[A] + uy * dR_dy[A] + uz * dR_dz[A];

      Residual[A] += gwts * ( NA * strong + tau * u_grad_NA * strong );

      for(int B=0; B<nLocBas; ++B)
      {
        const double NB = R[B];
        const double u_grad_NB = ux * dR_dx[B] + uy * dR_dy[B] + uz * dR_dz[B];

        Tangent[nLocBas * A + B] += gwts * (
            alpha_m * NA * NB
          + alpha_f * gamma * dt * NA * u_grad_NB
          + tau * u_grad_NA * ( alpha_m * NB + alpha_f * gamma * dt * u_grad_NB ) );
      }
    }
  }
}

void PLocAssem_Transport_VMS_NS_GenAlpha::Assem_Mass_Residual(
    const double * const &phi,
    const double * const &flow_sol,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z )
{
  elementv->buildBasis( quadv.get(), eleCtrlPts_x, eleCtrlPts_y, eleCtrlPts_z );
  Zero_Tangent_Residual();

  std::vector<double> R(nLocBas, 0.0), dR_dx(nLocBas, 0.0), dR_dy(nLocBas, 0.0), dR_dz(nLocBas, 0.0);

  for(int qua=0; qua<nqpv; ++qua)
  {
    double phi_x = 0.0, phi_y = 0.0, phi_z = 0.0;
    double ux = 0.0, uy = 0.0, uz = 0.0;
    double u_x = 0.0, u_y = 0.0, u_z = 0.0;
    double v_x = 0.0, v_y = 0.0, v_z = 0.0;
    double w_x = 0.0, w_y = 0.0, w_z = 0.0;

    elementv->get_R_gradR(qua, &R[0], &dR_dx[0], &dR_dy[0], &dR_dz[0]);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      phi_x += phi[ii] * dR_dx[ii];
      phi_y += phi[ii] * dR_dy[ii];
      phi_z += phi[ii] * dR_dz[ii];

      ux += flow_sol[4*ii+1] * R[ii];
      uy += flow_sol[4*ii+2] * R[ii];
      uz += flow_sol[4*ii+3] * R[ii];
      u_x += flow_sol[4*ii+1] * dR_dx[ii];
      u_y += flow_sol[4*ii+1] * dR_dy[ii];
      u_z += flow_sol[4*ii+1] * dR_dz[ii];
      v_x += flow_sol[4*ii+2] * dR_dx[ii];
      v_y += flow_sol[4*ii+2] * dR_dy[ii];
      v_z += flow_sol[4*ii+2] * dR_dz[ii];
      w_x += flow_sol[4*ii+3] * dR_dx[ii];
      w_y += flow_sol[4*ii+3] * dR_dy[ii];
      w_z += flow_sol[4*ii+3] * dR_dz[ii];
    }

    const double adv = ux * phi_x + uy * phi_y + uz * phi_z;
    const double tau_s = get_scalar_shear_stress(u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z);
    const double source = get_source(tau_s);
    const double strong = adv - source;
    const double tau = get_tau(1.0, elementv->get_invJacobian(qua), ux, uy, uz);
    const double gwts = elementv->get_detJac(qua) * quadv->get_qw(qua);

    for(int A=0; A<nLocBas; ++A)
    {
      const double u_grad_NA = ux * dR_dx[A] + uy * dR_dy[A] + uz * dR_dz[A];
      Residual[A] += gwts * ( R[A] * strong + tau * u_grad_NA * strong );

      for(int B=0; B<nLocBas; ++B)
        Tangent[nLocBas * A + B] += gwts * R[A] * R[B];
    }
  }
}
