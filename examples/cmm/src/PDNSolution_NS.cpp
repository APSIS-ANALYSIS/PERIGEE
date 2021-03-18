#include "PDNSolution_NS.hpp"

PDNSolution_NS::PDNSolution_NS( 
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print(isprint)
{
  if( pNode->get_dof() != 4 ) SYS_T::print_fatal("Error: PDNSolution_NS : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 1:
      Init_flow_parabolic( pNode, fNode_ptr, infbc );
      break;
    case 2:
      Init_pipe_parabolic( pNode, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_NS: No such type of initial condition.\n");
      break;
  }
}


PDNSolution_NS::PDNSolution_NS( 
    const APart_Node * const &pNode,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print( isprint )
{
  if( pNode->get_dof() != 4 ) SYS_T::print_fatal("Error: PDNSolution_NS : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_NS : No such type of initial condition.\n");
      break;
  }
}


// ==== WOMERSLEY CHANGES BEGIN ====
PDNSolution_NS::PDNSolution_NS( 
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const double &rho,
    const double &vis_mu,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print( isprint )
{
  if( pNode->get_dof() != 4 ) SYS_T::print_fatal("Error: PDNSolution_NS : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 3:
      Init_womersley( pNode, fNode_ptr, rho, vis_mu );
      break;
    case 4:
      Init_womersley_dot( pNode, fNode_ptr, rho, vis_mu );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_NS : No such type of initial condition.\n");
      break;
  }
}
// ==== WOMERSLEY CHANGES END ====


PDNSolution_NS::~PDNSolution_NS()
{}


void PDNSolution_NS::Init_zero(const APart_Node * const &pNode_ptr)
{
  int location[4];
  const double value[4] = {0.0, 0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  if( is_print )
  {
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = 0.0 \n");
    SYS_T::commPrint("                       velo_y = 0.0 \n");
    SYS_T::commPrint("                       velo_z = 0.0 \n");
  }
}


void PDNSolution_NS::Init_flow_parabolic(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  // Maximum speed formula is 
  //             2.0 x flow rate (1.0) / surface area
  // Here I use the unit flow rate, and the actual flow rate is adjusted
  // based on the CVFlowRate class.
  const double vmax = 2.0 / infbc->get_fularea();

  const double out_nx = infbc->get_outvec(0);
  const double out_ny = infbc->get_outvec(1);
  const double out_nz = infbc->get_outvec(2);

  // If there are inflow nodes, set their value to be parabolic flow
  if( infbc->get_Num_LD() > 0)
  {
    for(int ii=0; ii<nlocalnode; ++ii)
    {
      if( infbc->is_inLDN(pNode_ptr->get_node_loc(ii)) )
      {
        location[0] = pNode_ptr->get_node_loc(ii) * 4;
        location[1] = location[0] + 1;
        location[2] = location[0] + 2;
        location[3] = location[0] + 3;

        const double x = fNode_ptr->get_ctrlPts_x(ii);
        const double y = fNode_ptr->get_ctrlPts_y(ii);
        const double z = fNode_ptr->get_ctrlPts_z(ii);

        const double r = infbc->get_radius(x,y,z);
        const double vel = vmax * (1.0 - r*r);

        // -1.0 is multiplied to make the flow direction inward
        value[1] = vel * (-1.0) * out_nx;
        value[2] = vel * (-1.0) * out_ny;
        value[3] = vel * (-1.0) * out_nz;

        VecSetValues(solution, 4, location, value, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = parabolic \n");
    SYS_T::commPrint("                       velo_y = parabolic \n");
    SYS_T::commPrint("                       velo_z = parabolic \n");
    SYS_T::commPrint("                       flow rate 1.0 .\n");
    SYS_T::commPrint("                       max speed %e.\n", vmax);
    SYS_T::commPrint("                       active area is %e.\n", infbc->get_actarea() );
    SYS_T::commPrint("                       full area is %e.\n", infbc->get_fularea() );
    SYS_T::commPrint("                       direction [%e %e %e].\n", out_nx, out_ny, out_nz);
  }
}


void PDNSolution_NS::Init_pipe_parabolic(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  // Maximum speed formula is 
  //             2.0 x flow rate (1.0) / surface area
  // Here I use the unit flow rate, and the actual flow rate is adjusted
  // based on the CVFlowRate class.
  const double vmax = 2.0 * 100.0 / (3.14 * 4);

  const double out_nx = 0.0; 
  const double out_ny = 0.0;
  const double out_nz = 1.0;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    const double x = fNode_ptr->get_ctrlPts_x(ii);
    const double y = fNode_ptr->get_ctrlPts_y(ii);
    //const double z = fNode_ptr->get_ctrlPts_z(ii);

    const double r = std::sqrt(x*x + y*y) / 2.0;
    const double vel = vmax * (1.0 - r*r);

    // -1.0 is multiplied to make the flow direction inward
    value[1] = vel * (-1.0) * out_nx;
    value[2] = vel * (-1.0) * out_ny;
    value[3] = vel * (-1.0) * out_nz;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = parabolic \n");
    SYS_T::commPrint("                       velo_y = parabolic \n");
    SYS_T::commPrint("                       velo_z = parabolic \n");
    SYS_T::commPrint("                       flow rate 100.0 .\n");
    SYS_T::commPrint("                       direction [%e %e %e].\n", out_nx, out_ny, out_nz);
  }
}


// ==== WOMERSLEY CHANGES BEGIN ====
void PDNSolution_NS::Init_womersley(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const double &rho,
    const double &vis_mu )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  const double R     = 0.3;                                                  // pipe radius
  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);
  const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);
  const auto Omega   = std::sqrt(rho * omega / vis_mu) * R;                  // womersley number 
  const auto Lambda  = i1_1d5 * Omega;

  const double k0 = -21.0469;                                                // mean pressure gradient
  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed
  const std::complex<double> G1(0.829733473284180,      -0.374935589823809); // elasticity factor

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    const double x  = fNode_ptr->get_ctrlPts_x(ii);
    const double y  = fNode_ptr->get_ctrlPts_y(ii);
    const double z  = fNode_ptr->get_ctrlPts_z(ii);
    const double r  = std::sqrt(x*x + y*y);
    const auto   xi = Lambda * r / R;

    const auto bes0_xi     = sp_bessel::besselJ(0, xi);
    const auto bes1_xi     = sp_bessel::besselJ(1, xi);
    const auto bes0_Lambda = sp_bessel::besselJ(0, Lambda);

    // pressure
    const double pres = k0 * z + std::real( B1 * exp(-i1*omega*z/c1) );

    // radial velo
    const double u = std::real( i1 * omega * R * B1 / ( 2.0 * rho * c1 * c1 )
        * ( r / R - 2.0 * G1 * bes1_xi / (Lambda * bes0_Lambda) ) * exp(-i1*omega*z/c1) );

    // axial velocity
    const double w = k0 * (x*x + y*y - R*R) / (4.0*vis_mu)
        + std::real( B1 / (rho * c1) * (1.0 - G1 * bes0_xi / bes0_Lambda) * exp(-i1*omega*z/c1) );

    value[0] = pres;
    value[1] = u;
    value[2] = u;
    value[3] = w;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();
}


void PDNSolution_NS::Init_womersley_dot(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const double &rho,
    const double &vis_mu )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  const double R     = 0.3;                                                  // pipe radius
  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);
  const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);
  const auto Omega   = std::sqrt(rho * omega / vis_mu) * R;                  // womersley number 
  const auto Lambda  = i1_1d5 * Omega;

  // const double k0 = -21.0469;                                                // mean pressure gradient
  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed
  const std::complex<double> G1(0.829733473284180,      -0.374935589823809); // elasticity factor

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    const double x  = fNode_ptr->get_ctrlPts_x(ii);
    const double y  = fNode_ptr->get_ctrlPts_y(ii);
    const double z  = fNode_ptr->get_ctrlPts_z(ii);
    const double r  = std::sqrt(x*x + y*y);
    const auto   xi = Lambda * r / R;

    const auto bes0_xi     = sp_bessel::besselJ(0, xi);
    const auto bes1_xi     = sp_bessel::besselJ(1, xi);
    const auto bes0_Lambda = sp_bessel::besselJ(0, Lambda);

    // dot pressure
    const double dot_pres = std::real( i1 * omega * B1 * exp(-i1*omega*z/c1) );

    // dot radial velocity
    const double dot_u = std::real( -omega * omega * R * B1 / ( 2.0 * rho * c1 * c1 )
        * ( r / R - 2.0 * G1 * bes1_xi / (Lambda * bes0_Lambda) ) * exp(-i1*omega*z/c1) );

    // dot axial velocity
    const double dot_w = std::real( i1 * omega * B1 / (rho * c1)
        * (1.0 - G1 * bes0_xi / bes0_Lambda) * exp(-i1*omega*z/c1) );

    value[0] = dot_pres; 
    value[1] = dot_u;
    value[2] = dot_u;
    value[3] = dot_w;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();
}

// EOF
