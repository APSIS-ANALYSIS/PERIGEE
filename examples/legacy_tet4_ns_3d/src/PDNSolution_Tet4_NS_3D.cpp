#include "PDNSolution_Tet4_NS_3D.hpp"

PDNSolution_Tet4_NS_3D::PDNSolution_Tet4_NS_3D(  const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr, const int &type, 
    const bool &isprint ) : PDNSolution( pNode ), is_print( isprint )
{
  if( pNode->get_dof() != 4 ) SYS_T::print_fatal("Error: PDNSolution_Tet4_NS_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Tet4_NS_3D: No such type of initional condition. \n");
      break;
  }
}


PDNSolution_Tet4_NS_3D::PDNSolution_Tet4_NS_3D(  const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc,
    const double &in_vol_rate_Q,
    const int &type, 
    const bool &isprint ) : PDNSolution( pNode ), is_print( isprint )
{
  if( pNode->get_dof() != 4 ) SYS_T::print_fatal("Error: PDNSolution_Tet4_NS_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode, fNode_ptr );
      break;
    case 1:
      Init_inflow_parabolic_2cm( pNode, fNode_ptr, infbc, in_vol_rate_Q );
      break;
    case 2:
      Init_flow_parabolic( pNode, fNode_ptr, infbc, in_vol_rate_Q );
      break;
    case 3:
      Init_flow_plug( pNode, fNode_ptr, infbc, in_vol_rate_Q );
      break;
    case 5:
      Init_inflow_parabolic_z_0d6cm( pNode, fNode_ptr, infbc, in_vol_rate_Q );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Tet4_NS_3D: No such type of initional condition. \n");
      break;
  }
}


PDNSolution_Tet4_NS_3D::~PDNSolution_Tet4_NS_3D()
{}


void PDNSolution_Tet4_NS_3D::Init_zero(const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();
  
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print){
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = 0.0 \n");
    SYS_T::commPrint("                       velo_y = 0.0 \n");
    SYS_T::commPrint("                       velo_z = 0.0 \n");
  }
}

void PDNSolution_Tet4_NS_3D::Init_inflow_parabolic_2cm(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc,
    const double &vol_rate_Q )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();
  double x, y, r, vel;

  const double pi = MATH_T::PI;
  const double vmax = 2.0 * vol_rate_Q / (pi * 4.0);

  const double out_nx = infbc->get_outvec().x();
  const double out_ny = infbc->get_outvec().y();
  const double out_nz = infbc->get_outvec().z();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  if( infbc->get_Num_LD() > 0 )
  {
    for(int ii=0; ii<nlocalnode; ++ii)
    {
      if( infbc->is_inLDN(pNode_ptr->get_node_loc(ii)) )
      {
        location[0] = pNode_ptr->get_node_loc(ii) * 4;
        location[1] = location[0] + 1;
        location[2] = location[0] + 2;
        location[3] = location[0] + 3;

        x = fNode_ptr->get_ctrlPts_x(ii);
        y = fNode_ptr->get_ctrlPts_y(ii);

        r = sqrt(x*x + y*y);

        vel = vmax * (1.0 - r*r/4.0);

        // time -1 because out_nx points in the outward direction
        value[1] = vel * (-1.0) * out_nx;
        value[2] = vel * (-1.0) * out_ny;
        value[3] = vel * (-1.0) * out_nz;
        VecSetValues(solution, 4, location, value, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print){
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = parabolic \n");
    SYS_T::commPrint("                       velo_y = parabolic \n");
    SYS_T::commPrint("                       velo_z = parabolic \n");
    PetscPrintf(PETSC_COMM_WORLD,
        "                       flow rate %e cm^3/sec \n", vol_rate_Q);
  }
}


void PDNSolution_Tet4_NS_3D::Init_flow_parabolic(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc,
    const double &vol_rate_Q )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  // Maximum speed is 2.0 x flow rate / surface area
  const double vmax = 2.0 * vol_rate_Q / infbc->get_fularea();

  // Get the outward normal vector
  const double out_nx = infbc->get_outvec().x();
  const double out_ny = infbc->get_outvec().y();
  const double out_nz = infbc->get_outvec().z();

  if( infbc->get_Num_LD() > 0 )
  {
    for(int ii=0; ii<nlocalnode; ++ii)
    {
      if( infbc->is_inLDN(pNode_ptr->get_node_loc(ii)) )
      {
        location[0] = pNode_ptr->get_node_loc(ii) * 4;
        location[1] = location[0] + 1;
        location[2] = location[0] + 2;
        location[3] = location[0] + 3;

        const Vector_3 pt = fNode_ptr -> get_ctrlPts_xyz(ii);
        const double r =  infbc -> get_radius( pt );

        // Set zero pressure 
        value[0] = 0.0;

        // The infbc class does not contain the boundary points.
        // So we do not need to worry about setting them as zero. They
        // are not touched by this function.
        const double vel = vmax * (1.0 - r*r);

        // time -1 because out_nx points in the outward direction
        value[1] = vel * (-1.0) * out_nx;
        value[2] = vel * (-1.0) * out_ny;
        value[3] = vel * (-1.0) * out_nz;

        VecSetValues(solution, 4, location, value, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print){
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = parabolic \n");
    SYS_T::commPrint("                       velo_y = parabolic \n");
    SYS_T::commPrint("                       velo_z = parabolic \n");
    PetscPrintf(PETSC_COMM_WORLD,
        "                       flow rate %e (cm^3/sec).\n", vol_rate_Q);
    PetscPrintf(PETSC_COMM_WORLD,
        "                       max speed %e.\n", vmax);
    PetscPrintf(PETSC_COMM_WORLD,
        "                       active area is %e.\n", infbc->get_actarea() );
    PetscPrintf(PETSC_COMM_WORLD,
        "                       full area is %e.\n", infbc->get_fularea() );
    PetscPrintf(PETSC_COMM_WORLD,
        "                       direction [%e %e %e].\n", out_nx, out_ny, out_nz);
  }
}
void PDNSolution_Tet4_NS_3D::Init_flow_plug(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc,
    const double &vol_rate_Q )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();
  double x, y, z, r;

  const double pi = MATH_T::PI;
  const double vel = vol_rate_Q / (pi * 1.95 * 1.95);

  const double out_nx = infbc->get_outvec().x();
  const double out_ny = infbc->get_outvec().y();
  const double out_nz = infbc->get_outvec().z();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    r = sqrt(x*x + y*y);

    if( r < 1.95 && z > 29.9 )
    {
      // time -1 because out_nx points in the outward direction
      value[1] = vel * (-1.0) * out_nx;
      value[2] = vel * (-1.0) * out_ny;
      value[3] = vel * (-1.0) * out_nz;
    }
    else
    {
      value[1] = 0.0;
      value[2] = 0.0;
      value[3] = 0.0;
    }

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print){
    SYS_T::commPrint("===> Initial solution: pres   = 0 \n");
    SYS_T::commPrint("                       velo_x = plug \n");
    SYS_T::commPrint("                       velo_y = plug \n");
    SYS_T::commPrint("                       velo_z = plug \n");
    PetscPrintf(PETSC_COMM_WORLD,
        "                       flow rate %e cm^3/sec \n", vol_rate_Q);
  }
}

void PDNSolution_Tet4_NS_3D::Init_inflow_parabolic_z_0d6cm(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc,
    const double &vol_rate_Q )
{
  int location[4];
  double value[4] = {0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();
  double x, y, rr, vel;

  const double pi = MATH_T::PI;
  const double rd = 0.6;
  const double vmax = 2.0 * vol_rate_Q / (pi * rd*rd);

  const double out_nx = infbc->get_outvec().x();
  const double out_ny = infbc->get_outvec().y();
  const double out_nz = infbc->get_outvec().z();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  if( infbc->get_Num_LD() > 0 )
  {
    for(int ii=0; ii<nlocalnode; ++ii)
    {
      if( infbc->is_inLDN(pNode_ptr->get_node_loc(ii)) )
      {
        location[0] = pNode_ptr->get_node_loc(ii) * 4;
        location[1] = location[0] + 1;
        location[2] = location[0] + 2;
        location[3] = location[0] + 3;

        y = fNode_ptr->get_ctrlPts_y(ii);
        x = fNode_ptr->get_ctrlPts_x(ii);

        rr = x*x + y*y;

        vel = vmax * ( 1.0 - (rr/(rd*rd)) );

        // time -1 because out_nx points in the outward direction
        value[1] = vel * (-1.0) * out_nx;
        value[2] = vel * (-1.0) * out_ny;
        value[3] = vel * (-1.0) * out_nz;
        VecSetValues(solution, 4, location, value, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print){
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = parabolic \n");
    SYS_T::commPrint("                       velo_y = parabolic \n");
    SYS_T::commPrint("                       velo_z = parabolic \n");
    PetscPrintf(PETSC_COMM_WORLD,
        "                       flow rate %e cm^3/sec \n", vol_rate_Q);
  }
}

// EOF
