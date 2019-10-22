#include "PDNSolution_Mixed_U_Hyperelastic_3D.hpp"

PDNSolution_Mixed_U_Hyperelastic_3D::PDNSolution_Mixed_U_Hyperelastic_3D(
    const class APart_Node * const &pNode,
    const class IPLocAssem * const &locassem_ptr,
    const class FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution( pNode )
{
  if( pNode->get_dof() != 7 ) SYS_T::print_fatal("Error: PDNSolution_Mixed_U_HyperElastic_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case -3:
      Init_test_3( pNode, locassem_ptr, fNode_ptr );
      break;
    case -2:
      Init_test_2( pNode, locassem_ptr, fNode_ptr );
      break;
    case -1:
      Init_test_1( pNode, locassem_ptr, fNode_ptr );
      break;
    case 0:
      Init_zero( pNode, locassem_ptr, fNode_ptr );
      break;
    case 1:
      Init_zero_u_p_one_v( pNode, locassem_ptr, fNode_ptr );
      break;
    case 2:
      Init_vx_linear( pNode, locassem_ptr, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Mixed_U_HyperElastic_3D: No such type of initional condition. \n");
      break;
  }
}


PDNSolution_Mixed_U_Hyperelastic_3D::~PDNSolution_Mixed_U_Hyperelastic_3D()
{}


void PDNSolution_Mixed_U_Hyperelastic_3D::Init_test_1(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    value[0] = 1.0 * location[0];
    value[1] = 1.0 * location[1];
    value[2] = 1.0 * location[2];
    value[3] = 1.0 * location[3];
    value[4] = 1.0 * location[4];
    value[5] = 1.0 * location[5];
    value[6] = 1.0 * location[6];
    
    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = test \n");
  SYS_T::commPrint("                       disp_y = test \n");
  SYS_T::commPrint("                       disp_z = test \n");
  SYS_T::commPrint("                       pres   = test \n");
  SYS_T::commPrint("                       velo_x = test \n");
  SYS_T::commPrint("                       velo_y = test \n");
  SYS_T::commPrint("                       velo_z = test \n");
}


void PDNSolution_Mixed_U_Hyperelastic_3D::Init_test_2(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    value[0] = SYS_T::gen_randomD_open(0.0, 1.0); 
    value[1] = SYS_T::gen_randomD_open(0.0, 1.0);
    value[2] = SYS_T::gen_randomD_open(0.0, 1.0);
    value[3] = SYS_T::gen_randomD_open(0.0, 1.0);
    value[4] = SYS_T::gen_randomD_open(0.0, 1.0);
    value[5] = SYS_T::gen_randomD_open(0.0, 1.0);
    value[6] = SYS_T::gen_randomD_open(0.0, 1.0);
    
    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = random \n");
  SYS_T::commPrint("                       disp_y = random \n");
  SYS_T::commPrint("                       disp_z = random \n");
  SYS_T::commPrint("                       pres   = random \n");
  SYS_T::commPrint("                       velo_x = random \n");
  SYS_T::commPrint("                       velo_y = random \n");
  SYS_T::commPrint("                       velo_z = random \n");
}


void PDNSolution_Mixed_U_Hyperelastic_3D::Init_test_3(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();
  double x,y,z;
  const double pi = MATH_T::PI;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    value[4] = sin(pi*x) * sin(pi*y) * sin(pi*z);
    value[5] = value[4];
    value[6] = value[4];
   
    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }
  
  Assembly_GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
  SYS_T::commPrint("                       pres   = 0.0 \n");
  SYS_T::commPrint("                       velo_x = sin(pi*x) sin(pi*y) sin(pi*z)\n");
  SYS_T::commPrint("                       velo_y = sin(pi*x) sin(pi*y) sin(pi*z)\n");
  SYS_T::commPrint("                       velo_z = sin(pi*x) sin(pi*y) sin(pi*z)\n");
}


void PDNSolution_Mixed_U_Hyperelastic_3D::Init_zero(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
  SYS_T::commPrint("                       pres   = 0.0 \n");
  SYS_T::commPrint("                       velo_x = 0.0 \n");
  SYS_T::commPrint("                       velo_y = 0.0 \n");
  SYS_T::commPrint("                       velo_z = 0.0 \n");
}


void PDNSolution_Mixed_U_Hyperelastic_3D::Init_zero_u_p_one_v(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, MATH_T::PI, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
  SYS_T::commPrint("                       pres   = 0.0 \n");
  SYS_T::commPrint("                       velo_x = PI  \n");
  SYS_T::commPrint("                       velo_y = 0.0 \n");
  SYS_T::commPrint("                       velo_z = 0.0 \n");
}


void PDNSolution_Mixed_U_Hyperelastic_3D::Init_vx_linear(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    double z = fNode_ptr->get_ctrlPts_z(ii);
    
    value[4] = z * 10.0 / 6.0;
    
    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
  SYS_T::commPrint("                       disp_y = 0.0 \n");
  SYS_T::commPrint("                       disp_z = 0.0 \n");
  SYS_T::commPrint("                       pres   = 0.0 \n");
  SYS_T::commPrint("                       velo_x = 10z/6  \n");
  SYS_T::commPrint("                       velo_y = 0.0 \n");
  SYS_T::commPrint("                       velo_z = 0.0 \n");
}

// EOF
