#include "PDNSolution_U_Mixed_3D.hpp"

PDNSolution_U_Mixed_3D::PDNSolution_U_Mixed_3D(
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode, 3 ), is_print( isprint )
{
  if( pNode->get_dof() != 7 ) SYS_T::print_fatal("Error: PDNSolution_U_Mixed_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case -1:
      Init_test_1( pNode, fNode_ptr );
      break;
    case -2:
      Init_test_2( pNode, fNode_ptr );
      break;
    case 0:
      Init_zero( pNode, fNode_ptr );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_U_Mixed_3D: No such type of initional condition. \n");
      break;
  }
}

PDNSolution_U_Mixed_3D::~PDNSolution_U_Mixed_3D()
{}

void PDNSolution_U_Mixed_3D::Init_test_1(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[3];
  double value[3] = {0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    value[0] = location[0] * 1.0;
    value[1] = location[1] * 1.0;
    value[2] = location[2] * 1.0;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: disp_x = test \n");
    SYS_T::commPrint("                       disp_y = test \n");
    SYS_T::commPrint("                       disp_z = test \n");
  }
}

void PDNSolution_U_Mixed_3D::Init_test_2(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[3];
  double value[3] = {0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    value[0] = SYS_T::gen_randomD_open(0.0, 1.0); 
    value[1] = SYS_T::gen_randomD_open(0.0, 1.0); 
    value[2] = SYS_T::gen_randomD_open(0.0, 1.0); 

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if( is_print )
  {
    SYS_T::commPrint("===> Initial solution: disp_x = random \n");
    SYS_T::commPrint("                       disp_y = random \n");
    SYS_T::commPrint("                       disp_z = random \n");
  }
}

void PDNSolution_U_Mixed_3D::Init_zero(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[3];
  const double value[3] = {0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 3;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  Assembly_GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
    SYS_T::commPrint("                       disp_y = 0.0 \n");
    SYS_T::commPrint("                       disp_z = 0.0 \n");
  }
}

// EOF
