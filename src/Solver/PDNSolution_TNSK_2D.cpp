#include "PDNSolution_TNSK_2D.hpp"

PDNSolution_TNSK_2D::PDNSolution_TNSK_2D( 
    const class APart_Node * const &pNode,
    const class IALocal_BC * const &lbc,
    const class IPLocAssem * const &locassem_ptr,
    const class FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution(pNode)
{
  switch (type)
  {
    case 0:
      Init_uniformRhoThe_zeroVel( pNode, locassem_ptr, fNode_ptr );
      break;
    case 1:
      Init_twoBubble_zeroVel_The0d85( pNode, locassem_ptr, fNode_ptr );
      break;
    case 2:
      Init_twoBubble_zeroVel_The0d95( pNode, locassem_ptr, fNode_ptr );
      break;
    case 3:
      Init_randRho_zeroVel( pNode, locassem_ptr, fNode_ptr );
      break;
    case 4:
      Init_1bub_heatbc( pNode, locassem_ptr, fNode_ptr, lbc );
      break;
    case 5:
      Init_1bub_cBothTop( pNode, locassem_ptr, fNode_ptr );
      break;
    case 6:
      Init_boiling( pNode, locassem_ptr, fNode_ptr );
      break;
    case 7:
      Init_RB_convection( pNode, locassem_ptr, fNode_ptr );
      break;
    case 8:
      Init_RB_convection_2( pNode, locassem_ptr, fNode_ptr );
      break;
    case 9:
      Init_thermal_spinodal_decomp( pNode, locassem_ptr, fNode_ptr, lbc );
      break;
    case 10:
      Init_twophase_boiling( pNode, locassem_ptr, fNode_ptr );
      break;
    case 11:
      Init_twophase_boiling_rec( pNode, locassem_ptr, fNode_ptr );
      break;
    case 12:
      Init_multibub_cBothTop( pNode, locassem_ptr, fNode_ptr );
      break;
    case 13:
      Init_1bub_cBothTop_2( pNode, locassem_ptr, fNode_ptr );
      break;
    case 14:
      Init_twophase_boiling_rec_new( pNode, locassem_ptr, fNode_ptr );
      break;
    case 15:
      Init_twophase_boiling_sin( pNode, locassem_ptr, fNode_ptr );
      break;
    default:
      SYS_T::commPrint("ERROR: PDNSolution_TNSK_2D: No such type of initial solution. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}



PDNSolution_TNSK_2D::PDNSolution_TNSK_2D( 
    const class APart_Node * const &pNode,
    const class ALocal_NodalBC * const &lnbc,
    const class IPLocAssem * const &locassem_ptr,
    const class FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution(pNode)
{
  switch (type)
  {
    case -1:
      Init_full_random( pNode );
      break;
    case 0:
      Init_uniformRhoThe_zeroVel( pNode, locassem_ptr, fNode_ptr );
      break;
    case 1:
      Init_twoBubble_zeroVel_The0d85( pNode, locassem_ptr, fNode_ptr );
      break;
    case 2:
      Init_twoBubble_zeroVel_The0d95( pNode, locassem_ptr, fNode_ptr );
      break;
    case 3:
      Init_randRho_zeroVel( pNode, locassem_ptr, fNode_ptr );
      break;
    case 7:
      Init_RB_convection( pNode, locassem_ptr, fNode_ptr );
      break;
    case 11:
      Init_twophase_boiling_rec( pNode, locassem_ptr, fNode_ptr );
      break;
    case 15:
      Init_twophase_boiling_sin( pNode, locassem_ptr, fNode_ptr );
      break;
    case 16:
      Init_thermal_spinodal_decomp_2( pNode, locassem_ptr, fNode_ptr, lnbc );
      break;
    case 17:
      Init_newRB( pNode, locassem_ptr, fNode_ptr );
      break;
    case 18:
      Init_twophase_boiling_wn( pNode, locassem_ptr, fNode_ptr );
      break;
    case 19:
      Init_twophase_boiling_linearRand_temp( pNode, locassem_ptr, fNode_ptr );
      break;
    case 20:
      Init_twophase_boiling_homogeneous( pNode, locassem_ptr, fNode_ptr );
      break;
    case 21:
      Init_twophase_boiling_homogeneous_nopert( pNode, locassem_ptr, fNode_ptr );
      break;
    case 22:
      Init_twophase_boiling_homogeneous_nopert_0d78( pNode, locassem_ptr, fNode_ptr );
      break;
    case 23:
      Init_twophase_boiling_homogeneous_nopert_0d78_square( pNode, locassem_ptr, fNode_ptr );
      break;
    default:
      SYS_T::commPrint("Error: PDNSolution_TNSK_2D: No such type of initial solution. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}


PDNSolution_TNSK_2D::~PDNSolution_TNSK_2D()
{}



void PDNSolution_TNSK_2D::Init_full_random( const APart_Node * const &pNode_ptr )
{
  int location[5];
  double value[5];
  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    value[0] = SYS_T::gen_randomD_open(0,1);
    value[1] = SYS_T::gen_randomD_open(1,2);
    value[2] = SYS_T::gen_randomD_open(2,3);
    value[3] = SYS_T::gen_randomD_open(3,4);
    value[4] = SYS_T::gen_randomD_open(4,5);
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: fully random entries. \n");
}



void PDNSolution_TNSK_2D::Init_uniformRhoThe_zeroVel( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.85;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    value[0] = 0.6;
    value[3] = uniform_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = 0.6, u = 0.0, v = 0.0, theta = 0.85. \n");
}


void PDNSolution_TNSK_2D::Init_twoBubble_zeroVel_The0d85( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();

  double dist1, dist2, x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.85;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    dist1 = sqrt((x-0.4)*(x-0.4) + (y-0.5)*(y-0.5));
    dist2 = sqrt((x-0.75)*(x-0.75) + (y-0.5)*(y-0.5));

    value[0] = 0.1 + 0.25 * (tanh((dist1 - 0.25)/(2*Ca))
                +tanh((dist2-0.1)/(2*Ca)));

    value[3] = uniform_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = two bubble, u = 0.0, v = 0.0, theta = 0.85. \n");
}


void PDNSolution_TNSK_2D::Init_twoBubble_zeroVel_The0d95( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();

  double dist1, dist2, x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.95;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    dist1 = sqrt((x-0.29)*(x-0.29) + (y-0.5)*(y-0.5));
    dist2 = sqrt((x-0.71)*(x-0.71) + (y-0.5)*(y-0.5));

    value[0] = 0.1 + 0.25 * (tanh((dist1 - 0.2)/(2*Ca))
                +tanh((dist2-0.2)/(2*Ca)));

    value[3] = uniform_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = two bubble, u = 0.0, v = 0.0, theta = 0.95. \n");
}


void PDNSolution_TNSK_2D::Init_randRho_zeroVel( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.85;

  const double bot_temp = -1.0 / 0.90;
  const double top_temp = -1.0 / 0.80;

  srand(time(NULL));
  double discrete_rand, random_num;
  
  double y;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    y = fNode_ptr->get_ctrlPts_y(ii);
    
    discrete_rand = rand() % 1001;

    value[0] = 0.35 + (double) discrete_rand * 1.0e-6 - 0.0005;

    value[3] = uniform_temp;
    
    discrete_rand = rand() % 1001;
    random_num = ( (double) discrete_rand * 1.0e-4 - 5.0e-2 ) * 1.0;

    // no random perturbation on temperature bc values
    random_num = 0.0;

    if(y<1.0e-5)
      value[3] = bot_temp + random_num;
    
    if(y>0.99999)
      value[3] = top_temp + random_num;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = random, u = 0.0, v = 0.0, ");
  SYS_T::commPrint(" theta is 0.9 in bot 0.8 in top with randomness. \n ");
}


void PDNSolution_TNSK_2D::Init_1bub_heatbc( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr,
    const IALocal_BC * const &lbc )
{
  const double Ca = locassem_ptr->get_model_para_1();

  double dist1, x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.85;
  const double boundary_temp = -1.0 / 0.95;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));

    value[0] = 0.3545 + 0.2479 * tanh((dist1 - 0.25)/(2*Ca));

    value[3] = uniform_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  
  const int num = lbc->get_Num_LD(3); // get the dirichlet nodes for temperature
  int * bc_index = new int [num];
  double * bc_value = new double [num];

  for(int ii=0; ii<num; ++ii)
  {
    bc_index[ii] = lbc->get_LDN(3, ii) * 5 + 3;
    bc_value[ii] = boundary_temp;
  }
  
  VecSetValues(solution, num, bc_index, bc_value, INSERT_VALUES);
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
 
  delete [] bc_index; delete [] bc_value; 
  
  SYS_T::commPrint("===> Initial solution: rho = one bubble, u = 0.0, v = 0.0, ");
  PetscPrintf(PETSC_COMM_WORLD, 
      "theta = %e inside and %e on boundary \n", -1.0/uniform_temp, -1.0/boundary_temp);
}


void PDNSolution_TNSK_2D::Init_1bub_cBothTop( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();

  double dist1, x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();
  
  const double temp_u = -1.0 / 0.85;
  const double temp_h = -1.0 / 0.90;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));

    value[0] = 0.35 + 0.25 * tanh((dist1 - 0.2)/(2*Ca));

    if(y<0.99999)
      value[3] = temp_u;
    else 
      value[3] = temp_h;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = one bubble, u = 0.0, v = 0.0 ");
  SYS_T::commPrint("theta = 0.8 bottom 0.9 top. \n");
}


void PDNSolution_TNSK_2D::Init_boiling( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();
  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.990;

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  //srand(time(NULL));
  //double discrete_rand;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    value[0] = 0.37 - 0.31 * tanh((y - 0.5)/(2*Ca));

    /* 
       if(y<1.45)
       { 
       discrete_rand = rand() % 101;
       value[0] += (double) discrete_rand * 1.0e-5 - 0.0005;
       }
       */

    if(y < 0.01  && abs(x-0.5)<0.1)
      value[3] = temp2;
    else
      value[3] = temp1;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = bottom liquid top vapor ");
  SYS_T::commPrint("zero velocity, theta = 0.8 bottom 0.9 top. \n");
}


void PDNSolution_TNSK_2D::Init_RB_convection( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  //const double Ca = locassem_ptr->get_model_para_1();
  const double temp1 = -1.0 / 0.850;
  const double temp2 = -1.0 / 0.950;

  //const double cst1 = 0.5 * (temp1 + temp2);
  //const double cst2 = 0.5 * (temp1 - temp2);

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  double discrete_rand, random_num, random_num_2;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    value[0] = 0.6024;

    discrete_rand = rand() % 1001;

    value[0] += (double) discrete_rand * 1.0e-6 - 5.0e-4;

    discrete_rand = rand() % 1001;

    random_num = (double) discrete_rand * 1.0e-5 - 5.0e-3;

    discrete_rand = rand() % 1001;

    random_num_2 = (double) discrete_rand * 1.0e-5 - 5.0e-3; 

    value[3] = 2.0 * y * (temp1+random_num_2) + (1.0 - 2.0*y)*(temp2 + random_num);

    //value[3] = cst1 + 0.5 * random_num_2
    //  + (cst2 - 0.5 * random_num_2) * tanh( (y-(0.01+random_num))/(2.0*Ca) ); 

    //if(y>0.001)
    //  value[3] = temp1 + random_num;
    //else
    //  value[3] = temp2 + random_num_2;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = single phase ");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top. \n");
}


void PDNSolution_TNSK_2D::Init_RB_convection_2( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();
  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  double discrete_rand, random_num;

  double inter_y;
  //const double pi = 3.1415926536;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    // density
    //inter_y = 0.9 + 0.001 * sin(10*pi*x);
    inter_y = 0.9;

    value[0] = 0.366 - 0.2971 * tanh((y - inter_y)/(2*Ca));

    // temperature
    discrete_rand = rand() % 1001;

    random_num = ( (double) discrete_rand * 1.0e-5 - 5.0e-3 ) * 10.0;

    if( y < 1.0e-5 )
      value[3] = temp2 + random_num;
    else
      value[3] = temp1;


    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = single phase ");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top. \n");
}


// 9. spindoal decomposition under heated bc with random perturbations. 
void PDNSolution_TNSK_2D::Init_thermal_spinodal_decomp( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr,
    const IALocal_BC * const &lbc )
{
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.875;
  const double boundary_temp = -1.0 / 0.90;

  // loop over to insert density, velocity.
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    value[0] = 0.49;

    value[3] = uniform_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  // undate the temperature value on dirichlet nodes  
  const int num = lbc->get_Num_LD(3); // get the dirichlet nodes for temperature
  int * bc_index = new int [num];
  double * bc_value = new double [num];

  srand(time(NULL));
  double discrete_rand, random_num;

  for(int ii=0; ii<num; ++ii)
  {
    bc_index[ii] = lbc->get_LDN(3, ii) * 5 + 3;

    discrete_rand = rand() % 1001;
    random_num = ( (double) discrete_rand * 1.0e-4 - 5.0e-2 ) * 2.0;
    bc_value[ii] = boundary_temp + random_num;
  }

  VecSetValues(solution, num, bc_index, bc_value, INSERT_VALUES);
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  delete [] bc_index; delete [] bc_value; 

  SYS_T::commPrint("===> Initial solution: rho = random, u = 0.0, v = 0.0, ");
  PetscPrintf(PETSC_COMM_WORLD, 
      "theta = %e inside and %e on boundary \n", -1.0/uniform_temp, -1.0/boundary_temp);
  PetscPrintf(PETSC_COMM_WORLD, " Boundary temperature at Dir nodes has random perturbations. \n");
}



void PDNSolution_TNSK_2D::Init_twophase_boiling( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();
  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;
  const double temp3 = -1.0 / 0.850;

  const double tan_a = 0.5 * (temp1 + temp3);
  const double tan_b = 0.5 * (temp3 - temp1);

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  double discrete_rand, random_num;

  double inter_y;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    // interface location 
    inter_y = 0.70;

    value[0] = 0.33565 - 0.26675 * tanh((y - inter_y)/(2*Ca));

    // temperature
    discrete_rand = rand() % 1001;

    random_num = ( (double) discrete_rand * 1.0e-5 - 5.0e-3 ) * 10.0;

    if( y < 1.0e-5 )
      value[3] = temp2 + random_num;
    else if(y>0.99999)
      value[3] = temp1 + random_num * 0.1;
    else
      value[3] = tan_a - tan_b * tanh((y-inter_y)/(2*Ca));


    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = two phase interface y = 0.7 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top both with randomness. \n");
  SYS_T::commPrint(" temperature in liquid: 0.85, temperature in vapor: 0.775");
}


// 11. Boiling in [0,1] x [0, 0.5] 
void PDNSolution_TNSK_2D::Init_twophase_boiling_rec( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();
  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;
  const double temp3 = -1.0 / 0.775;

  const double tan_a = 0.5 * (temp1 + temp3);
  const double tan_b = 0.5 * (temp3 - temp1);

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  double discrete_rand, random_num;

  const double inter_y = 0.35;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    value[0] = 0.366 - 0.2971 * tanh((y - inter_y)/(2*Ca));

    // temperature
    discrete_rand = rand() % 1001;

    random_num = ( (double) discrete_rand * 1.0e-5 - 5.0e-3 ) * 10.0;

    if( y < 1.0e-5 )
      value[3] = temp2 + random_num;
    else if(y > 0.5 - 1.0e-4)
      value[3] = temp1 + random_num * 0.1;
    else
      value[3] = tan_a - tan_b * tanh((y-inter_y)/(2*Ca));


    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = two phase interface y = 0.4 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top both with randomness. \n");
  SYS_T::commPrint(" temperature in liquid: 0.775, temperature in vapor: 0.775");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}

void PDNSolution_TNSK_2D::Init_multibub_cBothTop( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();

  double dist1, x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double inv_2Ca = 0.5 / Ca;

  const double unif_temp = -1.0 / 0.85;
  const double top_temp = -1.0 / 0.87;

  // ----- Define the position and radius of the bubbles
  const int num_bub = 8;

  std::vector<double> xx, yy, rr;

  xx.push_back(0.1); yy.push_back(0.1); rr.push_back(0.05);
  xx.push_back(0.15); yy.push_back(0.7); rr.push_back(0.05);
  xx.push_back(0.2); yy.push_back(0.5); rr.push_back(0.05);
  xx.push_back(0.25); yy.push_back(0.8); rr.push_back(0.05);
  xx.push_back(0.3); yy.push_back(0.2); rr.push_back(0.05);
  xx.push_back(0.35); yy.push_back(0.6); rr.push_back(0.05);
  xx.push_back(0.4); yy.push_back(0.4); rr.push_back(0.05);
  xx.push_back(0.42); yy.push_back(0.8); rr.push_back(0.05);
  // ----- End of the Definition

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    value[0] = 0.0;
    for(int jj=0; jj<num_bub; ++jj)
    {
      dist1 = sqrt((x-xx[jj])*(x-xx[jj]) + (y-yy[jj])*(y-yy[jj]));
      value[0] += 0.25 * tanh((dist1 - rr[jj]) * inv_2Ca);
    }
    value[0] += 0.6 - 0.25 * num_bub;

    if(y<0.99999)
      value[3] = unif_temp;
    else
      value[3] = top_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = Multi bubble, u = 0.0, v = 0.0 ");
  SYS_T::commPrint("theta = 0.85, top 0.87. \n");
}


void PDNSolution_TNSK_2D::Init_1bub_cBothTop_2( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();

  double dist1, x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();
  
  const double temp_u = -1.0 / 0.85;
  const double temp_h = -1.0 / 0.87;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));

    value[0] = 0.35 + 0.25 * tanh((dist1 - 0.2)/(2*Ca));

    if(y<0.99999)
      value[3] = temp_u;
    else 
      value[3] = temp_h;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = one bubble, u = 0.0, v = 0.0 ");
  SYS_T::commPrint("theta = 0.85, top heated to 0.87. \n");
}


void PDNSolution_TNSK_2D::Init_twophase_boiling_rec_new( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();
  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;
  const double temp3 = -1.0 / 0.85;

  const double tan_a = 0.5 * (temp1 + temp3);
  const double tan_b = 0.5 * (temp3 - temp1);

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  double discrete_rand, random_num;

  const double inter_y = 0.35;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    value[0] = 0.33565 - 0.26675 * tanh((y - inter_y)/(2*Ca));

    // temperature
    discrete_rand = rand() % 1001;

    random_num = ( (double) discrete_rand * 1.0e-5 - 5.0e-3 ) * 10.0;

    if( y < 1.0e-5 )
      value[3] = temp2 + random_num;
    else if(y > 0.5 - 1.0e-4)
      value[3] = temp1 + random_num * 0.1;
    else
      value[3] = tan_a - tan_b * tanh((y-inter_y)/(2*Ca));


    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = two phase interface y = 0.35 \n");
  SYS_T::commPrint("zero velocity, theta: 0.95 bottom, 0.775 top both with randomness. \n");
  SYS_T::commPrint(" temperature in liquid: 0.85, temperature in vapor: 0.775");
  SYS_T::commPrint(" density in liquid: 0.6024, in vapor: 0.0689");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}


void PDNSolution_TNSK_2D::Init_twophase_boiling_sin( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double pi = 3.1415926535898;
  const double Ca = locassem_ptr->get_model_para_1();
  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;
  const double temp3 = -1.0 / 0.775;

  const double tan_a = 0.5 * (temp1 + temp3);
  const double tan_b = 0.5 * (temp3 - temp1);

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  //double discrete_rand;
  double random_num;

  const double inter_y = 0.35;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    value[0] = 0.366 - 0.2971 * tanh((y - inter_y)/(2*Ca));

    // temperature
    //discrete_rand = rand() % 1000;
    //random_num = ( (double) discrete_rand * 1.0e-5 - 5.0e-3 ) * 10.0;
    
    random_num = 5.0e-2 * sin(10*pi*x);

    if( y < 1.0e-5 )
      value[3] = temp2 + random_num;
    else if(y > 0.5 - 1.0e-4)
      value[3] = temp1 + random_num * 0.1;
    else
      value[3] = tan_a - tan_b * tanh((y-inter_y)/(2*Ca));


    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = two phase interface y = 0.35 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top both with randomness. \n");
  SYS_T::commPrint(" temperature in liquid: 0.775, temperature in vapor: 0.775");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}


// 16. spindoal decomposition under heated bc with random perturbations. 
void PDNSolution_TNSK_2D::Init_thermal_spinodal_decomp_2( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_NodalBC * const &lbc )
{
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.875;
  const double boundary_temp = -1.0 / 0.90;

  // loop over to insert density, velocity.
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    value[0] = 0.49;

    value[3] = uniform_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  // undate the temperature value on dirichlet nodes  
  const int num = lbc->get_Num_LD(3); // get the dirichlet nodes for temperature
  int * bc_index = new int [num];
  double * bc_value = new double [num];

  srand(time(NULL));
  double discrete_rand, random_num;

  for(int ii=0; ii<num; ++ii)
  {
    bc_index[ii] = lbc->get_LDN(3, ii) * 5 + 3;

    discrete_rand = rand() % 1001;
    random_num = ( (double) discrete_rand * 1.0e-4 - 5.0e-2 ) * 2.0;
    bc_value[ii] = boundary_temp + random_num;
  }

  VecSetValues(solution, num, bc_index, bc_value, INSERT_VALUES);
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  delete [] bc_index; delete [] bc_value; 

  SYS_T::commPrint("===> Initial solution: rho = random, u = 0.0, v = 0.0, ");
  PetscPrintf(PETSC_COMM_WORLD, 
      "theta = %e inside and %e on boundary \n", -1.0/uniform_temp, -1.0/boundary_temp);
  PetscPrintf(PETSC_COMM_WORLD, " Boundary temperature at Dir nodes has random perturbations. \n");
}


void PDNSolution_TNSK_2D::Init_newRB( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double pi = 3.1415926535898;
  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  double random_num;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    value[0] = 0.5631;

    random_num = 5.0e-2 * sin(10*pi*x);

    if( y < 1.0e-5 )
      value[3] = temp2 + random_num;
    else if(y > 0.5 - 1.0e-4)
      value[3] = temp1 + random_num * 0.1;
    else
      value[3] = temp1;


    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = 0.5631 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top both with randomness. \n");
  SYS_T::commPrint(" temperature in bulk: 0.775");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}


// 18. Two-phase boiling with white noise perturbation on temp. bc.
void PDNSolution_TNSK_2D::Init_twophase_boiling_wn( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();
  const double temp1 = -1.0 / 0.775;
  //const double temp2 = -1.0 / 0.950;
  const double temp3 = -1.0 / 0.775;

  const double tan_a = 0.5 * (temp1 + temp3);
  const double tan_b = 0.5 * (temp3 - temp1);

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  double discrete_rand, random_num;

  const double inter_y = 0.35;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    value[0] = 0.366 - 0.2971 * tanh((y - inter_y)/(2*Ca));

    // temperature
    discrete_rand = rand() % 1000;
    random_num = (double) discrete_rand * 1.0e-5 - 5.0e-3;
    
    if( y < 1.0e-5 )
      value[3] = (-1.0) / ( 0.950 + random_num );
    else if(y > 0.5 - 1.0e-5)
      value[3] = (-1.0) / ( 0.775 + random_num );
    else
      value[3] = tan_a - tan_b * tanh((y-inter_y)/(2*Ca));


    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = two phase interface y = 0.35 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top both with white noise. \n");
  SYS_T::commPrint(" temperature in liquid: 0.950, temperature in vapor: 0.775");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}



// 19. Two-phase boiling with white noise perturbation on linear bulk temp..
void PDNSolution_TNSK_2D::Init_twophase_boiling_linearRand_temp( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double Ca = locassem_ptr->get_model_para_1();
  
  //const double temp1 = -1.0 / 0.775;
  //const double temp2 = -1.0 / 0.950;
  //const double temp3 = -1.0 / 0.775;

  //const double tan_a = 0.5 * (temp1 + temp3);
  //const double tan_b = 0.5 * (temp3 - temp1);

  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  double discrete_rand, random_num;

  const double inter_y = 0.35;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    // two phase density
    value[0] = 0.366 - 0.2971 * tanh((y - inter_y)/(2*Ca));

    // temperature
    discrete_rand = rand() % 1000;
    random_num = (double) discrete_rand * 1.0e-5 - 5.0e-3;
    
    if( y < 1.0e-7 )
      value[3] = (-1.0) / 0.950;
    else if(y > 0.5 - 1.0e-7)
      value[3] = (-1.0) / 0.775;
    else
      value[3] = (-1.0) / ( 0.950 - 0.35 * y + random_num );


    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = two phase interface y = 0.35 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top. \n");
  SYS_T::commPrint(" temperature in bulk, linear with random perturbations. \n");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}



// 20. Two-phase boiling with buld liquid. 
void PDNSolution_TNSK_2D::Init_twophase_boiling_homogeneous( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));
  double discrete_rand, random_num;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    //value[0] = 0.6631;
    value[0] = 0.5000;

    // temperature
    discrete_rand = rand() % 1000;
    random_num = (double) discrete_rand * 1.0e-5 - 5.0e-3;
    
    value[3] = (-1.0) / (0.950 - 0.35 * y + random_num); 

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = 0.5 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top both with white noise. \n");
  SYS_T::commPrint(" temperature in liquid: 0.950, temperature in vapor: 0.775");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}



// 21. Two-phase boiling with buld liquid, no pert of temperature on bc. 
void PDNSolution_TNSK_2D::Init_twophase_boiling_homogeneous_nopert( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    //value[0] = 0.6631;
    value[0] = 0.5000;

    // temperature
    value[3] = (-1.0) / (0.950 - 0.35 * y); 

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = 0.5 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top with NO white noise. \n");
  SYS_T::commPrint(" temperature in liquid: 0.950, temperature in vapor: 0.775");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}



// 22. Two-phase boiling with buld liquid, no pert of temperature on bc, 
//     bottom 0.95, top 0.78 
void PDNSolution_TNSK_2D::Init_twophase_boiling_homogeneous_nopert_0d78( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  double x, y;
  int location[5];
  double value[5];

  const double rho_ave = 0.8;

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    // density
    value[0] = rho_ave;

    // temperature
    value[3] = (-1.0) / (0.950 - 0.34 * y); 
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  PetscPrintf(PETSC_COMM_WORLD, "===> Initial solution: rho = %e \n", rho_ave);
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top with NO white noise. \n");
  SYS_T::commPrint(" temperature in liquid: 0.950, temperature in vapor: 0.78");
  SYS_T::commPrint(" This initial solution is prepared for [0,1] by [0, 0.5] domain. \n");
}



// 23. Two-phase boiling with buld liquid, no pert of temperature on bc, 
//     bottom 0.95, top 0.78, geometry is [0,1]^2 
void PDNSolution_TNSK_2D::Init_twophase_boiling_homogeneous_nopert_0d78_square( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  double x, y;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; 
  value[3] = 0.0; value[4] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  srand(time(NULL));

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);

    // density
    value[0] = 0.5;

    // temperature
    value[3] = (-1.0) / (0.950 - 0.17 * y); 
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = 0.5 \n");
  SYS_T::commPrint("zero velocity, theta: hot bottom, cold top with NO white noise. \n");
  SYS_T::commPrint(" temperature in liquid: 0.950, temperature in vapor: 0.78");
  SYS_T::commPrint(" This initial solution is prepared for [0,1]x[0,1] domain. \n");
}


// EOF
