#include "PDNSolution_TNSK_3D.hpp"

PDNSolution_TNSK_3D::PDNSolution_TNSK_3D( const APart_Node * const &pNode,
    const IALocal_BC * const &bc,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution(pNode)
{
  switch (type)
  {
    case 0:
      Init_case0( pNode, locassem_ptr, fNode_ptr );
      break;
    case 1:
      Init_case1( pNode, locassem_ptr, fNode_ptr );
      break;
    case 2:
      Init_case2( pNode, locassem_ptr, fNode_ptr );
      break;
    case 3:
      Init_case3( pNode, locassem_ptr, fNode_ptr );
      break;
    case 4:
      Init_case4( pNode, locassem_ptr, fNode_ptr );
      break;
    case 5:
      Init_case5( pNode, locassem_ptr, fNode_ptr );
      break;
    case 6:
      Init_case6( pNode, locassem_ptr, fNode_ptr );
      break;
    default:
      SYS_T::commPrint("ERROR: PDNSolution_TNSK_3D: No such type of initial solution. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}


PDNSolution_TNSK_3D::PDNSolution_TNSK_3D(
    const class APart_Node * const &pNode,
    const class ALocal_NodalBC * const &lnbc,
    const class IPLocAssem * const &locassem_ptr,
    const class FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution(pNode)
{
  switch (type)
  {
    case 0:
      Init_case0( pNode, locassem_ptr, fNode_ptr );
      break;
    case 1:
      Init_case1( pNode, locassem_ptr, fNode_ptr );
      break;
    case 2:
      Init_case2( pNode, locassem_ptr, fNode_ptr );
      break;
    case 3:
      Init_case3( pNode, locassem_ptr, fNode_ptr );
      break;
    case 4:
      Init_case4( pNode, locassem_ptr, fNode_ptr );
      break;
    case 5:
      Init_case5( pNode, locassem_ptr, fNode_ptr );
      break;
    case 6:
      Init_case6( pNode, locassem_ptr, fNode_ptr );
      break;
    case 7:
      Init_case7( pNode, locassem_ptr, fNode_ptr );
      break;
    default:
      SYS_T::commPrint("ERROR: PDNSolution_TNSK_3D: No such type of initial solution. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}



PDNSolution_TNSK_3D::~PDNSolution_TNSK_3D()
{}


void PDNSolution_TNSK_3D::Init_case0( const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[6];
  double value[6];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0;
  value[3] = 0.0; value[4] = 0.0; value[5] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.85;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    value[0] = 0.6;
    value[4] = uniform_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 6;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: rho = 0.6, u = 0.0, v = 0.0, theta = 0.85. \n");
}


void PDNSolution_TNSK_3D::Init_case1( const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[6];
  double value[6];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0;
  value[3] = 0.0; value[4] = 0.0; value[5] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double uniform_temp = -1.0 / 0.85;
  const double Ca = locassem_ptr->get_model_para_1();
  const double inv_2Ca = 0.5 / Ca;
  double dist1, x, y, z;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));

    value[0] = 0.35 + 0.25 * (tanh((dist1 - 0.25)*inv_2Ca));
    value[4] = uniform_temp;

    location[0] = pNode_ptr->get_node_loc(ii) * 6;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: single bubble in center, \n");
  SYS_T::commPrint("     u = 0.0, v = 0.0, theta = 0.85. \n");
}


void PDNSolution_TNSK_3D::Init_case2( const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[6];
  double value[6];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0;
  value[3] = 0.0; value[4] = 0.0; value[5] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double Ca = locassem_ptr->get_model_para_1();
  const double inv_2Ca = 0.5 / Ca;
  double dist1, dist2, dist3, x, y, z;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    dist1 = sqrt((x-0.2)*(x-0.2) + (y-0.2)*(y-0.2) + (z-0.2)*(z-0.2));

    dist2 = sqrt((x-0.3)*(x-0.3) + (y-0.3)*(y-0.3) + (z-0.5)*(z-0.5));

    dist3 = sqrt((x-0.2)*(x-0.2) + (y-0.3)*(y-0.3) + (z-0.8)*(z-0.8));

    value[0] = -0.15 + 0.25 * (tanh((dist1 - 0.1)*inv_2Ca)) 
      + 0.25 * tanh((dist2 - 0.1) * inv_2Ca) + 0.25 * tanh((dist3-0.1)*inv_2Ca);

    value[4] = -1.0 / (0.8 + 0.1 * z);

    location[0] = pNode_ptr->get_node_loc(ii) * 6;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: three bubbles, \n");
  SYS_T::commPrint("     u = 0.0, v = 0.0, bottom 0.8 top 0.9 uniform graddient. \n");
  SYS_T::commPrint(" NOTE: Geometry should be [0,0.5]x[0,0.5]x[0,1]. \n");
}



void PDNSolution_TNSK_3D::Init_case3( const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[6];
  double value[6];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0;
  value[3] = 0.0; value[4] = 0.0; value[5] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double Ca = locassem_ptr->get_model_para_1();
  const double inv_2Ca = 0.5 / Ca;
  double dist1, x, y, z;

  const double temp_u = -1.0 / 0.85;

  const double temp_t = -1.0 / 0.87;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    dist1 = sqrt((x-0.25)*(x-0.25) + (y-0.25)*(y-0.25) + (z-0.3)*(z-0.3));

    value[0] = 0.35 + 0.25 * (tanh((dist1 - 0.2)*inv_2Ca)); 

    if(z<0.99999)
      value[4] = temp_u;
    else
      value[4] = temp_t;

    location[0] = pNode_ptr->get_node_loc(ii) * 6;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  SYS_T::commPrint("===> Initial solution: single bubbles, \n");
  SYS_T::commPrint("     u = 0.0, v = 0.0, bottom 0.85 top 0.87. \n");
  SYS_T::commPrint(" NOTE: Geometry should be [0,0.5]x[0,0.5]x[0,1]. \n");
}


void PDNSolution_TNSK_3D::Init_case4( const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[6];
  double value[6];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0;
  value[3] = 0.0; value[4] = 0.0; value[5] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double Ca = locassem_ptr->get_model_para_1();
  const double inv_2Ca = 0.5 / Ca;
  double z;

  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;
  const double temp3 = -1.0 / 0.850;

  const double tan_a = 0.5 * (temp1 + temp3);
  const double tan_b = 0.5 * (temp3 - temp1);

  srand(time(NULL));
  double discrete_rand, random_num;

  const double inter_z = 0.15;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    //x = fNode_ptr->get_ctrlPts_x(ii);
    //y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    // density
    value[0] = 0.33565 - 0.26675 * tanh((z - inter_z)*inv_2Ca);

    // temperature 
    discrete_rand = rand() % 1001;
    random_num = ( (double) discrete_rand * 1.0e-5 - 5.0e-3 ) * 10.0;

    if( z < 1.0e-5 )
      value[4] = temp2 + random_num;
    else if( z > 0.25 - 1.0e-4 )
      value[4] = temp1 + random_num * 0.1;
    else
      value[4] = tan_a - tan_b * tanh((z-inter_z)*inv_2Ca);

    // set location 
    location[0] = pNode_ptr->get_node_loc(ii) * 6;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: rho = two phase interface z = 0.15 \n");
  SYS_T::commPrint("zero velocity, theta: 0.95 bottom, 0.775 top both with randomness. \n");
  SYS_T::commPrint(" temperature in liquid: 0.85, temperature in vapor: 0.775");
  SYS_T::commPrint(" density in liquid: 0.6024, in vapor: 0.0689");
  SYS_T::commPrint(" NOTE: Geometry should be [0,1]x[0,0.5]x[0,0.25]. \n");
}


void PDNSolution_TNSK_3D::Init_case5( const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[6];
  double value[6];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0;
  value[3] = 0.0; value[4] = 0.0; value[5] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  const double Ca = locassem_ptr->get_model_para_1();
  const double inv_2Ca = 0.5 / Ca;
  double z;

  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;
  const double temp3 = -1.0 / 0.850;

  const double tan_a = 0.5 * (temp1 + temp3);
  const double tan_b = 0.5 * (temp3 - temp1);

  srand(time(NULL));
  double discrete_rand, random_num;

  const double inter_z = 0.35;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    //x = fNode_ptr->get_ctrlPts_x(ii);
    //y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    // density
    value[0] = 0.33565 - 0.26675 * tanh((z - inter_z)*inv_2Ca);

    // temperature 
    discrete_rand = rand() % 1001;
    random_num = ( (double) discrete_rand * 1.0e-5 - 5.0e-3 ) * 10.0;

    if( z < 1.0e-5 )
      value[4] = temp2 + random_num;
    else if( z > 0.25 - 1.0e-4 )
      value[4] = temp1 + random_num * 0.1;
    else
      value[4] = tan_a - tan_b * tanh((z-inter_z)*inv_2Ca);

    // set location 
    location[0] = pNode_ptr->get_node_loc(ii) * 6;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: rho = two phase interface z = 0.35 \n");
  SYS_T::commPrint("zero velocity, theta: 0.95 bottom, 0.775 top both with randomness. \n");
  SYS_T::commPrint(" temperature in liquid: 0.85, temperature in vapor: 0.775");
  SYS_T::commPrint(" density in liquid: 0.6024, in vapor: 0.0689");
  SYS_T::commPrint(" NOTE: Geometry should be [0,1]x[0,0.4]x[0,0.4]. \n");
}


void PDNSolution_TNSK_3D::Init_case6( const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[6];
  double value[6];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0;
  value[3] = 0.0; value[4] = 0.0; value[5] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  double z;

  const double temp1 = -1.0 / 0.775;
  const double temp2 = -1.0 / 0.950;
  const double temp3 = -1.0 / 0.850;

  srand(time(NULL));
  double discrete_rand;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    z = fNode_ptr->get_ctrlPts_z(ii);

    discrete_rand = rand() % 1001;

    // density
    value[0] = 0.5 + (double) discrete_rand * 1.0e-5 - 5e-3;

    // temperature 
    if( z < 1.0e-5 )
      value[4] = temp2;
    else if( z > 1.0 - 1.0e-4 )
      value[4] = temp1;
    else
      value[4] = temp3;

    // set location 
    location[0] = pNode_ptr->get_node_loc(ii) * 6;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: rho = 0.5 + randomness \n");
  SYS_T::commPrint("zero velocity, theta: 0.95 bottom, 0.775 top, uniformly. \n");
  SYS_T::commPrint(" theta in middle 0.85");
  SYS_T::commPrint(" NOTE: Geometry should be [0,1]x[0,1]x[0,1]. \n");
}



void PDNSolution_TNSK_3D::Init_case7( const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  const double rho_ave = 0.8;

  int location[6];
  double value[6];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0;
  value[3] = 0.0; value[4] = 0.0; value[5] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();
  double z;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    z = fNode_ptr->get_ctrlPts_z(ii);

    // density
    value[0] = rho_ave;

    // temperature 
    value[4] = (-1.0) / (0.95 - 0.34 * z);

    // set location 
    location[0] = pNode_ptr->get_node_loc(ii) * 6;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;

    VecSetValues(solution, 6, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  PetscPrintf(PETSC_COMM_WORLD, "===> Initial solution: rho = %e \n", rho_ave);
  SYS_T::commPrint("zero velocity, theta: 0.95 bottom, 0.78 top, linear profile in z-dir. \n");
  SYS_T::commPrint(" NOTE: Geometry should be [0,1]x[0,1]x[0,0.5]. \n");
}

// EOF
