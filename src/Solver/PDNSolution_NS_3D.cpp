#include "PDNSolution_NS_3D.hpp"

PDNSolution_NS_3D::PDNSolution_NS_3D(const APart_Node * const &pNode,
    const IALocal_BC * const &bc,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution(pNode)
{
  switch(type)
  {
    case 0:
      Init_test( pNode, locassem_ptr, fNode_ptr );
      break;
    case 1:
      Init_test_2( pNode, locassem_ptr, fNode_ptr );
      break;
    case 2:
      Init_3D_Driven_Cavity( pNode, locassem_ptr, fNode_ptr );
      break;
    case 3:
      Init_3D_Poiseuille_Pipe( pNode, locassem_ptr, fNode_ptr );
      break;
    case 4:
      Init_3D_Pressure_Pipe( pNode, locassem_ptr, fNode_ptr );
      break;
    case 5:
      Init_random( pNode, locassem_ptr, fNode_ptr );
      break;
    case 6:
      Init_3D_Pressure_Coronary_3patch( pNode, locassem_ptr, fNode_ptr );
      break;
    case 7:
      Init_3D_Dir_inflow_Coronary_3patch( pNode, locassem_ptr, fNode_ptr );
      break;
    case 8:
      Init_3D_Pressure_Coronary_107patch( pNode, locassem_ptr, fNode_ptr );
      break;
    default:
      SYS_T::commPrint("ERROR: PDNSolution_NS_3D: No such type of initial condition. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}

PDNSolution_NS_3D::~PDNSolution_NS_3D()
{}

void PDNSolution_NS_3D::Init_test(const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

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

  SYS_T::commPrint("===> Initial solution: Test Case: Full zero vector. \n");
}

void PDNSolution_NS_3D::Init_test_2(const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];
  //double x, y;

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    //x = fNode_ptr->get_ctrlPts_x(ii);
    //y = fNode_ptr->get_ctrlPts_y(ii);

    //value[3] = sin(3.1415926 * x) * sin(3.1415926 * y);

    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: Convergence Study. \n");
}


void PDNSolution_NS_3D::Init_3D_Driven_Cavity( 
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];
  double x, y, z, pow_x, pow_y;

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    if( z > 1.0 - 1.0e-4 )
    {
      pow_x = pow((2.0*x-1), 18);
      pow_y = pow((2.0*y-1), 18);
      value[0] = (1.0 - pow_x) * (1.0 - pow_y);
    }

    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: u = (1-(2x-1)^18)(1-(2y-1)^18) \n");
  SYS_T::commPrint("                       v = 0.0  \n");
  SYS_T::commPrint("                       w = 0.0  \n");
  SYS_T::commPrint("                       p = 0.0  \n");
  SYS_T::commPrint(" Note: Geometry is [0,1]^3.     \n");
}

void PDNSolution_NS_3D::Init_3D_Poiseuille_Pipe( 
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];
  double x, y, z, x2y2;

  const double radius2 = 0.05 * 0.05;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    x2y2 = x*x + y*y;
    
    if( z > 1.0 - 1.0e-6 )
      value[2] = 5.0 * (x2y2 - radius2);

    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);

  }
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: w = -1.0 on top boundary \n");
  SYS_T::commPrint("                       otherwise, vel = 0.0.     \n");
  SYS_T::commPrint("                       p = 0.0.     \n");
  SYS_T::commPrint("Note: geometry_lumen_original.txt \n"); 
}  


void PDNSolution_NS_3D::Init_3D_Pressure_Pipe( 
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];

  int nlocalnode = pNode_ptr->get_nlocalnode();

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    double z = fNode_ptr->get_ctrlPts_z(ii);
  
    value[3] = 1.0 * z;

    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);

  }
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: vel = 0.0 \n");
  SYS_T::commPrint("              p = 1.0 on top bc. \n");
  SYS_T::commPrint("              p = 0.0 on bot bc. \n");
  SYS_T::commPrint("Note: geometry_lumen_original.txt \n"); 
}  



void PDNSolution_NS_3D::Init_random(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];

  srand(time(NULL));

  int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    value[0] = rand() % 100001;
    value[1] = rand() % 100001;
    value[2] = rand() % 100001;
    value[3] = rand() % 100001;

    value[0] = value[0] * 1.0e-5;
    value[1] = value[1] * 1.0e-5;
    value[2] = value[2] * 1.0e-5;
    value[3] = value[3] * 1.0e-5;

    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: Test Case: Random vector. \n");
}




void PDNSolution_NS_3D::Init_3D_Pressure_Coronary_3patch(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];
  double x, y, z;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

  // n0 dot (x - x0) = 0 gives the bottom plane of patch 0
  const double n0[3] = {0.746313853044165, -0.537450240675482, 0.392635800140842};
  const double x0[3] = {39.1826992412, -14.0891355102, -12.5904822952};

  const double n1[3] = { 0.029420819084185, -0.855635931543001, 0.516741297030686 };
  const double x1[3] = {  52.622072791500003, -46.918228118899997, 6.468909591290000 };
 
  const double n2[3] = { 0.579364198916169, -0.115531629292201, 0.806839245232478 };
  const double x2[3] = { 62.654280889100001, -30.991856628400001, 8.888908530249999 }; 

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    double dist0 = std::abs(n0[0]*(x-x0[0]) + n0[1]*(y-x0[1]) + n0[2]*(z-x0[2]));
    double dist1 = std::abs(n1[0]*(x-x1[0]) + n1[1]*(y-x1[1]) + n1[2]*(z-x1[2]));
    double dist2 = std::abs(n2[0]*(x-x2[0]) + n2[1]*(y-x2[1]) + n2[2]*(z-x2[2]));

    value[3] = 0.0;

    if(dist0 < 1.0e-1)
      value[3] = 2.0;

    if(dist1 < 1.0e-1)
      value[3] = 1.0;

    if(dist2 < 1.0e-1)
      value[3] = 0.5;

    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);

  }
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: vel = 0.0 \n");
  SYS_T::commPrint("              p = 2.0 on patch 0 bot. \n");
  SYS_T::commPrint("              p = 1.0 on patch 1 top. \n");
  SYS_T::commPrint("              p = 1.0 on patch 2 bot. \n");
  SYS_T::commPrint("Note: Geometry is coronary 3-patch case. \n"); 
}


void PDNSolution_NS_3D::Init_3D_Dir_inflow_Coronary_3patch(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];
  double x, y, z;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

  // n0 dot (x - x0) = 0 gives the bottom plane of patch 0
  const double n0[3] = {0.746313853044165, -0.537450240675482, 0.392635800140842};
  const double x0[3] = {39.1826992412, -14.0891355102, -12.5904822952};

  const double n1[3] = { 0.029420819084185, -0.855635931543001, 0.516741297030686 };
  const double x1[3] = {  52.622072791500003, -46.918228118899997, 6.468909591290000 };
 
  const double n2[3] = { 0.579364198916169, -0.115531629292201, 0.806839245232478 };
  const double x2[3] = { 62.654280889100001, -30.991856628400001, 8.888908530249999 }; 

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    double dist0 = std::abs(n0[0]*(x-x0[0]) + n0[1]*(y-x0[1]) + n0[2]*(z-x0[2]));
    double dist1 = std::abs(n1[0]*(x-x1[0]) + n1[1]*(y-x1[1]) + n1[2]*(z-x1[2]));
    double dist2 = std::abs(n2[0]*(x-x2[0]) + n2[1]*(y-x2[1]) + n2[2]*(z-x2[2]));

    double distc = sqrt( (x-x0[0])*(x-x0[0]) + (y-x0[1])*(y-x0[1]) + (z-x0[2])*(z-x0[2]) );

    if(dist0 < 1.0e-2 && distc < 1.0 )
    {
      //value[0] = 1.0e2 * n0[0] * (1.0 - distc) * (1.0 - distc) * (0.1 - dist0);
      //value[1] = 1.0e2 * n0[1] * (1.0 - distc) * (1.0 - distc) * (0.1 - dist0);
      //value[2] = 1.0e2 * n0[2] * (1.0 - distc) * (1.0 - distc) * (0.1 - dist0);
      value[0] = 1.0e0 * n0[0];
      value[1] = 1.0e0 * n0[1];
      value[2] = 1.0e0 * n0[2];
    }
    else
    {
      value[0] = 0.0;
      value[1] = 0.0;
      value[2] = 0.0;
    }

    if(dist1 < 1.0e-1)
      value[3] = 0.0;

    if(dist2 < 1.0e-1)
      value[3] = 0.0;
    
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }
  
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: vel = -1.0 normal on patch 0 bottom \n");
  SYS_T::commPrint("     Pressure on patch 1 top    : 0.0 \n");
  SYS_T::commPrint("     Pressure on patch 2 bottom : 0.0 \n");
  SYS_T::commPrint("Note: Geometry is coronary 3-patch case. \n"); 
}


void PDNSolution_NS_3D::Init_3D_Pressure_Coronary_107patch(
    const APart_Node * const &pNode_ptr,
    const IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  int location[4];
  double value[4];
  double x, y, z;

  int nlocalnode = pNode_ptr->get_nlocalnode();

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

  // n0 dot (x - x0) = 0 gives the bottom plane of patch 0
  const double n0[3] = {0.407477, 0.4842211, 0.77426898};
  const double x0[3] = {-15.068475168, -1.189208589, 62.5963674870};

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    double dist0 = std::abs(n0[0]*(x-x0[0]) + n0[1]*(y-x0[1]) + n0[2]*(z-x0[2]));

    value[3] = 0.0;

    if(dist0 < 1.0e-1)
      value[3] = 7.99932e2;

    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);

  }
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);

  GhostUpdate();

  SYS_T::commPrint("===> Initial solution: vel = 0.0 \n");
  SYS_T::commPrint("              p = 799 on patch 0 bot. \n");
  SYS_T::commPrint("Note: Geometry is coronary 107-patch case. \n"); 
}


// EOF
