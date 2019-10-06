#include "PDNSolution_SINSK.hpp"

PDNSolution_SINSK::PDNSolution_SINSK( const class APart_Node * const &pNode,
    const class IALocal_BC * const &locbc,
    const class IPLocAssem * const &locassem_ptr,
    const class FEANode * const &fNode_ptr,
    const int &type ) : PDNSolution(pNode)
{
  switch (type)
  {
    case 0:
      Init_Zero();
      SYS_T::commPrint("===> Initial solution: Zero vector. \n");
      break;
    case 1:
      Init_OneStaticBubble(pNode, locassem_ptr, fNode_ptr);
      SYS_T::commPrint("===> Initial solution: One static vapor bubbles for isothermal problem.\n");
      SYS_T::commPrint("     Geometry setting: [0,1]^3 cube.\n");
      break;
    case 2:
      Init_TwoBubble(pNode, locassem_ptr, fNode_ptr);
      SYS_T::commPrint("===> Initial solution: Two vapor bubbles for isothermal problem. \n");
      SYS_T::commPrint("     Geometry setting: [0,1]^3 cube.\n");
      break;
    case 3:
      Init_OneDroplet(pNode, locassem_ptr, fNode_ptr);
      SYS_T::commPrint("===> Initial solution: One droplet for isothermal problem. \n");
      SYS_T::commPrint("     Geometry setting: [0,1]^3 cube.\n");
      break;
    case 4:
      Init_8Bubble_AnnularSect(pNode, locassem_ptr, fNode_ptr);
      SYS_T::commPrint("===> Initial solution: 8 Vapor bubbles for isothermal problem. \n");
      SYS_T::commPrint("     Geometry setting: Annular w out radius 2.0, in radius 0.5,\n");
      SYS_T::commPrint("     0<z<1.\n");
      break;
    case 5:
      Init_Debug(pNode, locassem_ptr, fNode_ptr);
      SYS_T::commPrint("===> Initial solution: Debug for isothermal problem. \n");
      SYS_T::commPrint("     Geometry setting: [0,1]^3 cube.\n");
      break;
    case 6:
      Init_ThreeBubble(pNode, locassem_ptr, fNode_ptr);
      SYS_T::commPrint("===> Initial solution: Three vapor bubbles for INSK. \n");
      SYS_T::commPrint("     Geometry setting: [0,1]^3 cube.\n");
      break;
    case 7:
      Init_convTube_InOutFL(pNode, locassem_ptr, fNode_ptr);
      SYS_T::commPrint("===> Initial solution: Inflow outflow in conv Tube for INSK. \n");
      SYS_T::commPrint("     Geometry setting: converging tube, radius 1, 0.3.\n");
      break;
    default:
      SYS_T::commPrint("ERROR: PDNSolution_INSK: No such type of initial solution. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1); 
  }
}


PDNSolution_SINSK::~PDNSolution_SINSK()
{}


void PDNSolution_SINSK::Init_Zero()
{
  VecSet(solution, 0.0);
  GhostUpdate();
}


void PDNSolution_SINSK::Init_TwoBubble( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )    
{
  double Ca = locassem_ptr->get_model_para_1();
  
  double dist1, dist2, x, y, z;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  value[4] = 0.0;

  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    dist1 = sqrt((x-0.4)*(x-0.4) + (y-0.5)*(y-0.5) + (z-0.6)*(z-0.6));
    dist2 = sqrt((x-0.75)*(x-0.75) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    value[0] = 0.1 + 0.25 * (tanh((dist1 - 0.25)/(2*Ca))
        +tanh((dist2-0.1)/(2*Ca)));

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}



void PDNSolution_SINSK::Init_OneStaticBubble( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )    
{
  double Ca = locassem_ptr->get_model_para_1();
  
  double dist1, x, y, z;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  value[4] = 0.0;

  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    value[0] = 0.3545 + 0.2475 * tanh((dist1 - 0.25)/(2*Ca));

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}


void PDNSolution_SINSK::Init_OneDroplet( const class APart_Node * const &pNode_ptr,
           const class IPLocAssem * const &locassem_ptr,
                  const FEANode * const &fNode_ptr )
{
  double Ca = locassem_ptr->get_model_para_1();
  
  double dist1, x, y, z;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  value[4] = 0.0;

  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + z*z);
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    value[0] = 0.35 - 0.25 * tanh((dist1 - 0.2)/(2*Ca));

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}


void PDNSolution_SINSK::Init_8Bubble_AnnularSect( const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  double Ca = locassem_ptr->get_model_para_1();
  double inv_2Ca = 0.5 / Ca;
  double x, y, z, d;
  int num_bub = 8;

  double * xx = new double [num_bub];
  double * yy = new double [num_bub];
  double * zz = new double [num_bub];
  double * rr = new double [num_bub];

  // define the bubble centre and radii
  xx[0] = 1.05; xx[1] = 0.56; xx[2] = -1.09; xx[3] = 0.09;
  xx[4] = 1.0; xx[5] = -1.0; xx[6] = -1.0; xx[7] = 1.0;

  yy[0] = 0.07; yy[1] = 0.85; yy[2] = 0.1; yy[3] = -1.02;
  yy[4] = 0.8; yy[5] = 1.0; yy[6] = -1.0; yy[7] = -1.0;

  zz[0] = 0.5; zz[1] = 0.47; zz[2] = 0.62; zz[3] = 0.52;
  zz[4] = 0.33; zz[5] = 0.5; zz[6] = 0.59; zz[7] = 0.45;

  rr[0] = 0.20; rr[1] = 0.2; rr[2] = 0.25; rr[3] = 0.21;
  rr[4] = 0.19; rr[5] = 0.23; rr[6] = 0.21; rr[7] = 0.32;
  
  int location[5];
  double value[5];
  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  value[4] = 0.0;

  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    
    value[0] = 0.0;
    for(int jj=0; jj<num_bub; ++jj)
    {
      d = sqrt((x-xx[jj])*(x-xx[jj]) + (y-yy[jj])*(y-yy[jj]) + (z-zz[jj])*(z-zz[jj]));
      value[0] += 0.25 * tanh(inv_2Ca*(d - rr[jj]));
    }

    value[0] = value[0] - 0.25 * num_bub + 0.6;

    // check value of density
    if(value[0] < 0.09)
      value[0] = 0.107;
    if(value[0] > 0.61)
      value[0] = 0.602;

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  delete [] xx; delete [] yy; delete [] zz; delete [] rr; 
}

void PDNSolution_SINSK::Init_Debug( const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  double Ca = locassem_ptr->get_model_para_1();
  double dist1, dist2, dist3, x, y, z;
  
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  value[4] = 0.0;

  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    dist1 = sqrt( (y-0.75)*(y-0.75) + (z-0.5)*(z-0.5));
    dist2 = sqrt( (y-0.25)*(y-0.25) + (z-0.5)*(z-0.5));
    dist3 = sqrt( (y-0.4)*(y-0.4) + (z-0.75)*(z-0.75));
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    value[0] = -0.15 + 0.25 *( tanh((dist1 - 0.1)/(2.0*Ca))
        + tanh((dist2 - 0.15)/(2.0*Ca)) + tanh((dist3 - 0.08)/(2.0*Ca)) );
    

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}

void PDNSolution_SINSK::Init_ThreeBubble( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )    
{
  double Ca = locassem_ptr->get_model_para_1();
  
  double dist1, dist2, dist3, x, y, z;
  int location[5];
  double value[5];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  value[4] = 0.0;

  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    dist1 = sqrt((x+0.2)*(x+0.2) + (y-0.0)*(y-0.0) + (z-0.5)*(z-0.5));
    dist2 = sqrt((x-0.1)*(x-0.1) + (y+0.3)*(y+0.3) + (z-0.4)*(z-0.4));
    dist3 = sqrt((x-0.2)*(x-0.2) + (y-0.2)*(y-0.2) + (z-0.6)*(z-0.6));
    
    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    value[0] = -0.15 + 0.25 * ( tanh((dist1 - 0.2)/(2*Ca))
        + tanh((dist2-0.15)/(2*Ca)) + tanh((dist3-0.1)/(2*Ca)) );

    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}


void PDNSolution_SINSK::Init_convTube_InOutFL( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  double x, y, z;
  int location[5];
  double value[5];
  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  value[4] = 0.0;

  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);

    location[0] = pNode_ptr->get_node_loc(ii) * 5;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;

    value[0] = 0.6;
    value[1] = 0.0;
    value[2] = 0.0;
    if(z < 0.01)
    {
      value[3] = 1.0 - x*x - y*y;
    }
    else if(z>3.99)
    {
      value[3] = (9.0 - 100.0 * x*x - 100.0 * y*y) / 9.0;
    }
    else
    {
      value[3] = 0.0;
    }
    VecSetValues(solution, 5, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}

// EOF
