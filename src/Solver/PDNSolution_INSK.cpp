#include "PDNSolution_INSK.hpp"

PDNSolution_INSK::PDNSolution_INSK( const class APart_Node * const &pNode,
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
    default:
      SYS_T::commPrint("ERROR: PDNSolution_INSK: No such type of initial solution. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1); 
  }
}


PDNSolution_INSK::~PDNSolution_INSK()
{}


void PDNSolution_INSK::Init_Zero()
{
  VecSet(solution, 0.0);
  GhostUpdate();
}


void PDNSolution_INSK::Init_TwoBubble( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )    
{
  double Ca = locassem_ptr->get_model_para_1();
  
  double dist1, dist2, x, y, z;
  int location[4];
  double value[4];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  
  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    dist1 = sqrt((x-0.4)*(x-0.4) + (y-0.5)*(y-0.5) + (z-0.6)*(z-0.6));
    dist2 = sqrt((x-0.75)*(x-0.75) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
    
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    value[0] = 0.1 + 0.25 * (tanh((dist1 - 0.25)/(2*Ca))
        +tanh((dist2-0.1)/(2*Ca)));

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}



void PDNSolution_INSK::Init_OneStaticBubble( 
    const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )    
{
  double Ca = locassem_ptr->get_model_para_1();
  
  double dist1, x, y, z;
  int location[4];
  double value[4];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  
  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
    
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    value[0] = 0.3545 + 0.2475 * tanh((dist1 - 0.25)/(2*Ca));

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}


void PDNSolution_INSK::Init_OneDroplet( const class APart_Node * const &pNode_ptr,
           const class IPLocAssem * const &locassem_ptr,
                  const FEANode * const &fNode_ptr )
{
  double Ca = locassem_ptr->get_model_para_1();
  
  double dist1, x, y, z;
  int location[4];
  double value[4];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  
  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + z*z);
    
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    value[0] = 0.35 - 0.25 * tanh((dist1 - 0.2)/(2*Ca));

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}


void PDNSolution_INSK::Init_8Bubble_AnnularSect( const class APart_Node * const &pNode_ptr,
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

  rr[0] = 0.16; rr[1] = 0.17; rr[2] = 0.18; rr[3] = 0.2;
  rr[4] = 0.20; rr[5] = 0.28; rr[6] = 0.18; rr[7] = 0.17;
  
  int location[4];
  double value[4];
  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;

  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    x = fNode_ptr->get_ctrlPts_x(ii);
    y = fNode_ptr->get_ctrlPts_y(ii);
    z = fNode_ptr->get_ctrlPts_z(ii);
    
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    
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

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  delete [] xx; delete [] yy; delete [] zz; delete [] rr; 
}


void PDNSolution_INSK::Init_Debug( const class APart_Node * const &pNode_ptr,
    const class IPLocAssem * const &locassem_ptr,
    const FEANode * const &fNode_ptr )
{
  //double Ca = locassem_ptr->get_model_para_1();
  
  //double dist1, x, y, z;
  int location[4];
  double value[4];

  value[0] = 0.0; value[1] = 0.0; value[2] = 0.0; value[3] = 0.0;
  
  for(int ii=0; ii<pNode_ptr->get_nlocalnode(); ++ii)
  {
    //x = fNode_ptr->get_ctrlPts_x(ii);
    //y = fNode_ptr->get_ctrlPts_y(ii);
    //z = fNode_ptr->get_ctrlPts_z(ii);
    
    //dist1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
    
    location[0] = pNode_ptr->get_node_loc(ii) * 4;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;

    //value[0] = 0.3545 + 0.2475 * tanh((dist1 - 0.25)/(2*Ca));
    value[0] = 0.6;

    VecSetValues(solution, 4, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
}

// EOF
