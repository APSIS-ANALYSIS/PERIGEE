// ==================================================================
// vis_laplace.cpp
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "QuadPts_vis_tet10_v2.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "VisDataPrep_Lap.hpp"
#include "VTK_Writer_Lap.hpp"

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string anode_mapping_file = "node_mapping.h5";
  const std::string pnode_mapping_file = "post_node_mapping.h5";
  const std::string part_file="postpart";
  const int dof = 1;

  const std::string sol_bname("SOL_900000000");
  const std::string out_bname = sol_bname;
  const bool isXML = true;
  
  // ===== Initialize MPI run =====
  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::print_fatal_if(size != PartBasic->get_cpu_size(), "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);
  
  IQuadPts * quad = nullptr;
  FEAElement * element = nullptr;

  if( GMIptr->get_elemType() == 501 )
  {
    quad = new QuadPts_vis_tet4();
    element = new FEAElement_Tet4( quad-> get_num_quadPts() );
  }
  else if( GMIptr->get_elemType() == 502 )
  {
    quad = new QuadPts_vis_tet10_v2();
    element = new FEAElement_Tet10_v2( quad-> get_num_quadPts() );
  }
  else SYS_T::print_fatal( "Error: unsupported element type \n" );

  quad -> print_info();

  IVisDataPrep * visprep = new VisDataPrep_Lap();

  visprep -> print_info();

  double ** solArrays = new double * [visprep->get_ptarray_size()];
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    solArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_ptarray_comp_length(ii)];

  VTK_Writer_Lap * vtk_w = new VTK_Writer_Lap( GMIptr->get_nElem(),
            GMIptr->get_nLocBas(), element_part_file );

  visprep->get_pointArray(sol_bname, anode_mapping_file, pnode_mapping_file,
      pNode, GMIptr->get_nFunc(), dof, solArrays);

  vtk_w->writeOutput( fNode, locIEN, locElem,
      visprep, element, quad, solArrays,
      rank, size, 0.0, sol_bname, out_bname, out_bname, isXML );

  MPI_Barrier(PETSC_COMM_WORLD);

  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] solArrays[ii];
  delete [] solArrays;

  delete fNode; delete locIEN; delete GMIptr; delete PartBasic; delete locElem;
  delete pNode; delete quad; delete element; delete visprep; delete vtk_w;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
