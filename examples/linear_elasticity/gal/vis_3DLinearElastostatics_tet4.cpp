// ==================================================================
// vis_3DLinearElastostatics_tet4.cpp
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "APart_Node.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "FEAElement_Tet4.hpp"
#include "VisDataPrep_LinearElastic_3D.hpp"
#include "VTK_Writer_LE_Tet4.hpp"

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string anode_mapping_file = "node_mapping.h5";
  const std::string pnode_mapping_file = "post_node_mapping.h5";

  // Solution name
  std::string sol_bname("SOL_");

  // Output base name
  std::string out_bname = sol_bname;

  // Degrees of freedom
  const int dof = 3;

  // vtk format
  bool isXML = true;

  // plot on ref. configuration
  bool isRef = false;

  // partition file base name
  std::string part_file("postpart");

  PetscMPIInt rank, size;
  // ===== PETSc Initialization =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-ref", isRef);

  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  if(isXML) PetscPrintf(PETSC_COMM_WORLD, "-xml: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-xml: false \n");

  if(isRef) PetscPrintf(PETSC_COMM_WORLD, "-ref: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-ref: false \n");

  SYS_T::commPrint("===> Reading mesh files ... ");

  // Control points
  FEANode * fNode = new FEANode(part_file, rank);

  // LIEN
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // Global mesh info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);

  // Partitioning info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  // Local element info
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // Node partitioning info
  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::commPrint("Done! \n");

  if(size != PartBasic->get_cpu_size()) SYS_T::print_fatal(
      "Error: number of processors does not match with prepost! \n");

  PetscPrintf(PETSC_COMM_WORLD, 
      "\n===> %d processor(s) are assigned for:", size);
  PetscPrintf(PETSC_COMM_WORLD, "Postprocessing - visualization.\n");

  SYS_T::commPrint("===> Build sampling points.");
  IQuadPts * quad = new QuadPts_vis_tet4();

  quad -> print_info();

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

  // ----------------------------------------------------------------
  // Start preparing and writing vtk files
  // ----------------------------------------------------------------
  IVisDataPrep * visprep = new VisDataPrep_LinearElastic_3D();
  visprep->print_info();

  double ** pointArrays = new double * [visprep->get_ptarray_size()];
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    pointArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_ptarray_comp_length(ii)];

  // Setup writer context
  VTK_Writer_LE_Tet4 * vtk_w_s = new VTK_Writer_LE_Tet4( 
      GMIptr->get_nElem(), element_part_file );  

  std::ostringstream time_index;

  const int time = 0;
  const double dt = 0.0;
  std::string name_to_read(sol_bname);
  std::string name_to_write(out_bname);
  time_index.str("");
  time_index<< 900000000 + time;
  name_to_read.append(time_index.str());
  name_to_write.append(time_index.str());

  PetscPrintf(PETSC_COMM_WORLD, "Time %d: Read %s and Write %s \n",
      time, name_to_read.c_str(), name_to_write.c_str() );

  visprep->get_pointArray(name_to_read, anode_mapping_file, pnode_mapping_file,
      pNode, GMIptr, dof, pointArrays);

  if( isRef )
    vtk_w_s->writeOutput_ref( fNode, locIEN, locElem,
        visprep, element, quad, pointArrays,
        rank, size, time * dt, sol_bname, out_bname, name_to_write, isXML );
  else
    vtk_w_s->writeOutput( fNode, locIEN, locElem,
        visprep, element, quad, pointArrays,
        rank, size, time * dt, sol_bname, out_bname, name_to_write, isXML );

  MPI_Barrier(PETSC_COMM_WORLD);

  // ==== Finalize ====
  delete vtk_w_s;
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] pointArrays[ii];
  delete [] pointArrays;
  delete visprep; delete element; delete quad;
  delete pNode; delete locElem; delete PartBasic; delete GMIptr;
  delete locIEN; delete fNode;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
