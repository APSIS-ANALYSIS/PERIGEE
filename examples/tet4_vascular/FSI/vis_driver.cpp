// ==================================================================
// vis_driver.cpp
// This is the visualization driver to output the whole FSI system
// as a single continuum.
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem_wTag.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "FEAElement_Tet4.hpp"
#include "VisDataPrep_Mixed_FSI_3D.hpp"
#include "VTK_Writer_FSI_Tet4.hpp"

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string anode_mapping_file = "node_mapping.h5";
  const std::string pnode_mapping_file = "post_node_mapping.h5";
  const std::string part_file = "postpart";

  std::string sol_bname("SOL_");
  std::string out_bname("VIS_FSI_");

  const int dof = 7;
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;
  double dt = 0.01;
  bool isXML = true;
  bool isClean = true;

  // ===== PETSc Initialization =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-clean", isClean);

  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  if(isXML) PetscPrintf(PETSC_COMM_WORLD, "-xml: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-xml: false \n");

  if(isClean) PetscPrintf(PETSC_COMM_WORLD, "-clean: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-clean: false \n");

  // If demand cleaning, remove all previous visualization files
  if( isClean )
  {
    SYS_T::execute("rm -rf VIS_FSI_*_p*.vtu");
    SYS_T::execute("rm -rf VIS_FSI_*.pvtu");
    SYS_T::execute("rm -rf VIS_FSI_.pvd");
  }

  FEANode * fNode = new FEANode(part_file, rank);
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  ALocal_Elem * locElem = new ALocal_Elem_wTag(part_file, rank);
  APart_Node * pNode = new APart_Node(part_file, rank);
  
  SYS_T::print_fatal_if(size != PartBasic->get_cpu_size(), "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);

  IQuadPts * quad = new QuadPts_vis_tet4();

  quad -> print_info();
  
  FEAElement * element = new FEAElement_Tet4( quad-> get_num_quadPts() );

  IVisDataPrep * visprep = new VisDataPrep_Mixed_FSI_3D();
  
  visprep->print_info();

  double ** pointArrays = new double * [visprep->get_ptarray_size()];
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    pointArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_ptarray_comp_length(ii)];

  VTK_Writer_FSI_Tet4 * vtk_w = new VTK_Writer_FSI_Tet4(
      GMIptr->get_nElem(), element_part_file );

  std::ostringstream time_index;

  for(int time = time_start; time<=time_end; time+= time_step)
  {
    std::string name_to_read(sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %d: Read %s and Write %s \n",
        time, name_to_read.c_str(), name_to_write.c_str() );

    visprep->get_pointArray(name_to_read, anode_mapping_file, pnode_mapping_file,
        pNode, GMIptr->get_nFunc(), dof, pointArrays);

    vtk_w->writeOutput( fNode, locIEN, locElem,
        visprep, element, quad, pointArrays, rank, size, 
        pNode -> get_ntotalnode(), 
        time * dt, sol_bname, out_bname, name_to_write, isXML );
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // Clean up memory
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] pointArrays[ii];
  delete [] pointArrays;
  delete visprep; delete element; delete quad;
  delete pNode; delete locElem; delete PartBasic; delete GMIptr;
  delete locIEN; delete fNode; delete vtk_w;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
