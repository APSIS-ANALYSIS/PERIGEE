// ==================================================================
// vis_3DElastoDynamics_tet4.cpp
// ==================================================================
#include "Sys_Tools.hpp" 
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "APart_Node.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "QuadPts_vis_tet10.hpp"
#include "MaterialModel_NeoHookean_ST91_Mixed.hpp"
#include "MaterialModel_GOH06_Incompressible_Mixed.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10.hpp"
#include "VisDataPrep_Mixed_U_Hyperelastic_3D.hpp"
#include "VTK_Writer_Solids_Tet4.hpp"

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
  const int dof = 7;

  // Solution time info
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;
  double dt = 0.01;

  // vtk format
  bool isXML = true;

  // plot on ref. configuration
  bool isRef = false;

  // flag to determine if previous vis files need to be removed
  bool isRestart = false;

  // partition file base name
  std::string part_file("postpart");

  PetscMPIInt rank, size;

  // ===== PETSc Initialization =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-ref", isRef);
  SYS_T::GetOptionBool("-restart", isRestart);

  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  
  if(isXML) PetscPrintf(PETSC_COMM_WORLD, "-xml: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-xml: false \n");

  if(isRef) PetscPrintf(PETSC_COMM_WORLD, "-ref: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-ref: false \n");
  
  if( !isRestart )
  {
    int sysret = system("rm -rf *_p*.vtu");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
    sysret = system("rm -rf *.pvtu");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
    sysret = system("rm -rf *_.pvd");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  }

  SYS_T::commPrint("===> Reading mesh files ... ");

  // Control points
  FEANode * fNode = new FEANode(part_file, rank);

  // LIEN
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // Global mesh info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file, rank);

  // Partitioning info
  APart_Basic_Info * PartBasic = new APart_Basic_Info( part_file, rank );

  // Local element info
  ALocal_Elem * locElem = new ALocal_Elem( part_file, rank );

  // Node partitioning info
  APart_Node * pNode = new APart_Node( part_file, rank );

  SYS_T::commPrint("Done! \n");

  if(size != PartBasic->get_cpu_size()) SYS_T::print_fatal(
      "Error: number of processors does not match with prepost! \n");

  PetscPrintf(PETSC_COMM_WORLD, 
      "===> %d processor(s) are assigned for:", size);
  PetscPrintf(PETSC_COMM_WORLD, "Postprocessing - visualization.\n");

  SYS_T::commPrint("===> Build sampling points.");
  IQuadPts * quad;

  if( GMIptr->get_elemType() == 531 )
    quad = new QuadPts_vis_tet4();
  else if( GMIptr->get_elemType() == 532 )
    quad = new QuadPts_vis_tet10();
  else
    SYS_T::print_fatal("Error: The element type is not implemented. \n");

  quad -> print_info();

  SYS_T::commPrint("===> Setup the Material model.\n");
  //const double fangle = 90 - 49.98;
  //IMaterialModel * matmodel = new MaterialModel_GOH06_Incompressible_Mixed(
  //    7.64, fangle, 0.0, (-1.0)*fangle, 0.0, 996.6, 524.6, 0.0 );

  const double mat_in_r = 1.0;
  const double mat_in_E = 2.40566e8;
  const double mat_in_nu = 0.4999;
  IMaterialModel * matmodel = new MaterialModel_NeoHookean_ST91_Mixed(
      mat_in_r, mat_in_E, mat_in_nu );

  matmodel->print_info();

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * element;
  
  if( GMIptr->get_elemType() == 531 )
    element = new FEAElement_Tet4( quad-> get_num_quadPts() );
  else if( GMIptr->get_elemType() == 532 )
    element = new FEAElement_Tet10( quad-> get_num_quadPts() );

  // ----------------------------------------------------------------
  // Start preparing and writing vtk files
  // ----------------------------------------------------------------
  IVisDataPrep * visprep = new VisDataPrep_Mixed_U_Hyperelastic_3D( isRef );
  visprep->print_info();

  double ** pointArrays = new double * [visprep->get_ptarray_size()];
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    pointArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_ptarray_comp_length(ii)];

  // Setup writer context
  VTK_Writer_Solids_Tet4 * vtk_w_s;
  
  if( GMIptr->get_elemType() == 531 )
    vtk_w_s = new VTK_Writer_Solids_Tet4(GMIptr->get_nElem(),element_part_file);  
  else if( GMIptr->get_elemType() == 532 )
    vtk_w_s = new VTK_Writer_Solids_Tet4(GMIptr->get_nElem(),element_part_file,10);  

  std::ostringstream time_index;

  for(int time = time_start; time<=time_end; time+= time_step)
  {
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
          visprep, matmodel, element, quad, pointArrays,
          rank, size, time * dt, out_bname, name_to_write, isXML );
    else 
      vtk_w_s->writeOutput_cur_isotropic( fNode, locIEN, locElem,
          visprep, matmodel, element, quad, pointArrays,
          rank, size, time * dt, out_bname, name_to_write, isXML );
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ==== Finalize ====
  delete vtk_w_s;
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] pointArrays[ii];
  delete [] pointArrays;
  delete visprep; delete element; delete matmodel; delete quad;
  delete pNode; delete locElem; delete PartBasic; delete GMIptr;
  delete locIEN; delete fNode;
  PetscFinalize();
  return 0;
}


// EOF
