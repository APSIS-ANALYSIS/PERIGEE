// ==================================================================
// vis_3DNLHeat.cpp
// ------------------------------------------------------------------
// This is the visualization driver to visualize 3DNLHeat solution
// in VTK format.
//
// This is a parallel routine. Users have to run Pre_postprocess to 
// obtain a parallel mesh partition.
//
// Date: Dec. 12 2013
// ==================================================================
#include <cmath>
#include <sstream>

#include "Sys_Tools.hpp"
#include "QuadPts_vis.hpp"
#include "HDF5_PartReader.hpp"
#include "BernsteinBasis_Array.hpp"
#include "FEANode.hpp"
#include "AExtractor_3D_NURBS_xyz.hpp"
#include "FEAElement_NURBS_3D_der0_v3.hpp"
#include "FEAElement_NURBS_3D_der1_v3.hpp"
#include "APart_Node.hpp"
#include "AGlobal_Mesh_Info_1Patch_NURBS_3D.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_meshSize_3D_NURBS.hpp"
#include "APart_Basic_Info.hpp"
#include "PostVectSolution.hpp"
#include "Post_error.hpp"
#include "IVisDataPrep.hpp"
#include "VisDataPrep_3DNLHeat.hpp"
#include "VTK_Writer.hpp"

using namespace std;

int main( int argc, char * argv[] )
{
  string element_part_file = "epart.h5";
  string anode_mapping_file = "node_mapping.h5";
  string pnode_mapping_file = "post_node_mapping.h5";

  // Number of visualization sampling points in each finite element in
  // each direction
  int samx = 3; int samy = 3; int samz = 3;

  // Solution base name, default SOL
  string sol_bname("SOL_");

  // Output base name, default is the same as sol_bname
  string out_bname = sol_bname;

  // Degree of freedom, it will be read from hdf5
  int dof;

  // Solution time
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;
  double dt = 0.01;

  // vtk format
  bool isXML = true;

  // partition file base name
  string part_file("postpart");

  PetscMPIInt rank, size;
  // ====== PETSc Initialize =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ====== Read command line =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  
  SYS_T::GetOptionInt("-samx", samx);
  SYS_T::GetOptionInt("-samy", samy);
  SYS_T::GetOptionInt("-samz", samz);
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionBool("-xml", isXML);

  samx = samx > 1 ? samx : 2;
  samy = samy > 1 ? samy : 2;
  samz = samz > 1 ? samz : 2;

  SYS_T::cmdPrint("-samx:", samx);
  SYS_T::cmdPrint("-samy:", samy);
  SYS_T::cmdPrint("-smaz:", samz);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  if(isXML == true)
    PetscPrintf(PETSC_COMM_WORLD, "-xml: true \n");
  else
    PetscPrintf(PETSC_COMM_WORLD, "-xml: false \n");

  // ===== Read Partition file =====
  HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);
  SYS_T::commPrint("===> Reading mesh files ... \n");
  // get global mesh info
  IAGlobal_Mesh_Info * gInfo_ptr = new AGlobal_Mesh_Info_1Patch_NURBS_3D(h5reader);

  // get node partition info
  APart_Node * pNode = new APart_Node(part_file, rank);  

  // get degree-of-freedom from partition files
  h5reader->get_GMI_dofNum( dof );

  // get control points
  FEANode * fNode = new FEANode(part_file, rank);

  // get mesh size
  IALocal_meshSize * locmSize = new ALocal_meshSize_3D_NURBS(h5reader);

  // get extraction operator
  IAExtractor * fExt = new AExtractor_3D_NURBS_xyz(h5reader);

  // get LIEN
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // get local element info
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // get partition info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  SYS_T::synPrint(" done\t", rank);
  delete h5reader;

  if(size != PartBasic->get_cpu_size())
    MPI_Abort(PETSC_COMM_WORLD, 1);

  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for:", size);
  PetscPrintf(PETSC_COMM_WORLD, "Postprocessing - visualization.\n");

  // ===== Prepare for building sampling points in Bernstein basis 
  SYS_T::commPrint("\n===> Build sampling points ... \n");
  IQuadPts * quad_z = new QuadPts_vis(samz);
  IQuadPts * quad_y = new QuadPts_vis(samy);
  IQuadPts * quad_x = new QuadPts_vis(samx);

  SYS_T::commPrint("===> Build univariate Bezier elements ... \n\n");
  BernsteinBasis_Array Bena_x(gInfo_ptr->get_xdegree(), quad_x);
  BernsteinBasis_Array Bena_y(gInfo_ptr->get_ydegree(), quad_y);
  BernsteinBasis_Array Bena_z(gInfo_ptr->get_zdegree(), quad_z);

  // ----------------------------------------------------------------
  //                --- Start Post-processing ---
  // ----------------------------------------------------------------
  // Data name and size setting
  IVisDataPrep * visprep = new VisDataPrep_3DNLHeat();
  PetscPrintf(PETSC_COMM_WORLD, "======================================= \n");
  PetscPrintf(PETSC_COMM_WORLD, "Data to be visualized: \n");
  PetscPrintf(PETSC_COMM_WORLD, "-- %d type(s) of data. \n", visprep->get_arrayCompSize());
  for(int ii=0; ii<visprep->get_arrayCompSize(); ++ii)
  {
    string visprep_temp_name = visprep->get_arrayNames(ii);
    PetscPrintf(PETSC_COMM_WORLD, "-- %s \t", visprep_temp_name.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "with size %d \n", visprep->get_arraySizes(ii));
  }
  PetscPrintf(PETSC_COMM_WORLD, "======================================= \n");

  // Allocate memory for pointArrays
  double ** pointArrays = new double * [visprep->get_arrayCompSize()];
  for(int ii=0; ii<visprep->get_arrayCompSize(); ++ii)
    pointArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_arraySizes(ii)];

  // Allocate the class of vtk writer
  VTK_Writer * vtk_w = new VTK_Writer(gInfo_ptr);

  ostringstream time_index;
  // Loop over time to write solutions
  for(int time = time_start; time<=time_end; time = time + time_step)
  {
    // prepare solution full names
    string name_to_read(sol_bname);
    string name_to_write(out_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());
    PetscPrintf(PETSC_COMM_WORLD, "Time %d: Read %s and Write %s \n",
        time, name_to_read.c_str(), name_to_write.c_str() );

    // use visprep to prepare pointArrays
    visprep->get_pointArray(name_to_read, anode_mapping_file, pnode_mapping_file,
        pNode, gInfo_ptr, dof, pointArrays);

    // call writer to write vtk files
    vtk_w->writeOutput(element_part_file, gInfo_ptr, fNode, locmSize, fExt, locIEN,
        locElem, visprep, &Bena_x, &Bena_y, &Bena_z, pointArrays, rank, size,
        time * dt, sol_bname, out_bname, name_to_write, isXML  );
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  // ===== PETSc Finalize =====
  delete pNode;
  delete gInfo_ptr;
  delete fNode;
  delete locmSize;
  delete fExt;
  delete locIEN;
  delete locElem;
  delete PartBasic;
  delete quad_x; delete quad_y; delete quad_z;
  delete vtk_w;
  for(int ii=0; ii<visprep->get_arrayCompSize(); ++ii)
    delete [] pointArrays[ii];
  delete [] pointArrays;
  delete visprep;
  PetscFinalize();
  return 0;
}

// EOF
