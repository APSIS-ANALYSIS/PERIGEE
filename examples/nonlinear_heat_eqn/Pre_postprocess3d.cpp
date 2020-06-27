// ==================================================================
// Pre_postprocess.cpp
// ------------------------------------------------------------------
// Objective:
// This routine provides a mesh partition for all postprocessers. This
// routine is needed because the postprocessors should be run in 
// parallel. But we do not require the postprocessor use the same mesh
// partition as the analysis part, since most of the time, postprocess
// requires less cpu.
//
// Usage: We requires the input of the mesh infomation from command
//         line, including:
//         -elemx, -elemy, -elemz
//         -addSDegree, -addTDegree, -addUDegree
//         -geo_file
//        should be the same as the input in the preprocessor.
//        The partition parameters, including:
//         -cpu_size, -in_ncommon, -METIS_isDualGraph
//        can be chosen differently.
//
//        The partition file name -part_file should be chosen as a 
//        different name to avoid error in postprocess main drivers.
//
// Output:
// 1. Partitioned mesh infomation for each processor;
// 2. post_epart.h5 post_npart.h5 files that describes the element
//    and node partition in postprocessors;
// 3. post_node_mapping.h5: the old_2_new and new_2_old mapping for
//    the postprocessors.
//
// Date: Dec. 10th 2013
// ==================================================================
#include <cassert>

#include "Sys_Tools.hpp"
#include "NURBS_Tools.hpp"
#include "kRefinement.hpp"
#include "Mesh_NURBS_1Patch_3D.hpp"
#include "IEN_NURBS_1Patch_3D.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Map_Node_Index.hpp"
#include "Part_NURBS_1Patch_3D_METIS.hpp"
#include "HDF5_PartReader.hpp"
using namespace std;

int main(int argc, char *argv[])
{
  vector<double> sKnots, tKnots, uKnots;
  int sDegree, tDegree, uDegree, numCPts;
  vector<double> ctrlPts;
  const int cpDim = 4;
  const int spatialDim = 3;

  // problem description
  //  -- I donot like this part
  int probDim = 3;
  int dofNum = 1;
  int elemType = 0;

  // degree to be added
  int addSDegree = 1, addTDegree = 1, addUDegree = 1;

  // number of knots to be inserted
  int num_inserted_x = 8, num_inserted_y = 8, num_inserted_z = 8;

  // geometry file
  char * char_home_dir = getenv("HOME");
  string geo_file(char_home_dir);
  geo_file.append("/PERIGEE/input/geometry_3d_cube.txt");

  // partition file base_name
  string part_file("postpart");

  // partition parameters
  int cpu_size = 1;
  int in_ncommon = 9;
  bool isDualGraph = true;

  // flag to check the command line argument
  bool isread_part = true;

  PetscMPIInt rank, size;
  // ======= PETSc Initialization ======
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // This partition has to be serial
  if(size != 1)
  {
    PetscPrintf(PETSC_COMM_WORLD, "ERROR: This partition code has to be run in serial! \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  SYS_T::GetOptionInt("-elemx", num_inserted_x);
  SYS_T::GetOptionInt("-elemy", num_inserted_y);
  SYS_T::GetOptionInt("-elemz", num_inserted_z);
  SYS_T::GetOptionInt("-addSDegree", addSDegree);
  SYS_T::GetOptionInt("-addTDegree", addTDegree);
  SYS_T::GetOptionInt("-addUDegree", addUDegree);

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);

  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);

  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-part_file", part_file);

  // Print command line on screen
  cout<<"==== /Command Line Arguments ===="<<endl;
  cout<<" -elemx: "<<num_inserted_x<<endl;
  cout<<" -elemy: "<<num_inserted_y<<endl;
  cout<<" -elemz: "<<num_inserted_z<<endl;
  cout<<" -addSDegree: "<<addSDegree<<endl;
  cout<<" -addTDegree: "<<addTDegree<<endl;
  cout<<" -addUDegree: "<<addUDegree<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph)
    cout<<" -isDualGraph: true \n";
  else
    cout<<" -isDualGraph: false \n";
  
  cout<<"----------------------------------\n";
  cout<<"probDim: "<<probDim<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // Read in the geometry file
  ifstream infile( geo_file.c_str(), ifstream::in );
  if(infile.is_open() == false)
  {
    cerr<<"ERROR: Can not find file: "<<geo_file<<endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  else
  {
    SYS_T::readFile(infile, sKnots, tKnots, uKnots,
        sDegree, tDegree, uDegree, numCPts, ctrlPts);
  }
  infile.close();

  NURBS_T::projectUp(ctrlPts, spatialDim);

  int num_knotspan_s = knotVec_check(sKnots, sDegree);
  vector<int> num_insert_s;
  for(int ii=0; ii<num_knotspan_s; ++ii)
    num_insert_s.push_back(num_inserted_x);

  int num_knotspan_t = knotVec_check(tKnots, tDegree);
  vector<int> num_insert_t;
  for(int ii=0; ii<num_knotspan_t; ++ii)
    num_insert_t.push_back(num_inserted_y);

  int num_knotspan_u = knotVec_check(uKnots, uDegree);
  vector<int> num_insert_u;
  for(int ii=0; ii<num_knotspan_u; ++ii)
    num_insert_u.push_back(num_inserted_z);

  vector<double> insertKnots_s, insertKnots_t, insertKnots_u;

  double hs_max, ht_max, hu_max, hs_min, ht_min, hu_min;
  hRefine_newKnot_Generator(sKnots, insertKnots_s, num_insert_s, hs_max, hs_min );
  hRefine_newKnot_Generator(tKnots, insertKnots_t, num_insert_t, ht_max, ht_min );
  hRefine_newKnot_Generator(uKnots, insertKnots_u, num_insert_u, hu_max, hu_min );


  kRefinement(addSDegree, addTDegree, addUDegree, insertKnots_s, insertKnots_t,
      insertKnots_u, sKnots, tKnots, uKnots, ctrlPts, cpDim, sDegree, tDegree, uDegree );

  cout<<endl<<"=== New knot vectors generated. \n";

  // Generate extraction operator 
  vector<NURBS_T::BezierElem*> seg_x, seg_y, seg_z;
  vector<double> fake_ctrl, fake_bCtrlPts;
  NURBS_T::decomposeCurve(sKnots, sDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_x, false );
  NURBS_T::decomposeCurve(tKnots, tDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_y, false );
  NURBS_T::decomposeCurve(uKnots, uDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_z, false );

  cout<<endl<<"=== Extraction operators generated. \n";

  // Project down control points
  NURBS_T::projectDown(ctrlPts, spatialDim);

  // Global Mesh Generation
  IMesh * Mesh = new Mesh_NURBS_1Patch_3D( sDegree, tDegree, uDegree,
      hs_max, ht_max, hu_max, hs_min, ht_min, hu_min, sKnots, tKnots, uKnots );
  IIEN * IEN = new IEN_NURBS_1Patch_3D(Mesh);

  if(isread_part)
  {
    cout<<"\n=== Check the compatibility ...\n";
    cout<<"-- Read from part_p00000.h5 \n";
    string partfilebasename = "part";
    HDF5_PartReader * check_reader = new HDF5_PartReader(partfilebasename, 0);

    int check_sd, check_td, check_ud;
    check_reader->get_GMI_degree(check_sd, check_td, check_ud);
    assert(check_sd == sDegree);
    assert(check_td == tDegree);
    assert(check_ud == uDegree);

    cout<<"-- Degrees are compatible... \n";

    int check_ne, check_nex, check_ney, check_nez;
    check_reader->get_GMI_nElem(check_ne, check_nex, check_ney, check_nez);
    assert(check_ne == Mesh->get_nElem());
    assert(check_nex == Mesh->get_nElem_x());
    assert(check_ney == Mesh->get_nElem_y());
    assert(check_nez == Mesh->get_nElem_z());

    cout<<"-- nElem, nElem_x/y/z are compatible. "<<endl;
    delete check_reader;
  }

  Mesh->print_info();

  // Generate global mesh partition files
  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, Mesh, IEN, "post_epart", "post_npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( Mesh, "post_epart", "post_npart" );
  else
  {
    cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<endl;
    exit(1);
  }

  // node reorder numbering
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, Mesh->get_nFunc());

  mnindex->write_hdf5("post_node_mapping");

  // Generate partition files
  cout<<"\n=== Start Partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_NURBS_1Patch_3D_METIS(Mesh, global_part, mnindex, IEN,
        ctrlPts, seg_x, seg_y, seg_z, proc_rank, proc_size, dofNum, elemType,
        isPrintPartInfo );

    part->write(part_file.c_str());

    delete part;
  }

  // Clean objects in memory 
  cout<<endl<<"=== Clean memory. \n";
  vector<NURBS_T::BezierElem*>::iterator it_bezier;
  for( it_bezier = seg_x.begin(); it_bezier != seg_x.end(); ++it_bezier )
    delete *it_bezier;
  for( it_bezier = seg_y.begin(); it_bezier != seg_y.end(); ++it_bezier )
    delete *it_bezier;
  for( it_bezier = seg_z.begin(); it_bezier != seg_z.end(); ++it_bezier )
    delete *it_bezier;

  delete IEN;
  delete Mesh;
  delete global_part;
  delete mnindex;

  PetscFinalize();
  return 0;
}


// EOF
