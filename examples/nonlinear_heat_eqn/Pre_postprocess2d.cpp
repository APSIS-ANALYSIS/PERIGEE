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
//         A. mesh parameter:
//         -elemx, -elemy
//         -addSDegree, -addTDegree
//         -geo_file
//        
//        They should be the same as the input in the preprocessor.
//        
//
//        B. Partition parameter:
//         -cpu_size, -in_ncommon, -METIS_isDualGraph
//        
//        They can be chosen differently from analysis preprocessor.
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
// Date: April 18 2014
// ==================================================================
#include "Sys_Tools.hpp"
#include "kRefinement.hpp"
#include "IMesh.hpp"
#include "Mesh_NURBS_1Patch_2D.hpp"
#include "IEN_NURBS_1Patch_2D.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_NURBS_1Patch_2D_METIS.hpp"
#include "HDF5_PartReader.hpp"
#include "petscsys.h"
using namespace std;

int main(int argc, char * argv[])
{
  vector<double> sKnots, tKnots, uKnots;
  int sDegree, tDegree, uDegree, numCPts;
  vector<double> ctrlPts;
  const int cpDim = 4;
  const int spatialDim = 3;

  // these three are not needed, just to fill in function call parameters
  int probDim = 2;
  int dofNum  = 1;
  int elemType = 0;

  // refinement
  int addSDegree = 1, addTDegree = 1;

  int num_inserted_x = 8, num_inserted_y = 8;

  // geometry file
  char * char_home_dir = getenv("HOME");
  string geo_file(char_home_dir);
  geo_file.append("/PERIGEE/input/geometry_2d_square.txt");

  // partition file basename
  string part_file("postpart");

  // partition parameter
  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;
  
  // flag to check compatibility of mesh
  bool isread_part = true;

  PetscMPIInt rank, size;

  // ======= PETSc Initialization =======
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // This partition has to be serial
  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

  // Read in command line arguments
  SYS_T::GetOptionInt("-elemx", num_inserted_x);
  SYS_T::GetOptionInt("-elemy", num_inserted_y);
  SYS_T::GetOptionInt("-addSDegree", addSDegree);
  SYS_T::GetOptionInt("-addTDegree", addTDegree);

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);

  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-part_file", part_file);

  SYS_T::GetOptionBool("-isread_part", isread_part);

  // Print command line argument on screen
  cout<<"==== /Command Line Arguments ===="<<endl;
  cout<<" -elemx: "<<num_inserted_x<<endl;
  cout<<" -elemy: "<<num_inserted_y<<endl;
  cout<<" -addSDegree: "<<addSDegree<<endl;
  cout<<" -addTDegree: "<<addTDegree<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph)
    cout<<" -isDualGraph: true \n";
  else
    cout<<" -isDualGraph: false \n";
  if(isread_part)
    cout<<" -isread_part: true \n";
  else
    cout<<" -isread_part: false \n";
  cout<<"----------------------------------\n";
  cout<<"probDim: "<<probDim<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // ======= Read in geometry file =======
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

  if(uDegree != 0 || uKnots.size() > 1 )
  {
    cerr<<"ERROR: The input geometry file is not 2D. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  NURBS_T::projectUp(ctrlPts, spatialDim);

  // refinement
  int num_knotspan_s = knotVec_check(sKnots, sDegree);
  vector<int> num_insert_s;
  for(int ii=0; ii<num_knotspan_s; ++ii)
    num_insert_s.push_back(num_inserted_x);

  int num_knotspan_t = knotVec_check(tKnots, tDegree);
  vector<int> num_insert_t;
  for(int ii=0; ii<num_knotspan_t; ++ii)
    num_insert_t.push_back(num_inserted_y);

  vector<double> insertKnots_s, insertKnots_t;
  double hs_max, ht_max, hs_min, ht_min;
  hRefine_newKnot_Generator(sKnots, insertKnots_s, num_insert_s, hs_max, hs_min );
  hRefine_newKnot_Generator(tKnots, insertKnots_t, num_insert_t, ht_max, ht_min );

  kRefinement(addSDegree, addTDegree, insertKnots_s, insertKnots_t,
      sKnots, tKnots, ctrlPts, cpDim, sDegree, tDegree );

  cout<<endl<<"=== New knot vector generated. \n";

  // Generate extraction operators
  vector<NURBS_T::BezierElem*> seg_x, seg_y;
  vector<double> fake_ctrl, fake_bCtrlPts;
  NURBS_T::decomposeCurve(sKnots, sDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_x, false );
  NURBS_T::decomposeCurve(tKnots, tDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_y, false );

  cout<<endl<<"=== Extraction operators generated. \n";

  // Project down the control points
  NURBS_T::projectDown( ctrlPts, spatialDim );

  // Mesh object generation
  IMesh * Mesh;
  IIEN * IEN;

  Mesh = new Mesh_NURBS_1Patch_2D( sDegree, tDegree,
      hs_max, ht_max, hs_min, ht_min, sKnots, tKnots );

  IEN = new IEN_NURBS_1Patch_2D(Mesh);

  if(isread_part)
  {
    cout<<"\n=== Check the compatibility of mesh ...\n";
    cout<<"-- Read from part_p00000.h5 \n";
    string partfilebasename = "part";
    HDF5_PartReader * check_reader = new HDF5_PartReader(partfilebasename, 0);

    int check_sd, check_td;
    check_reader->get_GMI_degree(check_sd, check_td);
    assert(check_sd == sDegree);
    assert(check_td == tDegree);

    cout<<"-- Degrees are compatible... \n";

    int check_ne, check_nex, check_ney;
    check_reader->get_GMI_nElem(check_ne, check_nex, check_ney);
    assert(check_ne == Mesh->get_nElem());
    assert(check_nex == Mesh->get_nElem_x());
    assert(check_ney == Mesh->get_nElem_y());

    cout<<"-- nElem, nElem_x/y are compatible. "<<endl;
    delete check_reader;
  }

  Mesh->print_info();

  // Global mesh partition 
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

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, Mesh->get_nFunc());
  mnindex->write_hdf5("post_node_mapping");

  cout<<"\n=== Start mesh partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_NURBS_1Patch_2D_METIS(Mesh, global_part, mnindex, IEN,
        ctrlPts, seg_x, seg_y, proc_rank, proc_size, dofNum, elemType,
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

  delete Mesh;
  delete IEN;
  delete global_part;
  delete mnindex;

  PetscFinalize();
  return 0;
}

// EOF
