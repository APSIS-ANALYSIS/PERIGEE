// =======================================================================
// PREPROCESS3D_MAIN.CPP
// -----------------------------------------------------------------------
// Objective:
// This preprocess code is for mesh generation and partition for 3d
// isogeometric analysis. 
//
// Output: 
// 1. Partitioned mesh information for each processor.
// 2. vtk file for the geometry for initial visualization
//
// Note: 
// 1. Preprocess code is a serial code.
// 2. This one should be used for nonlinear heat equation only.
//    The template for preprocess should be found in 
//    /trunk/preprocess3d_main.cpp
//
// Date: Jan 05 2014
// =======================================================================
#include "Sys_Tools.hpp"
#include "NURBS_Bezier.hpp"
#include "NURBS_Tools.hpp"
#include "NURBS_FileIO.hpp"
#include "HDF5_Writer.hpp"
#include "kRefinement.hpp"
#include "Mesh_NURBS_1Patch_3D.hpp"
#include "IEN_NURBS_1Patch_3D.hpp"
#include "BoundaryCond.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Map_Node_Index.hpp"
#include "BC_Partition.hpp"
#include "Part_NURBS_1Patch_3D_METIS.hpp"

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  std::vector<double> sKnots, tKnots, uKnots;
  int sDegree, tDegree, uDegree, numCPts;
  std::vector<double> ctrlPts;
  const int cpDim = 4;
  const int spatialDim = 3;
  
  // problem description
  int probDim  = 3;  // the problem dimension
  int dofNum   = 1;  // the number of degree of freedoms
  int elemType = 0; // types of elements

  // Degree to be added in s,t,u directions
  int addSDegree=1, addTDegree = 1, addUDegree = 1;
  
  // Num of knots inserted in spatial discretizations
  int num_inserted_x = 8, num_inserted_y = 8, num_inserted_z = 8;

  char * char_home_dir = getenv("HOME");
  std::string geo_file(char_home_dir);
  geo_file.append("/PERIGEE/input/geometry_3d_cube.txt");

  // partition file name 
  std::string part_file("part");

  // partition parameters
  int cpu_size = 1;
  int in_ncommon = 9;
  bool isDualGraph = true;
  bool isWriteCNet = true;

  PetscMPIInt rank, size;
  // ====== PETSc Initialize ======
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Enforce serial run.
  if(size != 1)
  {
    cout<<"ERROR: Given processor number is greater than 1."; 
    cout<<"Preprocess code has to be serial!"<<endl;
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
  SYS_T::GetOptionBool("-isWriteControlNet", isWriteCNet);

  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-part_file", part_file);

  // Print command line options on screen
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
  if(isWriteCNet)
    cout<<" -isWriteControlNet: true \n";
  else
    cout<<" -isWriteControlNet: false \n";
  cout<<"----------------------------------\n";
  cout<<"probDim: "<<probDim<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // Read in Geometry file
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

  // k-refinement
  int num_knotspan_s = knotVec_check(sKnots, sDegree);
  std::vector<int> num_insert_s;
  for(int ii=0; ii<num_knotspan_s; ++ii)
    num_insert_s.push_back(num_inserted_x);

  int num_knotspan_t = knotVec_check(tKnots, tDegree);
  std::vector<int> num_insert_t;
  for(int ii=0; ii<num_knotspan_t; ++ii)
    num_insert_t.push_back(num_inserted_y);

  int num_knotspan_u = knotVec_check(uKnots, uDegree);
  std::vector<int> num_insert_u;
  for(int ii=0; ii<num_knotspan_u; ++ii)
    num_insert_u.push_back(num_inserted_z);

  std::vector<double> insertKnots_s, insertKnots_t, insertKnots_u;

  double hs_max, ht_max, hu_max, hs_min, ht_min, hu_min;
  hRefine_newKnot_Generator(sKnots, insertKnots_s, num_insert_s, hs_max, hs_min );
  hRefine_newKnot_Generator(tKnots, insertKnots_t, num_insert_t, ht_max, ht_min );
  hRefine_newKnot_Generator(uKnots, insertKnots_u, num_insert_u, hu_max, hu_min );


  kRefinement(addSDegree, addTDegree, addUDegree, insertKnots_s, insertKnots_t,
      insertKnots_u, sKnots, tKnots, uKnots, ctrlPts, cpDim, sDegree, tDegree, uDegree );

  cout<<endl<<"=== New knot vectors generated. \n";

  // Generate extraction operators
  std::vector<NURBS_T::BezierElem*> seg_x, seg_y, seg_z;
  std::vector<double> fake_ctrl, fake_bCtrlPts;
  NURBS_T::decomposeCurve(sKnots, sDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_x, false );
  NURBS_T::decomposeCurve(tKnots, tDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_y, false );
  NURBS_T::decomposeCurve(uKnots, uDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_z, false );

  cout<<endl<<"=== Extraction operators generated. \n";

  // Project down control points
  NURBS_T::projectDown(ctrlPts, spatialDim);

  // write out a control net for the volume
  if(isWriteCNet)
  {
    NURBS_T::writeControlNet( "volume_control_net", sKnots, tKnots, uKnots,
        sDegree, tDegree, uDegree, ctrlPts, spatialDim );
  }

  //Global Mesh Generation
  IMesh * Mesh;
  Mesh = new Mesh_NURBS_1Patch_3D( sDegree, tDegree, uDegree, 
      hs_max, ht_max, hu_max, hs_min, ht_min, hu_min, sKnots, tKnots, uKnots ); 

  Mesh->print_info();

  // IEN array
  IIEN * IEN;
  IEN = new IEN_NURBS_1Patch_3D(Mesh);

  // generate partition of element and node
  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, Mesh, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( Mesh, "epart", "npart" );
  else
  {
    cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<endl;
    exit(1);
  }

  // Node reorder due to partition
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, Mesh->get_nFunc());

  mnindex->write_hdf5("node_mapping");

  // Boundary conditions
  std::vector<BoundaryCond *> BC_list;

  BoundaryCond * bc_1 = new BoundaryCond(Mesh, 1);


  BC_list.push_back(bc_1); 

  if( int(BC_list.size()) != dofNum )
  {
    cerr<<"ERROR: The Boundary conditions does not match the number of unknowns. \n";
    exit(1);
  }

  cout<<"\n=== Start Partition ... \n";

  // Generate local informations
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_NURBS_1Patch_3D_METIS(Mesh, global_part, mnindex, IEN,
        ctrlPts, seg_x, seg_y, seg_z, proc_rank, proc_size, dofNum, elemType,
        isPrintPartInfo );

    part->write(part_file.c_str());    

    BC_Partition * bcpart = new BC_Partition(part, mnindex, BC_list);

    bcpart->write_hdf5(part_file.c_str());

    delete part;
    delete bcpart;
  }

  // Clean objects in memory
  cout<<endl<<"=== Clean memory. \n";
  std::vector<NURBS_T::BezierElem*>::iterator it_bezier;
  for( it_bezier = seg_x.begin(); it_bezier != seg_x.end(); ++it_bezier )
    delete *it_bezier;
  for( it_bezier = seg_y.begin(); it_bezier != seg_y.end(); ++it_bezier )
    delete *it_bezier;
  for( it_bezier = seg_z.begin(); it_bezier != seg_z.end(); ++it_bezier )
    delete *it_bezier;

  std::vector<BoundaryCond *>::iterator it_bc;
  for(it_bc = BC_list.begin(); it_bc != BC_list.end(); ++it_bc)
    delete *it_bc;


  delete IEN;
  delete Mesh;
  delete global_part;
  delete mnindex;
  PetscFinalize();
  return 0;
}
//EOF
