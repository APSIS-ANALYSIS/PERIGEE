// ==================================================================
// preprocess2D_main.cpp
// ------------------------------------------------------------------
// Objective:
// This is the preprocess code for mesh generation and partition for 
// 2d isogeometric/finite element analysis.
//
// Output:
// 1. Partitioned mesh information for each analysis processor;
// 2. (optional) vtk file for the input geometry
//
// Note:
// 1. This is a serial code.
// 2. The user may have to modify the code for different problems, 
//    including: 
//         the number of degree-of-freedom(dof);
//         each dof's boundary condition.
//
// -----------------------------------------------------------------
// Date: Dec 30 2013
// ==================================================================
#include "NURBS_FileIO.hpp"
#include "kRefinement.hpp"
#include "Mesh_NURBS_1Patch_2D.hpp"
#include "IEN_NURBS_1Patch_2D.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_NURBS_1Patch_2D_METIS.hpp"
#include "BC_Partition2D.hpp"

using namespace std;

int main(int argc, char * argv[])
{
  vector<double> sKnots, tKnots, ctrlPts;
  int sDegree, tDegree, numCPts;
  const int cpDim = 4;
  const int spatialDim = 3;
  const int probDim = 2;

  // The dof of the finite element problem
  const int dofNum = 1;
  
  // The types of finite element
  int elemType = 0;

  // Degrees to be added in s t directions
  int addSDegree = 1, addTDegree = 1;

  // Number of knots to be inserted in spatial discretizations
  int num_inserted_x = 8, num_inserted_y = 8;

  // Geometry file
  char * char_home_dir = getenv("HOME");
  string geo_file(char_home_dir);
  geo_file.append("/IsoPETSc3D/trunk/input/geometry_2d_square.txt");

  // partition file name
  string part_file("part");

  // partition parameters
  int cpu_size = 1;
  int in_ncommon = 9;
  bool isDualGraph = true;
  bool isWriteCNet = true;

  PetscMPIInt rank, size;
  
  // ====== PETSc Initialization ======
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
 
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

  SYS_T::GetOptionInt("-elemx", num_inserted_x);
  SYS_T::GetOptionInt("-elemy", num_inserted_y);
  SYS_T::GetOptionInt("-addSDegree", addSDegree);
  SYS_T::GetOptionInt("-addTDegree", addTDegree);
  
  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);
  SYS_T::GetOptionBool("-isWriteControlNet", isWriteCNet);

  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-part_file", part_file);

  // Print command line arguments
  cout<<"==== /Command Line Arguments ===="<<endl;
  cout<<" -elemx: "<<num_inserted_x<<endl;
  cout<<" -elemy: "<<num_inserted_y<<endl;
  cout<<" -addSDegree: "<<addSDegree<<endl;
  cout<<" -addTDegree: "<<addTDegree<<endl;
  cout<<" -dof: "<<dofNum<<endl;
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
  cout<<"elemType: "<<elemType<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // ====== Read Geometry Files ======
  int uDegree;
  vector<double> uKnots;
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

  if(uDegree != 0 || uKnots.size() > 1)
  {
    cerr<<"ERROR: The inpute geometry file is not 2D. \n";
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  NURBS_T::projectUp(ctrlPts, spatialDim);

  // k-refinement
  int num_knotspan_s = knotVec_check(sKnots, sDegree);
  vector<int> num_insert_s;
  for(int ii=0; ii<num_knotspan_s; ++ii)
    num_insert_s.push_back(num_inserted_x);

  int num_knotspan_t = knotVec_check(tKnots, tDegree);
  vector<int> num_insert_t;
  for(int ii=0; ii<num_knotspan_t; ++ii)
    num_insert_t.push_back(num_inserted_y);

  // generate new knot
  vector<double> insertKnots_s, insertKnots_t;
  double hs_max, ht_max, hs_min, ht_min;
  hRefine_newKnot_Generator(sKnots, insertKnots_s, num_insert_s, hs_max, hs_min );
  hRefine_newKnot_Generator(tKnots, insertKnots_t, num_insert_t, ht_max, ht_min );

  // k-refinement for the mesh
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

  // project down the control points
  NURBS_T::projectDown( ctrlPts, spatialDim );

  // write control net for the problem
  if( isWriteCNet )
  {
    NURBS_T::writeControlNet( "surface_control_net", sKnots, tKnots, uKnots,
        sDegree, tDegree, uDegree, ctrlPts, spatialDim );
    cout<<endl<<"=== surface_control_net.vtk file written on disk. \n";
  }

  // Global Mesh Info generator
  IMesh * Mesh;
  Mesh = new Mesh_NURBS_1Patch_2D( sDegree, tDegree,
      hs_max, ht_max, hs_min, ht_min, sKnots, tKnots );

  Mesh->print_mesh_info();

  // IEN array
  IIEN * IEN;
  IEN = new IEN_NURBS_1Patch_2D(Mesh);

  // Generate partition of element and node
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

  // Node reordering due to partition
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, Mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");


  // Boundary conditions
  BoundaryCond2D * bc_1 = new BoundaryCond2D(Mesh, 1);

  vector<BoundaryCond2D *> BC_list;

  BC_list.push_back(bc_1);

  if( int(BC_list.size()) != dofNum )
  {
    cerr<<"ERROR: The number of imposed BC's does not match the dofNum parameter.";
    cerr<<"\n       BC_list.size() =  "<<BC_list.size();
    cerr<<"\t       dofNum =  "<<dofNum<<endl;
    exit(1);
  }


  // ------------------------------------------------------
  // Partition the 2D mesh
  // ------------------------------------------------------
  cout<<"\n=== Start mesh partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_NURBS_1Patch_2D_METIS(Mesh, global_part, mnindex, IEN,
        ctrlPts, seg_x, seg_y, proc_rank, proc_size, dofNum, elemType,
        isPrintPartInfo );

    part->write(part_file.c_str());

    BC_Partition2D * bcpart = new BC_Partition2D(part, mnindex, BC_list);

    bcpart->write_hdf5(part_file.c_str());

    delete part;
    delete bcpart;
  }

  // Free memory allocations
  cout<<endl<<"=== Clean memory. \n";

  vector<NURBS_T::BezierElem*>::iterator it_bezier;
  for( it_bezier = seg_x.begin(); it_bezier != seg_x.end(); ++it_bezier )
    delete *it_bezier;
  for( it_bezier = seg_y.begin(); it_bezier != seg_y.end(); ++it_bezier )
    delete *it_bezier;

  vector<BoundaryCond2D *>::iterator it_bc;
  for(it_bc = BC_list.begin(); it_bc != BC_list.end(); ++it_bc)
    delete *it_bc;

  delete Mesh;
  delete IEN;
  delete global_part;
  delete mnindex;

  PetscFinalize();
  return 0;
}

// EOF
