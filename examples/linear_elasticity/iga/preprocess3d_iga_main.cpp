// ============================================================================
// preprocess3d_main.cpp
// ----------------------------------------------------------------------------
// This preprocess code is used for mesh generation and partitioning for 3D
// hyper-elastic solids. The mesh is 3-dimensional B-splines or NURBS.
//
// Date: June 21 2016
// ============================================================================
#include "NURBS_FileIO.hpp"
#include "kRefinement.hpp"
#include "Mesh_NURBS_1Patch_3D.hpp"
#include "Mesh_NURBS_1Patch_3D_nze.hpp"
#include "IEN_NURBS_1Patch_3D_wPtr.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Map_Node_Index.hpp"
#include "NodalBC_3D_Linearelastic.hpp"
#include "ElemBC_3D_Linearelastic.hpp"
#include "Part_NURBS_1Patch_3D_METIS.hpp"
#include "NBC_Partition_3D.hpp"
#include "EBC_Partition_3D.hpp"

using namespace std;

int main(int argc, char * argv[])
{
  vector<double> sKnots, tKnots, uKnots, ctrlPts;
  int sDegree, tDegree, uDegree, numCPts;
  const int cpDim = 4;
  const int spatialDim = 3;
  
  const int probDim = 3;
  const int dofNum = 3; // displacement in x-y-z directions
  const int elemType = 0;

  int addSDegree = 1, addTDegree = 1, addUDegree = 1;
  int num_inserted_x = 8, num_inserted_y = 8, num_inserted_z = 8;

  char * char_home_dir = getenv("HOME");
  string geo_file(char_home_dir);
  geo_file.append("/PERIGEE/input/geometry_3d_cube.txt");

  string part_file("part");

  int cpu_size = 1;
  int in_ncommon = 9;
  bool isDualGraph = true;
  bool isWriteCNet = true;

  PetscMPIInt rank, size;

  // =========== PETSc Initialize ============
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

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
  if(isDualGraph) cout<<" -isDualGraph: true \n";
  else cout<<" -isDualGraph: false \n";
  if(isWriteCNet) cout<<" -isWriteControlNet: true \n";
  else cout<<" -isWriteControlNet: false \n";
  cout<<"----------------------------------\n";
  cout<<"probDim: "<<probDim<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // Read Geometry
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
  vector<int> num_insert_s;
  for(int ii=0; ii<num_knotspan_s; ++ii) num_insert_s.push_back(num_inserted_x);

  int num_knotspan_t = knotVec_check(tKnots, tDegree);
  vector<int> num_insert_t;
  for(int ii=0; ii<num_knotspan_t; ++ii) num_insert_t.push_back(num_inserted_y);

  int num_knotspan_u = knotVec_check(uKnots, uDegree);
  vector<int> num_insert_u;
  for(int ii=0; ii<num_knotspan_u; ++ii) num_insert_u.push_back(num_inserted_z);

  vector<double> insertKnots_s, insertKnots_t, insertKnots_u;

  double hs_max, ht_max, hu_max, hs_min, ht_min, hu_min;
  hRefine_newKnot_Generator(sKnots, insertKnots_s, num_insert_s, hs_max, hs_min );
  hRefine_newKnot_Generator(tKnots, insertKnots_t, num_insert_t, ht_max, ht_min );
  hRefine_newKnot_Generator(uKnots, insertKnots_u, num_insert_u, hu_max, hu_min );

  kRefinement(addSDegree, addTDegree, addUDegree, insertKnots_s, insertKnots_t,
      insertKnots_u, sKnots, tKnots, uKnots, ctrlPts, cpDim, sDegree, tDegree, uDegree );

  cout<<"\n===> New knot vectors generated. \n";

  // Extraction operator
  vector<NURBS_T::BezierElem*> seg_x, seg_y, seg_z;
  vector<double> fake_ctrl, fake_bCtrlPts;
  NURBS_T::decomposeCurve(sKnots, sDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_x, false );
  NURBS_T::decomposeCurve(tKnots, tDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_y, false );
  NURBS_T::decomposeCurve(uKnots, uDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_z, false );

  cout<<"\n===> Extraction operators generated. \n";

  NURBS_T::projectDown(ctrlPts, spatialDim);

  if(isWriteCNet)
  {
    NURBS_T::writeControlNet( "volume_control_net", sKnots, tKnots, uKnots,
        sDegree, tDegree, uDegree, ctrlPts, spatialDim );

    cout<<"\n===> Control net is written. \n";
  }
  
  IMesh * mesh_wz = new Mesh_NURBS_1Patch_3D( sDegree, tDegree, uDegree,
      hs_max, ht_max, hu_max, hs_min, ht_min, hu_min, sKnots, tKnots, uKnots );

  mesh_wz->print_mesh_info();

  IIEN * IEN_wz = new IEN_NURBS_1Patch_3D_wPtr( mesh_wz );

  IIEN * IEN = new IEN_NURBS_1Patch_3D_wPtr( IEN_wz, mesh_wz );

  delete IEN_wz;

  IMesh * mesh = new Mesh_NURBS_1Patch_3D_nze(mesh_wz);
  
  delete mesh_wz;

  mesh->print_mesh_info();

  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "epart", "npart" );
  else
  {
    cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<endl;
    exit(1);
  }

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");


  // Nodal BC specification
  vector<INodalBC *> NBC_list;
  NBC_list.clear();
  NBC_list.resize(dofNum);
  NBC_list[0] = new NodalBC_3D_Linearelastic(mesh->get_nFunc(), 
      mesh->get_nFunc_x(), mesh->get_nFunc_y(), mesh->get_nFunc_z(), 7);

  NBC_list[1] = new NodalBC_3D_Linearelastic(mesh->get_nFunc(), 
      mesh->get_nFunc_x(), mesh->get_nFunc_y(), mesh->get_nFunc_z(), 7);

  NBC_list[2] = new NodalBC_3D_Linearelastic(mesh->get_nFunc(), 
      mesh->get_nFunc_x(), mesh->get_nFunc_y(), mesh->get_nFunc_z(), 7);

  // Elemental BC specification 
  IElemBC * ebc = new ElemBC_3D_Linearelastic(mesh->get_nElem_x(), 
      mesh->get_nElem_y(), mesh->get_nElem_z(), 2);

  // Partition the mesh
  cout<<"\n=== Start mesh partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;

  vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  vector<double> list_ratio_g2l;

  int sum_nghostnode = 0;

  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_NURBS_1Patch_3D_METIS( mesh, global_part, mnindex, IEN,
        ctrlPts, seg_x, seg_y, seg_z, proc_rank, proc_size, dofNum, elemType,
        isPrintPartInfo );

    part->write(part_file.c_str());

    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);

    nbcpart->write_hdf5(part_file.c_str());

    IEBC_Partition * ebcpart = new EBC_Partition_3D(part, ebc);

    ebcpart->write_hdf5(part_file.c_str());

    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();

    delete part;
    delete ebcpart;
    delete nbcpart;
  }

  VEC_T::write_int_h5("NumLocalNode","nln", list_nlocalnode);

  cout<<endl<<"----- Partition Quality: "<<endl;
  cout<<"The largest ghost / local node ratio is: ";
  cout<<*max_element(&list_ratio_g2l[0], &list_ratio_g2l[cpu_size-1])<<endl;

  cout<<"The smallest ghost / local node ratio is: ";
  cout<<*min_element(&list_ratio_g2l[0], &list_ratio_g2l[cpu_size-1])<<endl;

  cout<<"The summation of the number of ghost nodes is: "<<sum_nghostnode<<endl;

  cout<<"The maximum badnode number is: ";
  cout<<*max_element(&list_nbadnode[0], &list_nbadnode[cpu_size-1])<<endl;

  const int maxpart_nlocalnode = *max_element(&list_nlocalnode[0],
      &list_nlocalnode[cpu_size-1]);
  const int minpart_nlocalnode = *min_element(&list_nlocalnode[0],
      &list_nlocalnode[cpu_size-1]);

  cout<<"The maximum and minimum local node numbers are ";
  cout<<maxpart_nlocalnode<<"\t";
  cout<<minpart_nlocalnode<<endl;
  cout<<"The maximum / minimum of local node is: ";
  cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<endl;


  // Free memory allocations
  vector<NURBS_T::BezierElem*>::iterator it_bezier;
  for( it_bezier = seg_x.begin(); it_bezier != seg_x.end(); ++it_bezier )
    delete *it_bezier;
  for( it_bezier = seg_y.begin(); it_bezier != seg_y.end(); ++it_bezier )
    delete *it_bezier;
  for( it_bezier = seg_z.begin(); it_bezier != seg_z.end(); ++it_bezier )
    delete *it_bezier;

  vector<INodalBC *>::iterator it_nbc;
  for(it_nbc= NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;

  delete ebc; 
  delete mnindex; delete global_part; delete IEN; delete mesh;
  PetscFinalize();
  return 0;
}

// EOF
