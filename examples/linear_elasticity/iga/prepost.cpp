// ==================================================================
// prepost.cpp
// 
// This is the partitioning routine for all parallel postprocessors.
//
// Date: July 8 2016
// ==================================================================
#include "kRefinement.hpp"
#include "Mesh_NURBS_1Patch_3D.hpp"
#include "Mesh_NURBS_1Patch_3D_nze.hpp"
#include "IEN_NURBS_1Patch_3D_wPtr.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_NURBS_1Patch_3D_METIS.hpp"
#include "HDF5_PartReader.hpp"

using std::endl;
using std::cerr;

int main( int argc, char * argv[] )
{
  std::vector<double> sKnots, tKnots, uKnots, ctrlPts;
  int sDegree, tDegree, uDegree, numCPts;
  const int cpDim = 4;
  const int spatialDim = 3;

  const int probDim = 3;
  const int dofNum = 3; // displacement in x-y-z directions
  const int elemType = 0;

  int addSDegree = 1, addTDegree = 1, addUDegree = 1;
  int num_inserted_x = 8, num_inserted_y = 8, num_inserted_z = 8;

  char * char_home_dir = getenv("HOME");
  std::string geo_file(char_home_dir);
  geo_file.append("/PERIGEE/input/geometry_3d_cube.txt");

  std::string part_file("postpart");

  int cpu_size = 1;
  int in_ncommon = 9;
  bool isDualGraph = true;
  bool isread_part = true;

  PetscMPIInt rank, size;

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

  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-part_file", part_file);

  SYS_T::GetOptionBool("-isread_part", isread_part);

  std::cout<<"==== /Command Line Arguments ===="<<endl;
  std::cout<<" -elemx: "<<num_inserted_x<<endl;
  std::cout<<" -elemy: "<<num_inserted_y<<endl;
  std::cout<<" -elemz: "<<num_inserted_z<<endl;
  std::cout<<" -addSDegree: "<<addSDegree<<endl;
  std::cout<<" -addTDegree: "<<addTDegree<<endl;
  std::cout<<" -addUDegree: "<<addUDegree<<endl;
  std::cout<<" -geo_file: "<<geo_file<<endl;
  std::cout<<" -part_file: "<<part_file<<endl;
  std::cout<<" -cpu_size: "<<cpu_size<<endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) std::cout<<" -METIS_isDualGraph: true \n";
  else std::cout<<" -METIS_isDualGraph: false \n";
  if(isread_part) std::cout<<" -isread_part: true \n";
  else std::cout<<" -isread_part: false \n";
  std::cout<<"----------------------------------\n";
  std::cout<<"probDim: "<<probDim<<endl;
  std::cout<<"dofNum: "<<dofNum<<endl;
  std::cout<<"elemType: "<<elemType<<endl;
  std::cout<<"====  Command Line Arguments/ ===="<<endl;

  std::ifstream infile( geo_file.c_str(), std::ifstream::in );
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
  std::vector<int> num_insert_s;
  for(int ii=0; ii<num_knotspan_s; ++ii) num_insert_s.push_back(num_inserted_x);

  int num_knotspan_t = knotVec_check(tKnots, tDegree);
  std::vector<int> num_insert_t;
  for(int ii=0; ii<num_knotspan_t; ++ii) num_insert_t.push_back(num_inserted_y);

  int num_knotspan_u = knotVec_check(uKnots, uDegree);
  std::vector<int> num_insert_u;
  for(int ii=0; ii<num_knotspan_u; ++ii) num_insert_u.push_back(num_inserted_z);

  std::vector<double> insertKnots_s, insertKnots_t, insertKnots_u;

  double hs_max, ht_max, hu_max, hs_min, ht_min, hu_min;
  hRefine_newKnot_Generator(sKnots, insertKnots_s, num_insert_s, hs_max, hs_min );
  hRefine_newKnot_Generator(tKnots, insertKnots_t, num_insert_t, ht_max, ht_min );
  hRefine_newKnot_Generator(uKnots, insertKnots_u, num_insert_u, hu_max, hu_min );

  kRefinement(addSDegree, addTDegree, addUDegree, insertKnots_s, insertKnots_t,
      insertKnots_u, sKnots, tKnots, uKnots, ctrlPts, cpDim, sDegree, tDegree, uDegree );

  std::cout<<"\n===> New knot vectors generated. \n";

  std::vector<NURBS_T::BezierElem*> seg_x, seg_y, seg_z;
  std::vector<double> fake_ctrl, fake_bCtrlPts;
  NURBS_T::decomposeCurve(sKnots, sDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_x, false );
  NURBS_T::decomposeCurve(tKnots, tDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_y, false );
  NURBS_T::decomposeCurve(uKnots, uDegree, fake_ctrl, cpDim,
      fake_bCtrlPts, seg_z, false );

  std::cout<<"\n===> Extraction operators generated. \n";

  NURBS_T::projectDown(ctrlPts, spatialDim);

  IMesh * Mesh_wz = new Mesh_NURBS_1Patch_3D( sDegree, tDegree, uDegree,
      hs_max, ht_max, hu_max, hs_min, ht_min, hu_min, sKnots, tKnots, uKnots );
  IIEN * IEN_wz = new IEN_NURBS_1Patch_3D_wPtr( Mesh_wz );

  IIEN * IEN = new IEN_NURBS_1Patch_3D_wPtr( IEN_wz, Mesh_wz );

  delete IEN_wz;

  IMesh * Mesh = new Mesh_NURBS_1Patch_3D_nze( Mesh_wz );

  delete Mesh_wz;

  if(isread_part)
  {
    std::cout<<"\n=== Check the compatibility ...\n";
    std::cout<<"-- Read from part_p00000.h5 \n";
    std::string partfilebasename = "part";
    HDF5_PartReader * check_reader = new HDF5_PartReader(partfilebasename, 0);

    int check_sd, check_td, check_ud;
    check_reader->get_GMI_degree(check_sd, check_td, check_ud);
    assert(check_sd == sDegree);
    assert(check_td == tDegree);
    assert(check_ud == uDegree);

    std::cout<<"-- Degrees are compatible... \n";

    int check_ne, check_nex, check_ney, check_nez;
    check_reader->get_GMI_nElem(check_ne, check_nex, check_ney, check_nez);
    assert(check_ne == Mesh->get_nElem());
    assert(check_nex == Mesh->get_nElem_x());
    assert(check_ney == Mesh->get_nElem_y());
    assert(check_nez == Mesh->get_nElem_z());

    std::cout<<"-- nElem, nElem_x/y/z are compatible. "<<endl;
    delete check_reader;
  }

  Mesh->print_info();

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

  std::cout<<"\n=== Start Partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_NURBS_1Patch_3D_METIS(Mesh, global_part, mnindex, IEN,
        ctrlPts, seg_x, seg_y, seg_z, proc_rank, proc_size, dofNum, elemType,
        isPrintPartInfo );
    part->write(part_file.c_str());
    delete part;
  }

  std::cout<<endl<<"=== Clean memory. \n";
  std::vector<NURBS_T::BezierElem*>::iterator it_bezier;
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
