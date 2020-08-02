// ==================================================================
// unit_cube_structural_mesher.cpp
//
// This is a code that generate tetrahedron mesh based on a structural
// grid. The idea is to divide a unit cube by sub-dividing it into
// 20 tetrahedrons. The volume mesh and the boundary mesh in tetrahedron
// and triangle will be writen into .vtu and .vtp files separately.
// 
//             4 ____________________ 5
//            /|                     /|
//           / |                    / |
//          /  |                   /  |         u
//         /   |                  /   |         ^
//        6______________________7    |         |
//        |    |                |     |         |
//        |    |                |     |         |
//        |    |0 --------------|---- | 2       -------> t
//        |    /                |    /         /
//        |   /                 |   /         /
//        |  /                  |  /        |/_
//        | /                   | /         s
//        |/____________________|/
//        1                     3
//
// 
// Date Created: Jan 29 2017
// ==================================================================
#include "NURBS_FileIO.hpp"
#include "kRefinement.hpp"
#include "Mesh_NURBS_1Patch_3D.hpp"
#include "Mesh_NURBS_1Patch_3D_nze.hpp"
#include "IEN_NURBS_1Patch_3D_wPtr.hpp"
#include "Tet_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"

using std::cout;
using std::endl;

int get_node_pos(const IIEN * const &inien, 
    const std::vector<int> &ingnode, 
    const int &inelem, const int &locindex )
{
  int nn = inien -> get_IEN( inelem, locindex );
  int pos = VEC_T::get_pos(ingnode, nn);
  assert(pos >= 0);
  return pos;
}


void append_pts(const std::vector<double> &all,
    const int &loc, std::vector<double> &out )
{
  out.push_back( all[4*loc] );
  out.push_back( all[4*loc+1] );
  out.push_back( all[4*loc+2] );
}


int main( int argc, char * argv[] )
{
  std::vector<double> sKnots, tKnots, uKnots, ctrlPts;
  int sDegree, tDegree, uDegree, numCPts;
  const int cpDim = 4;
  const int spatialDim = 3;

  int addSDegree = 0, addTDegree = 0, addUDegree = 0;
  int num_inserted_x = 0, num_inserted_y = 0, num_inserted_z = 0;

  char * char_home_dir = getenv("HOME");
  std::string geo_file(char_home_dir);
  geo_file.append("/PERIGEE/input/geometry_3d_cube_0d1.txt");

  PetscMPIInt rank, size;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

  SYS_T::GetOptionInt("-elemx", num_inserted_x);
  SYS_T::GetOptionInt("-elemy", num_inserted_y);
  SYS_T::GetOptionInt("-elemz", num_inserted_z);
  SYS_T::GetOptionString("-geo_file", geo_file);

  cout<<"==== /Command Line Arguments ===="<<endl;
  cout<<" -elemx: "<<num_inserted_x<<endl;
  cout<<" -elemy: "<<num_inserted_y<<endl;
  cout<<" -elemz: "<<num_inserted_z<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;
  
  ifstream infile( geo_file.c_str(), ifstream::in );
  if(infile.is_open() == false)
  {
    std::cerr<<"ERROR: Can not find file: "<<geo_file<<endl;
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

  cout<<"\n===> New knot vectors generated. \n";

  NURBS_T::projectDown(ctrlPts, spatialDim);

  IMesh * mesh_wz = new Mesh_NURBS_1Patch_3D( sDegree, tDegree, uDegree,
      hs_max, ht_max, hu_max, hs_min, ht_min, hu_min, sKnots, tKnots, uKnots );

  mesh_wz->print_info();

  IIEN * IEN_wz = new IEN_NURBS_1Patch_3D_wPtr( mesh_wz );

  IIEN * IEN = new IEN_NURBS_1Patch_3D_wPtr( IEN_wz, mesh_wz );

  delete IEN_wz;

  IMesh * mesh = new Mesh_NURBS_1Patch_3D_nze(mesh_wz);

  delete mesh_wz;

  mesh->print_info();

  // Get useful hex mesh quantities
  const int nFunc = mesh -> get_nFunc();
  const int nFunc_x = mesh -> get_nFunc_x();
  const int nFunc_y = mesh -> get_nFunc_y();
  const int nFunc_z = mesh -> get_nFunc_z();
  
  const int nElem = mesh -> get_nElem();
  const int nElem_x = mesh -> get_nElem_x();
  const int nElem_y = mesh -> get_nElem_y();
  const int nElem_z = mesh -> get_nElem_z();
  
  const int nElem_tet = nElem * 6;

  std::vector<int> tetIEN; tetIEN.resize(nElem_tet * 4);
  std::vector<int> temp_IEN; temp_IEN.resize(8);

  // Decompose a cube into tets: IEN setting
  for(int ee=0; ee<nElem; ++ee)
  {
    for(int ii=0; ii<8; ++ii) temp_IEN[ii] = IEN->get_IEN(ee, ii);

    int tetee_offset = ee * 6;

    tetIEN[tetee_offset*4]   = temp_IEN[0];
    tetIEN[tetee_offset*4+1] = temp_IEN[1];
    tetIEN[tetee_offset*4+2] = temp_IEN[2];
    tetIEN[tetee_offset*4+3] = temp_IEN[6];

    tetee_offset += 1;
    tetIEN[tetee_offset*4]   = temp_IEN[0];
    tetIEN[tetee_offset*4+1] = temp_IEN[4];
    tetIEN[tetee_offset*4+2] = temp_IEN[1];
    tetIEN[tetee_offset*4+3] = temp_IEN[6];

    tetee_offset += 1;
    tetIEN[tetee_offset*4]   = temp_IEN[1];
    tetIEN[tetee_offset*4+1] = temp_IEN[4];
    tetIEN[tetee_offset*4+2] = temp_IEN[5];
    tetIEN[tetee_offset*4+3] = temp_IEN[6];

    tetee_offset += 1;
    tetIEN[tetee_offset*4]   = temp_IEN[1];
    tetIEN[tetee_offset*4+1] = temp_IEN[3];
    tetIEN[tetee_offset*4+2] = temp_IEN[2];
    tetIEN[tetee_offset*4+3] = temp_IEN[6];

    tetee_offset += 1;
    tetIEN[tetee_offset*4]   = temp_IEN[1];
    tetIEN[tetee_offset*4+1] = temp_IEN[3];
    tetIEN[tetee_offset*4+2] = temp_IEN[6];
    tetIEN[tetee_offset*4+3] = temp_IEN[7];

    tetee_offset += 1;
    tetIEN[tetee_offset*4]   = temp_IEN[1];
    tetIEN[tetee_offset*4+1] = temp_IEN[6];
    tetIEN[tetee_offset*4+2] = temp_IEN[5];
    tetIEN[tetee_offset*4+3] = temp_IEN[7];
  }

  std::vector<double> tetctrlPts;
  tetctrlPts.resize(3*nFunc);
  for(int ii=0; ii<nFunc; ++ii)
  {
    tetctrlPts[3*ii]   = ctrlPts[4*ii];
    tetctrlPts[3*ii+1] = ctrlPts[4*ii+1];
    tetctrlPts[3*ii+2] = ctrlPts[4*ii+2];
  }
  
  cout<<"Tet nElem: "<<nElem_tet<<endl;
  cout<<"Tet nFunc: "<<nFunc<<endl;

  IIEN * tIEN = new IEN_Tetra_P1(nElem_tet, tetIEN);

  cout<<"\n===> Checking the tet4 mesh data... "; 
  TET_T::Tet4 * teton = new TET_T::Tet4();

  teton -> reset( tetctrlPts, tIEN, 0 );
  double teton_max_vol = teton -> get_volume();
  double teton_min_vol = teton_max_vol; 
  for(int ee = 0; ee<nElem_tet; ++ee)
  {
    // read in the ee-th element
    teton->reset( tetctrlPts, tIEN, ee );

    // Calculate the min / max edge length and return the element index if the
    // ratio is greater than the tolerance 
    if( teton->get_aspect_ratio() > 3.0 ) 
      std::cout<<ee<<'\t'<<teton->get_aspect_ratio()<<std::endl;

    // Check if there are elements that has distorted numbering of nodes
    double teton_ee_vol = teton -> get_volume();
    if( teton_ee_vol < 0.0 ) 
      std::cout<<"Element "<<ee<<" is distorted! \n";

    if( teton_max_vol < teton_ee_vol) teton_max_vol = teton_ee_vol;
    if( teton_min_vol > teton_ee_vol) teton_min_vol = teton_ee_vol;
  }
  cout<<" done. \n";
  cout<<"- maximum tetrahedron volume : "<<teton_max_vol<<endl;
  cout<<"- minimum tetrahedron volume : "<<teton_min_vol<<endl;

  cout<<"===> Write vol mesh file: ";
  std::string vol_file("vol");
  clock_t timer_vol_mesh = clock();
  // Write the tetrahedral volume element info on disk in a .vtu file
  TET_T::write_tet_grid(vol_file, nFunc, nElem_tet, tetctrlPts, tetIEN );
  timer_vol_mesh = clock() - timer_vol_mesh;
  cout<<(double)timer_vol_mesh/(double)CLOCKS_PER_SEC<<" sec.\n";

  // Prepare for boundary surfaces
  std::vector<double> tripts;
  std::vector<int> triien, gnode, gelem;
  tripts.clear(); triien.clear(); gnode.clear(); gelem.clear();
  
  int tri_numpts = nFunc_x * nFunc_y;
  int tri_numcls = 2 * nElem_x * nElem_y;

  for(int fy=0; fy<nFunc_y; ++fy)
  {
    for(int fx=0; fx<nFunc_x; ++fx)
    {
      int func = fy * nFunc_x + fx;
      append_pts( ctrlPts, func, tripts );
      gnode.push_back(func);
    }
  }

  for(int ex=0; ex<nElem_x; ++ex)
  {
    for(int ey=0; ey<nElem_y; ++ey)
    {
      int elem = ey * nElem_x + ex;
      gelem.push_back(6*elem);

      triien.push_back( get_node_pos(IEN, gnode, elem, 0) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 1) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 2) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem );
      teton->get_face_id(IEN->get_IEN(elem,0), 
          IEN->get_IEN(elem,1), IEN->get_IEN(elem,2));
#endif

      gelem.push_back(6*elem+3);

      triien.push_back( get_node_pos(IEN, gnode, elem, 1) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 2) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 3) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+3 );
      teton->get_face_id(IEN->get_IEN(elem,1), 
          IEN->get_IEN(elem,2), IEN->get_IEN(elem,3));
#endif
    }
  }

  TET_T::write_triangle_grid( "sur.1", tri_numpts, tri_numcls,
      tripts, triien, gnode, gelem );

  // Top face
  tripts.clear(); triien.clear(); gnode.clear(); gelem.clear();
  for(int fy=0; fy<nFunc_y; ++fy)
  {
    for(int fx=0; fx<nFunc_x; ++fx)
    {
      int func = fy * nFunc_x + fx + (nFunc_z -1) * nFunc_x * nFunc_y;
      append_pts( ctrlPts, func, tripts );
      gnode.push_back(func);
    }
  }
  
  for(int ex=0; ex<nElem_x; ++ex)
  {
    for(int ey=0; ey<nElem_y; ++ey)
    {
      int elem = ey * nElem_x + ex + (nElem_z-1)*nElem_y * nElem_x;
      gelem.push_back(6*elem+2);

      triien.push_back( get_node_pos(IEN, gnode, elem, 4) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 5) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 6) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+2 );
      teton->get_face_id(IEN->get_IEN(elem,4), 
          IEN->get_IEN(elem,5), IEN->get_IEN(elem,6));
#endif

      gelem.push_back(6*elem+5);

      triien.push_back( get_node_pos(IEN, gnode, elem, 5) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 6) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 7) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+5 );
      teton->get_face_id(IEN->get_IEN(elem,5), 
          IEN->get_IEN(elem,6), IEN->get_IEN(elem,7));
#endif
    }
  }
  TET_T::write_triangle_grid( "sur.2", tri_numpts, tri_numcls,
      tripts, triien, gnode, gelem );

  // Back face
  tripts.clear(); triien.clear(); gnode.clear(); gelem.clear();
  tri_numpts = nFunc_y * nFunc_z;
  tri_numcls = 2 * nElem_y * nElem_z;

  for(int fy=0; fy<nFunc_y; ++fy)
  {
    for(int fz=0; fz<nFunc_z; ++fz)
    {
      int func = fz * nFunc_x * nFunc_y + fy * nFunc_x;
      append_pts( ctrlPts, func, tripts );
      gnode.push_back(func);
    }
  }

  for(int ez=0; ez<nElem_z; ++ez)
  {
    for(int ey=0; ey<nElem_y; ++ey)
    {
      int elem = ez * nElem_x * nElem_y + ey * nElem_x;

      gelem.push_back(6*elem+1);
      
      triien.push_back( get_node_pos(IEN, gnode, elem, 0) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 4) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 6) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+1 );
      teton->get_face_id(IEN->get_IEN(elem,0), 
          IEN->get_IEN(elem,4), IEN->get_IEN(elem,6));
#endif
      gelem.push_back(6*elem+0);
      
      triien.push_back( get_node_pos(IEN, gnode, elem, 0) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 2) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 6) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+0 );
      teton->get_face_id(IEN->get_IEN(elem,0), 
          IEN->get_IEN(elem,2), IEN->get_IEN(elem,6));
#endif
    }
  }
  TET_T::write_triangle_grid( "sur.3", tri_numpts, tri_numcls,
      tripts, triien, gnode, gelem );

  // Front face
  tripts.clear(); triien.clear(); gnode.clear(); gelem.clear();

  for(int fy=0; fy<nFunc_y; ++fy)
  {
    for(int fz=0; fz<nFunc_z; ++fz)
    {
      int func = fz * nFunc_x * nFunc_y + fy * nFunc_x + nFunc_x - 1;
      append_pts( ctrlPts, func, tripts );
      gnode.push_back(func);
    }
  }

  for(int ez=0; ez<nElem_z; ++ez)
  {
    for(int ey=0; ey<nElem_y; ++ey)
    {
      int elem = ez * nElem_x * nElem_y + ey * nElem_x + nElem_x - 1;

      gelem.push_back(6*elem+4);
      
      triien.push_back( get_node_pos(IEN, gnode, elem, 1) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 3) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 7) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+4 );
      teton->get_face_id(IEN->get_IEN(elem,1), 
          IEN->get_IEN(elem,3), IEN->get_IEN(elem,7));
#endif
      gelem.push_back(6*elem+5);
      
      triien.push_back( get_node_pos(IEN, gnode, elem, 1) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 5) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 7) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+5 );
      teton->get_face_id(IEN->get_IEN(elem,1), 
          IEN->get_IEN(elem,5), IEN->get_IEN(elem,7));
#endif
    }
  }
  TET_T::write_triangle_grid( "sur.4", tri_numpts, tri_numcls,
      tripts, triien, gnode, gelem );

  // Left face
  tripts.clear(); triien.clear(); gnode.clear(); gelem.clear();
  tri_numpts = nFunc_x * nFunc_z;
  tri_numcls = 2 * nElem_x * nElem_z;

  for(int fx=0; fx<nFunc_x; ++fx)
  {
    for(int fz=0; fz<nFunc_z; ++fz)
    {
      int func = fz * nFunc_x * nFunc_y + fx;
      append_pts( ctrlPts, func, tripts );
      gnode.push_back(func);
    }
  }

  for(int ex=0; ex<nElem_x; ++ex)
  {
    for(int ez=0; ez<nElem_z; ++ez)
    {
      int elem = ez * nElem_x * nElem_y + ex;

      gelem.push_back(6*elem+1);

      triien.push_back( get_node_pos(IEN, gnode, elem, 0) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 1) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 4) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+1 );
      teton->get_face_id(IEN->get_IEN(elem,0),
          IEN->get_IEN(elem,1), IEN->get_IEN(elem,4));
#endif
      gelem.push_back(6*elem+2);

      triien.push_back( get_node_pos(IEN, gnode, elem, 1) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 4) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 5) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+2 );
      teton->get_face_id(IEN->get_IEN(elem,1),
          IEN->get_IEN(elem,4), IEN->get_IEN(elem,5));
#endif
    }
  }
  TET_T::write_triangle_grid( "sur.5", tri_numpts, tri_numcls,
      tripts, triien, gnode, gelem );

  // Right face
  tripts.clear(); triien.clear(); gnode.clear(); gelem.clear();

  for(int fx=0; fx<nFunc_x; ++fx)
  {
    for(int fz=0; fz<nFunc_z; ++fz)
    {
      int func = fz * nFunc_x * nFunc_y + (nFunc_y-1)*nFunc_x + fx;
      append_pts( ctrlPts, func, tripts );
      gnode.push_back(func);
    }
  }

  for(int ex=0; ex<nElem_x; ++ex)
  {
    for(int ez=0; ez<nElem_z; ++ez)
    {
      int elem = ez * nElem_x * nElem_y + ex + (nElem_y-1)*nElem_x;

      gelem.push_back(6*elem+3);

      triien.push_back( get_node_pos(IEN, gnode, elem, 2) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 3) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 6) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+3 );
      teton->get_face_id(IEN->get_IEN(elem,2),
          IEN->get_IEN(elem,3), IEN->get_IEN(elem,6));
#endif
      gelem.push_back(6*elem+4);

      triien.push_back( get_node_pos(IEN, gnode, elem, 3) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 6) );
      triien.push_back( get_node_pos(IEN, gnode, elem, 7) );

#ifdef ENABLE_TEST
      teton->reset( tetctrlPts, tIEN, 6*elem+4 );
      teton->get_face_id(IEN->get_IEN(elem,3),
          IEN->get_IEN(elem,6), IEN->get_IEN(elem,7));
#endif
    }
  }
  TET_T::write_triangle_grid( "sur.6", tri_numpts, tri_numcls,
      tripts, triien, gnode, gelem );


  delete teton; delete tIEN; delete mesh; delete IEN;
  PetscFinalize();
  return EXIT_SUCCESS;
}


// EOF
