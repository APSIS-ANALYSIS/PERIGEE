// ==================================================================
// preprocess_tets.cpp
// ------------------------------------------------------------------
// This preprocess code is used for handling the 3D geometry described
// by tetrahedral elements. 
//
// Date: Dec. 18 2016
// ==================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet4.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Tet.hpp"
#include "NodalBC_3D_vtp.hpp"
#include "ElemBC_3D_tet4.hpp"
#include "NBC_Partition_3D.hpp"
#include "EBC_Partition_vtp.hpp"

int main( int argc, char * argv[] )
{
  const int probDim = 3;
  const int dofNum = 3;
  const int elemType = 501;

  char * char_home_dir = getenv("HOME");
  std::string geo_file(char_home_dir);
  geo_file.append("/PERIGEE/input/cube_tet_jl/coarse/vol.vtu");

  std::string sur_file(char_home_dir);
  sur_file.append("/PERIGEE/input/cube_tet_jl/coarse/sur.");

  std::string part_file("part");

  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;

  PetscMPIInt rank, size;

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(size != 1) SYS_T::print_fatal("ERROR: preprocessor is a serial program! \n");

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);
  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-sur_file", sur_file);

  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -geo_file: "<<geo_file<<std::endl;
  std::cout<<" -sur_file: "<<sur_file<<std::endl;
  std::cout<<" -part_file: "<<part_file<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  if(isDualGraph) std::cout<<" -isDualGraph: true \n";
  else std::cout<<" -isDualGraph: false \n";
  std::cout<<"----------------------------------\n";
  std::cout<<" probDim: "<<probDim<<std::endl;
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<"====  Command Line Arguments/ ===="<<std::endl;

  // ----- Write the input argument into a HDF5 file
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", 
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);
 
  cmdh5w->write_intScalar("cpu_size", cpu_size); 
  cmdh5w->write_intScalar("in_ncommon", in_ncommon); 
  cmdh5w->write_intScalar("dofNum", dofNum); 
  cmdh5w->write_intScalar("probDim", probDim); 
  cmdh5w->write_intScalar("elemType", elemType); 
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("sur_file", sur_file);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; 
  H5Fclose(cmd_file_id);

  // Read file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN);

  if(int(vecIEN.size()) != nElem * 4) SYS_T::print_fatal("Error: the IEN from geo_file does not match the given number of element. \n");

  if(int(ctrlPts.size()) != nFunc * 3) SYS_T::print_fatal("Error: the ctrlPts from geo_file does not match the given number of nodes. \n");

  std::cout<<"nElem: "<<nElem<<std::endl;
  std::cout<<"nFunc: "<<nFunc<<std::endl;

  IIEN * IEN = new IEN_Tetra_P1(nElem, vecIEN);

  VEC_T::clean( vecIEN );

  // ----------------------------------------------------------------
  // check the tet mesh
  cout<<"\n===> Checking the tet4 mesh data... "; 
  TET_T::Tet4 * teton = new TET_T::Tet4();

  teton -> reset( ctrlPts, IEN, 0 );
  double teton_max_vol = teton -> get_volume();
  double teton_min_vol = teton_max_vol; 
  double teton_max_h = teton -> get_diameter();
  double teton_min_h = teton_max_h;
  for(int ee = 0; ee<nElem; ++ee)
  {
    // read in the ee-th element
    teton->reset( ctrlPts, IEN, ee );

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
    
    double teton_he = teton -> get_diameter();

    if( teton_max_h < teton_he ) teton_max_h = teton_he;

    if( teton_min_h > teton_he ) teton_min_h = teton_he;
  
  }
  delete teton;
  cout<<" done. \n";
  cout<<"- maximum tetrahedron volume : "<<teton_max_vol<<endl;
  cout<<"- minimum tetrahedron volume : "<<teton_min_vol<<endl;
  cout<<"- maximum tetrahedron diameter : "<<teton_max_h<<endl;
  cout<<"- minimum tetrahedron diameter : "<<teton_min_h<<endl;
  // ----------------------------------------------------------------

  IMesh * mesh = new Mesh_Tet4(nFunc, nElem);
  mesh -> print_mesh_info();
  
  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "epart", "npart" );
  else
  {
    cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<endl;
    exit(EXIT_FAILURE);
  }

  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");

  // ----------------------------------------------------------------
  // Setup boundary condition
  std::vector<INodalBC *> NBC_list;
  NBC_list.clear();
  NBC_list.resize( dofNum );

  std::string file4nbc1 = sur_file + SYS_T::to_string(1)+ ".vtp";
  std::string file4nbc2 = sur_file + SYS_T::to_string(2)+ ".vtp";
  std::string file4nbc3 = sur_file + SYS_T::to_string(3)+ ".vtp";
  std::string file4nbc4 = sur_file + SYS_T::to_string(4)+ ".vtp";
  std::string file4nbc5 = sur_file + SYS_T::to_string(5)+ ".vtp";
  std::string file4nbc6 = sur_file + SYS_T::to_string(6)+ ".vtp";
  
  /*
  std::vector<std::string> n0list, n1list, n2list;
  n0list.push_back( file4nbc1 );

  n1list.push_back( file4nbc2 );
  n1list.push_back( file4nbc3 );

  n2list.push_back( file4nbc2 );
  n2list.push_back( file4nbc5 );
  */
  
  std::vector<std::string> nlist;
  nlist.push_back( file4nbc1 );
  nlist.push_back( file4nbc2 );

  NBC_list[0] = new NodalBC_3D_vtp( nlist, nFunc );
  NBC_list[1] = new NodalBC_3D_vtp( nlist, nFunc );
  NBC_list[2] = new NodalBC_3D_vtp( nlist, nFunc );
  
  std::vector<std::string> ebclist;
  //ebclist.push_back( file4nbc1 );
  //ebclist.push_back( file4nbc2 );
  ebclist.push_back( file4nbc5 );
  ebclist.push_back( file4nbc6 );
  ebclist.push_back( file4nbc4 );
  ebclist.push_back( file4nbc3 );

  ElemBC * ebc = new ElemBC_3D_tet4( ebclist );
  
  // Correct the triangle's IEN so that we can get the outward normal 
  // direction easily. See the document for this reset function
  ebc -> resetTriIEN_outwardnormal( IEN );
  // ----------------------------------------------------------------

  const bool isPrintPartInfo = true;
  const int proc_size = cpu_size;
  
  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0;

  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    IPart * part = new Part_Tet( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, proc_size, dofNum, elemType,
        isPrintPartInfo );

    part -> write( part_file.c_str() );

    part -> print_part_loadbalance_edgecut();

    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
    nbcpart -> write_hdf5(part_file.c_str());

    IEBC_Partition * ebcpart = new EBC_Partition_vtp(part, mnindex, ebc);

    ebcpart -> write_hdf5(part_file.c_str());

    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete ebcpart;
  }

  VEC_T::write_int_h5("NumLocalNode","nln", list_nlocalnode);

  cout<<endl<<"----- Partition Quality: "<<endl;
  cout<<"The largest ghost / local node ratio is: ";
  cout<<*std::max_element(&list_ratio_g2l[0], &list_ratio_g2l[cpu_size-1])<<endl;

  cout<<"The smallest ghost / local node ratio is: ";
  cout<<*std::min_element(&list_ratio_g2l[0], &list_ratio_g2l[cpu_size-1])<<endl;

  cout<<"The summation of the number of ghost nodes is: "<<sum_nghostnode<<endl;

  cout<<"The maximum badnode number is: ";
  cout<<*std::max_element(&list_nbadnode[0], &list_nbadnode[cpu_size-1])<<endl;

  const int maxpart_nlocalnode = *std::max_element(&list_nlocalnode[0],
      &list_nlocalnode[cpu_size-1]);
  const int minpart_nlocalnode = *std::min_element(&list_nlocalnode[0],
      &list_nlocalnode[cpu_size-1]);

  cout<<"The maximum and minimum local node numbers are ";
  cout<<maxpart_nlocalnode<<"\t";
  cout<<minpart_nlocalnode<<endl;
  cout<<"The maximum / minimum of local node is: ";
  cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<endl;

  // Free memory
  delete ebc;
  std::vector<INodalBC *>::iterator it_nbc;
  for(it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;
  
  delete mnindex; delete global_part; delete mesh; delete IEN; 
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
