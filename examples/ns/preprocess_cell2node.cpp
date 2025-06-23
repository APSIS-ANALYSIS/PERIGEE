#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "NodalBC.hpp"
#include "NBC_Partition.hpp"
#include "VTK_Tools.hpp"
#include "yaml-cpp/yaml.h"

void write_ansys_node_data(const IPart * const &part,
  const Map_Node_Index * const &mnindex,
  const std::vector<double> &pts, const std::string &FileName,
  const std::string &GroupName, const int &tt);

void write_ansys_cell_data(const IPart * const &part,
    std::vector<double> &SV_P, std::vector<double> &SV_U,
    std::vector<double> &SV_V, std::vector<double> &SV_W,
    const std::string &FileName, const std::string &GroupName, const int &tt);

std::string get_ansys_result_name(const std::string &pre_fix, 
    const int &time, const std::string &suffix);

int main( int argc, char * argv[] )
{
  using namespace std;

  // Set number of threads and  print info of OpenMP
  SYS_T::print_omp_info();
  SYS_T::set_omp_num_threads();

  // Clean the potentially pre-existing hdf5 files in the job folder
  SYS_T::execute("rm -rf part_p*.h5");
  SYS_T::execute("rm -rf preprocessor_c2n_pla.h5");

  // Define basic problem settins
  constexpr int dofNum = 4; // degree-of-freedom for the physical problem
  constexpr int dofMat = 4; // degree-of-freedom in the matrix problem

  // Yaml options
  const std::string yaml_file("cell2node_preprocess.yml");

  // Check if the yaml file exist on disk
  SYS_T::file_check(yaml_file);

  YAML::Node paras = YAML::LoadFile( yaml_file );

  const std::string elemType_str      = paras["elem_type"].as<std::string>();
  const std::string geo_file          = paras["geo_file"].as<std::string>();
  const std::string part_file         = paras["part_file"].as<std::string>();
  const int cpu_size                  = paras["cpu_size"].as<int>();
  const int in_ncommon                = paras["in_ncommon"].as<int>();
  const bool isDualGraph              = paras["is_dualgraph"].as<bool>();
  const FEType elemType               = FE_T::to_FEType(elemType_str);

  // Print the command line arguments
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -elem_type: "<<elemType_str<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -isDualGraph: true \n";
  else cout<<" -isDualGraph: false \n";
  cout<<"---- Problem definition ----\n";
  cout<<" dofNum: "<<dofNum<<endl;
  cout<<" dofMat: "<<dofMat<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // Check if the vtu geometry files exist on disk
  SYS_T::file_check(geo_file); cout<<geo_file<<" found. \n";

  // Record the problem setting into a HDF5 file: preprocessor_cmd.h5
  hid_t cmd_file_id = H5Fcreate("preprocessor_c2n_pla.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_string("elemType", elemType_str);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);

  // Read the volumetric mesh file from the vtu file: geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;
  
  VTK_T::read_vtu_grid(geo_file, nFunc, nElem, ctrlPts, vecIEN);

  IIEN * IEN = new IEN_FEM(nElem, vecIEN);
  VEC_T::clean( vecIEN ); // clean the vector
  
  const int nLocBas = FE_T::to_nLocBas(elemType);

  SYS_T::print_fatal_if( IEN->get_nLocBas() != nLocBas, "Error: the nLocBas from the Mesh %d and the IEN %d classes do not match. \n", nLocBas, IEN->get_nLocBas() );

  // Call METIS to partition the mesh 
  IGlobal_Part * global_part = nullptr;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, nElem, nFunc, nLocBas, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( nElem, nFunc, "epart", "npart" );
  else SYS_T::print_fatal("ERROR: wrong cpu_size: %d \n", cpu_size);

  // Generate the new nodal numbering
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, nFunc);
  mnindex->write_hdf5("node_mapping");
  
  // Just use NodalBC to set up mapped ID
  std::vector<INodalBC *> NBC_list( dofMat, nullptr );

  std::vector<std::string> dir_list {}; // Empty

  NBC_list[0] = new NodalBC( nFunc );
  NBC_list[1] = new NodalBC( dir_list, nFunc );
  NBC_list[2] = new NodalBC( dir_list, nFunc );
  NBC_list[3] = new NodalBC( dir_list, nFunc );

  // // Reading cell data from ansys result
  // std::cout<<"=== Reading SV_P=== \n";
  // std::vector<double> SV_P = VTK_T::read_double_CellData( geo_file, "SV_P" );

  // std::cout<<"=== Reading SV_U === \n";
  // std::vector<double> SV_U = VTK_T::read_double_CellData( geo_file, "SV_U" );

  // std::cout<<"=== Reading SV_V === \n";
  // std::vector<double> SV_V = VTK_T::read_double_CellData( geo_file, "SV_V" );

  // std::cout<<"=== Reading SV_W === \n";
  // std::vector<double> SV_W = VTK_T::read_double_CellData( geo_file, "SV_W" );

  auto mytimer = SYS_T::make_unique<SYS_T::Timer>();

  std::string prefix = "FFF-1-01000-1-";
  std::string suffix_1 = ".cas.h5";
  std::string suffix_2 = ".dat.h5";

  for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();

    auto part = SYS_T::make_unique<Part_FEM>( 
        nElem, nFunc, nLocBas, global_part, mnindex, IEN,
        ctrlPts, proc_rank, cpu_size, elemType, 
        Field_Property(0, dofNum, true, "PLA") );
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    // write the part hdf5 file
    part -> write( part_file );

    part -> print_part_loadbalance_edgecut();

    // Partition Nodal BC and write to h5 file
    auto nbcpart = SYS_T::make_unique<NBC_Partition>(part.get(), mnindex, NBC_list);
    
    nbcpart -> write_hdf5( part_file );

    for(int tt = 500; tt < 18000; tt += 500)
    {
      std::string fileName = get_ansys_result_name(prefix, tt, suffix_1);

      hid_t f_id = H5Fopen( fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

      auto h5reader = SYS_T::make_unique<HDF5_Reader>(f_id);

      std::string GName_n("/meshes/1/nodes/zoneTopology");
      int nNode = h5reader->read_intScalar(GName_n.c_str(), "maxId");

      std::string GName_c("/meshes/1/nodes/coords");
      int n_col = 3;
      std::vector<double> coords = h5reader->read_doubleMatrix(GName_c.c_str(), "10", nNode, n_col);

      write_ansys_node_data(part.get(), mnindex, coords, part_file, "ctrlPts_loc", tt);

      H5Fclose( f_id );

      std::string fName = get_ansys_result_name(prefix, tt, suffix_2);

      hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

      auto h5r = SYS_T::make_unique<HDF5_Reader>(file_id);

      std::string GName_P("/results/1/phase-1/cells/SV_U");
      std::vector<double> SV_P = h5r->read_doubleVector(GName_P.c_str(), "1");

      std::string GName_U("/results/1/phase-1/cells/SV_U");
      std::vector<double> SV_U = h5r->read_doubleVector(GName_U.c_str(), "1");

      std::string GName_V("/results/1/phase-1/cells/SV_V");
      std::vector<double> SV_V = h5r->read_doubleVector(GName_V.c_str(), "1");

      std::string GName_W("/results/1/phase-1/cells/SV_W");
      std::vector<double> SV_W = h5r->read_doubleVector(GName_W.c_str(), "1");

      write_ansys_cell_data(part.get(), SV_P, SV_U, SV_V, SV_W, part_file, "Local_Elem", tt);

      H5Fclose( file_id );
    }
  }

  return EXIT_SUCCESS;
}

void write_ansys_node_data(const IPart * const &part,
  const Map_Node_Index * const &mnindex,
  const std::vector<double> &pts, const std::string &FileName,
  const std::string &GroupName, const int &tt)
{
  int nlocghonode = part->get_nlocghonode();
  std::vector<double> ctrlPts_x_loc(nlocghonode, 0.0);
  std::vector<double> ctrlPts_y_loc(nlocghonode, 0.0);
  std::vector<double> ctrlPts_z_loc(nlocghonode, 0.0);

  PERIGEE_OMP_PARALLEL_FOR
  for(int ii=0; ii<nlocghonode; ++ii)
  {
    int aux_index = part->get_local_to_global(ii); // new global index
    aux_index = mnindex->get_new2old(aux_index); // back to old global index
    ctrlPts_x_loc[ii] = pts[3*aux_index + 0];
    ctrlPts_y_loc[ii] = pts[3*aux_index + 1];
    ctrlPts_z_loc[ii] = pts[3*aux_index + 2];
  }

  const int cpu_rank = part->get_cpu_rank();

  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gopen(file_id, GroupName.c_str(), H5P_DEFAULT); 

  std::string subgroup_name = std::to_string(tt);

  hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(),
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  auto h5w = SYS_T::make_unique<HDF5_Writer>( file_id );

  h5w -> write_doubleVector( group_id, "ctrlPts_x_loc", ctrlPts_x_loc );

  h5w -> write_doubleVector( group_id, "ctrlPts_y_loc", ctrlPts_y_loc );

  h5w -> write_doubleVector( group_id, "ctrlPts_z_loc", ctrlPts_z_loc );

  H5Gclose( group_id );H5Gclose( g_id ); H5Fclose( file_id );
}

void write_ansys_cell_data(const IPart * const &part,
    std::vector<double> &SV_P,
    std::vector<double> &SV_U,
    std::vector<double> &SV_V,
    std::vector<double> &SV_W,
    const std::string &FileName,
    const std::string &GroupName,
    const int &tt)
{
  const int nlocalele = part->get_nlocalele();
  std::vector<double> local_p (nlocalele, 0.0);
  std::vector<double> local_u (nlocalele, 0.0);
  std::vector<double> local_v (nlocalele, 0.0);
  std::vector<double> local_w (nlocalele, 0.0);

  PERIGEE_OMP_PARALLEL_FOR
  for(int ii=0; ii<nlocalele; ++ii)
  {
    const int elem_index = part->get_elem_loc(ii);
    local_p[ii] = SV_P[elem_index];
    local_u[ii] = SV_U[elem_index];
    local_v[ii] = SV_V[elem_index];
    local_w[ii] = SV_W[elem_index];
  }

  const int cpu_rank = part->get_cpu_rank();

  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gopen(file_id, GroupName.c_str(), H5P_DEFAULT); 

  std::string subgroup_name = std::to_string(tt);

  hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(),
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  auto h5w = SYS_T::make_unique<HDF5_Writer>( file_id );

  h5w -> write_doubleVector( group_id, "local_p", local_p );

  h5w -> write_doubleVector( group_id, "local_u", local_u );

  h5w -> write_doubleVector( group_id, "local_v", local_v );

  h5w -> write_doubleVector( group_id, "local_w", local_w );

  H5Gclose( group_id );H5Gclose( g_id ); H5Fclose( file_id );
}

std::string get_ansys_result_name(const std::string &pre_fix, 
    const int &time, const std::string &suffix)
{
  std::string name = pre_fix;

  std::ostringstream time_index;

  std::string pre_time;

  if(time < 1000)
  {
    pre_time = "00";
  }
  else if(time >= 1000 && time < 10000)
  {
    pre_time = "0";
  }
  else if(time >= 10000)
  {
    pre_time = "";
  }
  else
   SYS_T::print_fatal("Error, wrong time index");

  time_index.str("");
  time_index<<time;
  pre_time += time_index.str();
   
  name += pre_time;

  name += suffix;

  return name;
}

// EOF