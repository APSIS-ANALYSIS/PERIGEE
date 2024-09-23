#include <chrono>
#include <thread>
#include <unistd.h>
#include <vector>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
#include "Tensor2_3D.hpp"
#include "SymmTensor4_3D.hpp"
#include "IEN_FEM.hpp"
#include "Mesh_Tet.hpp"
#include "Mesh_FEM.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "NodalBC.hpp"
#include "NodalBC_3D_inflow.hpp"
#include "ElemBC_3D.hpp"
#include "ElemBC_3D_outflow.hpp"
#include "ElemBC_3D_turbulence_wall_model.hpp"
#include "EBC_Partition_outflow.hpp"
#include "EBC_Partition_turbulence_wall_model.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "IQuadPts.hpp"
#include "QuadPts_Gauss_1D.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_UserDefined_Triangle.hpp"
#include "QuadPts_debug.hpp"
#include "QuadPts_vis_hex27.hpp"
#include "QuadPts_vis_hex8.hpp"
#include "QuadPts_vis_quad4.hpp"
#include "QuadPts_vis_quad9.hpp"
#include "QuadPts_vis_tet10.hpp"
#include "QuadPts_vis_tet10_v2.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "QuadPts_vis_tri6.hpp"
#include "FE_Tools.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include <iomanip>
#include "MaterialModel_Mixed_Elasticity.hpp"
#include "MaterialModel_ich_NeoHookean.hpp"
#include "MaterialModel_vol_Incompressible.hpp"
#include "MaterialModel_vol_ST91.hpp"
#include "MaterialModel_vol_M94.hpp"
#include "MaterialModel_ich_GOH06.hpp"
#include "MaterialModel_ich_GOH14.hpp"
#include "MaterialModel_ich_StVenant_Kirchhoff.hpp"
#include "HDF5_Reader.hpp"
#include "HDF5_Tools.hpp"
#include "yaml-cpp/yaml.h"

int main(int argc, char *argv[])
{
	// file_id
	hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

	// Yaml options
	const std::string yaml_file("ns_preprocess.yml");

	// Check if the yaml file exist on disk
	SYS_T::file_check(yaml_file);

	YAML::Node paras = YAML::LoadFile( yaml_file );

	const std::string geo_file          = paras["geo_file"].as<std::string>();
	const std::string sur_file_in_base  = paras["sur_file_in_base"].as<std::string>();
	const std::string part_file         = paras["part_file"].as<std::string>();
	const int num_inlet                 = paras["num_inlet"].as<int>();


	// group_id: file
	hid_t group_id_2 = H5Gcreate( cmd_file_id, "/ebc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ); 
	hid_t group_id_1 = H5Gcreate( cmd_file_id, "/Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ); 

	// HDF5_WRITE
	// write_string
	cmdh5w -> write_string("geo_file", geo_file);
	cmdh5w -> write_string("part_file", part_file);
	cmdh5w -> write_string(group_id_1, "sur_file_in_base", sur_file_in_base);

	// write_intScalar
	cmdh5w->write_intScalar("num_inlet", num_inlet);
	cmdh5w->write_intScalar(group_id_1, "num_inlet", num_inlet);

	// write_Tensor2_3D
	const Tensor2_3D foo {};
	cmdh5w -> write_Tensor2_3D(group_id_2, "foo", foo);
	cmdh5w -> write_Tensor2_3D("foo", foo);

	// HDF5_READER

	HDF5_Reader * cmdh5r = new HDF5_Reader(cmd_file_id);
	Tensor2_3D goo_1 = cmdh5r -> read_Tensor2_3D("/ebc", "foo");
	goo_1.print();

	H5Gclose( group_id_2); H5Gclose( group_id_1); 
	delete cmdh5r; 
	delete cmdh5w; H5Fclose( cmd_file_id );

	return EXIT_SUCCESS;
}
// EOF
