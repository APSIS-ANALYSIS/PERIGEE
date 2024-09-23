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
	// HDF5_reader: new.h5
	hid_t prepcmd_file = H5Fopen("new.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

	// read_intScalar
	const int elemType = cmd_h5r -> read_intScalar("/Global_Mesh_Info","elemType");
	std::cout << "read_intScalar: " << '\n';
	std::cout << "elemType: " << elemType << '\n';
	std::cout << '\n';

	// check_data
	const bool isTagged = cmd_h5r -> check_data("/Global_Mesh_Info/elemType");
	std::cout << "check_data: " << '\n';
	std::cout << isTagged << '\n';

	// read_doubleScalar
	const double Inflow_full_area = cmd_h5r -> read_doubleScalar("/inflow/nbcid_0", "Inflow_full_area");
	std::cout << "read_doubleScalar: " << '\n';
	std::cout << "Inflow_full_area: " << std::setprecision(17) << Inflow_full_area << '\n';
	std::cout << '\n';

	// read_doubleVector
	const std::vector<double> doubleVector = cmd_h5r -> 
														read_doubleVector("/inflow/nbcid_0", "Outward_normal_vector");
	std::cout << "read_doubleVector: " << '\n';
	VEC_T::print(doubleVector, '\t');
	for (const double num : doubleVector)
		std::cout << std::setprecision(1) << std::fixed <<  num << '\t';
	std::cout << '\n' << '\n';

	// read_intVector
	const std::vector<int> intVector = cmd_h5r -> 
												 read_intVector("/inflow/nbcid_0", "LDN");
	std::cout << "read_intVector: " << '\n';
	for(const int num : intVector)
	{
		std::cout << num << '\t';
		if(num == 82)
			break;
	}
	std::cout << '\n' << '\n';

	// read_Vector_3
	std::cout << "read_Vector_3: " << '\n';
	const Vector_3	foo = cmd_h5r -> 
												read_Vector_3("/inflow/nbcid_0", "Outward_normal_vector");

	foo.print();
	std::cout << '\n' << '\n';

	// read_intMatrix
	std::cout << "read_intMatrix:\n ";
  int num_row, num_col;
  const std::vector<int> LIEN_vec = cmd_h5r -> read_intMatrix("LIEN", "LIEN", num_row, num_col);
	VEC_T::print(LIEN_vec, '\t');

	// read_doubleMatrix
	std::cout << "read_doubleMatrix:\n ";
  int num_row_d, num_col_d;
  const std::vector<double> LIEN_vec_d = cmd_h5r -> read_doubleMatrix("LIEN", "LIEN", num_row_d, num_col_d);
	VEC_T::print(LIEN_vec_d, '\t');

	// read_intScalar
	const int output = HDF5_T::read_intScalar("new.h5", "/Global_Mesh_Info", "elemType");
	std::cout << "HDF_T::read_intScalar: " << '\n';
	std::cout << output << '\n';

	// read_intVector
	const std::vector<int> foo_3 = HDF5_T::read_intVector("new.h5", "/inflow/nbcid_0", "LDN");
	std::cout << "HDF_T::read_intVector: " << '\n';
	for(const int num : intVector)
	{
		std::cout << num << '\t';
		if(num == 82)
			break;
	}
	std::cout << '\n';

	delete cmd_h5r; H5Fclose( prepcmd_file );

	// HDF5_writer: preprocessor_cmd.h5
	
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

	// write_string
  cmdh5w -> write_string("geo_file", geo_file);
	
  // group_id: file
  hid_t group_id = H5Gcreate( cmd_file_id, "/Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT ); 
	cmdh5w -> write_string(group_id, "sur_file_in_base", sur_file_in_base);

  H5Gclose( group_id );

	delete cmdh5w; H5Fclose( cmd_file_id );

  return EXIT_SUCCESS;
}

// EOF
