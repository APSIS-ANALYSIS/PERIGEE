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
#include "FE_Tools.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "AGlobal_Mesh_Info.hpp"
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
	SYS_T::execute("rm -rf part_p*.h5");
	SYS_T::execute("rm -rf preprocess_cmd.h5");
	
	//*******************************************************HDF5_WRITER
	hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t group_id_1  = H5Gcreate(cmd_file_id, "/Info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

	const std::string geo_file { "SUSTech" };
	const int num_inlet { 1 };

	// write_string
	cmdh5w -> write_string("geo_file", geo_file);
	cmdh5w -> write_string(group_id_1, "geo_file", "geo_file");

	// write_intScalar
	cmdh5w -> write_intScalar("num_inlet", num_inlet);
	cmdh5w -> write_intScalar(group_id_1, "num_inlet", num_inlet);

	// write_doubleScalar
	const double doubleScalar = 1.2;
	cmdh5w -> write_doubleScalar("doubleScalar", doubleScalar);
	cmdh5w -> write_doubleScalar(group_id_1, "doubleScalar", doubleScalar);

	// write_intVector
	const std::vector<int> intVector {1, 2, 3};
	cmdh5w -> write_intVector("intVector", intVector);

	// write_doubleVector
	const std::vector<double> doubleVector {1.0, 2.0, 3.0};
	cmdh5w -> write_doubleVector("doubleVector", doubleVector);

	// write_Vector_3
	const Vector_3 vector_3 {};
	cmdh5w -> write_Vector_3("Vector_3", vector_3);
	cmdh5w -> write_Vector_3(group_id_1, "Vector_3", vector_3);

	// write_Tensor2_3D
	Tensor2_3D tensor2_3D {};
	cmdh5w -> write_Tensor2_3D("Tensor2_3D", tensor2_3D);
	cmdh5w -> write_Tensor2_3D(group_id_1, "Tensor2_3D", tensor2_3D);

	// write_intMatrix
	int num_row = 3, num_col = 3;
	const std::vector<int> intMatrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	cmdh5w -> write_intMatrix(group_id_1, "intMatrix", intMatrix, num_row, num_col);

	// write_doubleMatrix
	std::vector<double> doubleMatrix = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
	cmdh5w -> write_doubleMatrix(group_id_1, "doubleMatrix", doubleMatrix, num_row, num_col);
	

	delete cmdh5w; H5Gclose( group_id_1 ); H5Fclose( cmd_file_id );

	// *************************************************HDF5_READER
	hid_t file_id = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

	HDF5_Reader * cmdh5r = new HDF5_Reader(file_id);

	// read_intScalar
	const	int intScalar = cmdh5r -> read_intScalar("/Info", "num_inlet");
	std::cout << "read_intScalar: " << intScalar << '\n';

	// read_doubleScalar
	const double doubleScalar_2 = cmdh5r -> read_doubleScalar("/", "doubleScalar");
	std::cout << "read_doubleScalar: " << doubleScalar_2 << '\n';

	// read_intVector
	const std::vector<int> intVector_2 = cmdh5r -> read_intVector("/", "intVector");
	std::cout << "read_intVector: ";
	VEC_T::print(intVector_2);

	// read_doubleVector
	const std::vector<double> doubleVector_2 = cmdh5r -> read_doubleVector("/", "doubleVector");
	std::cout << "read_doubleVector: "; 
	VEC_T::print(doubleVector_2);

	// read_Vector_3
	const Vector_3 foo = cmdh5r -> read_Vector_3("/", "Vector_3");
	std::cout << "read_Vector_3: ";
	foo.print();

	// read_Tensor2_3D
	const Tensor2_3D foo_2 = cmdh5r -> read_Tensor2_3D("/", "Tensor2_3D");
	std::cout << "read_Tensor2_3D: ";
	foo_2.print();

	// check_data
	const bool foo_3 = cmdh5r -> check_data("Vector_3");
	std::cout << "check_data: " << foo_3 << '\n'; 

	// read_intMatrix
	const std::vector<int> intMatrix_2 = cmdh5r -> read_intMatrix("/Info", "intMatrix", num_row, num_col);
	VEC_T::print(intMatrix_2);

	// read_doubleMatrix
	const std::vector<double> doubleMatrix_2 = cmdh5r -> read_doubleMatrix("/Info", "doubleMatrix", num_row, num_col);
	VEC_T::print(doubleMatrix_2);

	// read_string
	const std::string sustech = cmdh5r -> read_string("/", "geo_file");
	std::cout << "read_string: " << sustech << '\n';

	delete cmdh5r; 	H5Fclose( file_id );

	//**********************************************HDF5_Tools.cpp
	const int foo_4 = HDF5_T::read_intScalar("preprocessor_cmd.h5", "/", "num_inlet");
	std::cout << "HDF5_T::read_intScalar: " << foo_4 << '\n';

	const std::vector<int> foo_5 = HDF5_T::read_intVector("preprocessor_cmd.h5", "/", "Vector_3");
	std::cout << "HDF5_T::read_intVector: ";
	VEC_T::print(foo_5);

	return EXIT_SUCCESS;
}
// EOF
