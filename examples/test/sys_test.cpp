#include <vector>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "SymmTensor4_3D.hpp"
#include "HDF5_Reader.hpp"
#include "HDF5_Writer.hpp"
#include "HDF5_Tools.hpp"

void print_array_double(const double* const array_d, const int length)
{
	for (int ii = 0; ii < length; ++ii)
		std::cout << array_d[ii] << " ";
	std::cout << std::endl;
}

void print_array_int(const int* const array_i, const int length)
{
	for (int ii = 0; ii < length; ++ii)
		std::cout << array_i[ii] << " ";
	std::cout << std::endl;
}

int main(int argc, char *argv[])
{
	SYS_T::execute("rm -rf *.h5");
    
	hid_t infile_id   = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t group_id  = H5Gcreate(infile_id, "/dir", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	HDF5_Writer  test_h5w(infile_id);

	// write strings into h5 file
	const std::string str_1 { "SUSTech" };
	const std::string str_2 { "MAE" };
	
	test_h5w.write_string("string_1", str_1);
	test_h5w.write_string(group_id, "string_2", str_2);
	std::cout << "write string: " << str_1 << std::endl;
	std::cout << "write string into group: " << str_2 << std::endl;

	// write int scalars into h5 file
	const int intScalar_1 =  1;
	const int intScalar_2 = -2;
	
	test_h5w.write_intScalar("intScalar_1", intScalar_1);
	test_h5w.write_intScalar(group_id, "intScalar_2", intScalar_2);
	std::cout << "write intScalar: " << intScalar_1 << std::endl;
	std::cout << "write intScalar into group: " << intScalar_2 << std::endl;

	// write double scalars into h5 file
	const double doubleScalar_1 =  1.5;
	const double doubleScalar_2 = -9.3;

	test_h5w.write_doubleScalar("doubleScalar_1", doubleScalar_1);
	test_h5w.write_doubleScalar(group_id, "doubleScalar_2", doubleScalar_2);
	std::cout << "write doubleScalar: " << doubleScalar_1 << std::endl;
	std::cout << "write doubleScalar into group: " << doubleScalar_2 << std::endl << std::endl;

	// write int array into h5 file
	int intVector_1[5] {1, 2, 3, 4, 5};
	int intVector_2[7] {-1, -3, -5, -7, -9, -11, -13};
	const std::vector<int> intVector_3 {10, 11, 12};
	const std::vector<int> intVector_4 {-20, -21, -22};

	test_h5w.write_intVector("intVector_1", intVector_1, 5);
	test_h5w.write_intVector(group_id, "intVector_2", intVector_2, 7);
	test_h5w.write_intVector("intVector_3", intVector_3);
	test_h5w.write_intVector(group_id, "intVector_4", intVector_4);

	std::cout << "write intVecotor_1: " << " ";
	print_array_int(intVector_1, 5);
	std::cout << "write intVector_2 into group: " << " ";
	print_array_int(intVector_2, 7);
	std::cout << "write intVector_3: " << " ";
	print_array_int(intVector_3.data(), intVector_3.size());
	std::cout << "write intVector_4 into group: " << " ";
	print_array_int(intVector_4.data(), intVector_4.size());

	// write double array into h5 file
	double doubleVector_1[6] {0.1, 0.2, 8.9, -9.8, 7.8, 2.3};
	double doubleVector_2[4] {-7.1, -5.6, -11.12, -98.7};
	const std::vector<double> doubleVector_3 {-21.21, 58.72, -92.11};
	const std::vector<double> doubleVector_4 {-14.14, 7.4, -32.1};

	test_h5w.write_doubleVector("doubleVector_1", doubleVector_1, 6);
	test_h5w.write_doubleVector(group_id, "doubleVector_2", doubleVector_2, 4);
	test_h5w.write_doubleVector("doubleVector_3", doubleVector_3);
	test_h5w.write_doubleVector(group_id, "doubleVector_4", doubleVector_4);

	std::cout << "write doubleVecotor_1: " << " ";
	print_array_double(doubleVector_1, 6);
	std::cout << "write doubleVecotor_2 into group: " << " ";
	print_array_double(doubleVector_2, 4);
	std::cout << "write doubleVecotor_3: " << " ";
	print_array_double(doubleVector_3.data(), doubleVector_3.size());
	std::cout << "write doubleVecotor_4 into group: " << " ";
	print_array_double(doubleVector_4.data(), doubleVector_4.size());
	std::cout << std::endl;

	H5Fclose(infile_id);
	H5Gclose(group_id);

	hid_t outfile_id = H5Fopen("test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

	HDF5_Reader test_h5r(outfile_id);

	//read strings from h5 file
	const std::string  str_1_r = test_h5r.read_string("/", "string_1");
	std::cout << "read string: " << str_1_r << std::endl;

	const std::string  str_2_r = test_h5r.read_string("/dir", "string_2");
	std::cout << "read string from group: " << str_2_r << std::endl;

	//read int scalars from h5 file
	const int  intScalar_1_r = test_h5r.read_intScalar("/", "intScalar_1");
	std::cout << "read intScalar: " << intScalar_1_r << std::endl;

	const int  intScalar_2_r = test_h5r.read_intScalar("/dir", "intScalar_2");
	std::cout << "read intScalar from group: " << intScalar_2_r << std::endl;

	//read double scalars from h5 file
	const double  doubleScalar_1_r = test_h5r.read_doubleScalar("/", "doubleScalar_1");
	std::cout << "read doubleScalar: " << doubleScalar_1_r << std::endl;

	const double  doubleScalar_2_r = test_h5r.read_doubleScalar("/dir", "doubleScalar_2");
	std::cout << "read doubleScalar from group: " << doubleScalar_2_r << std::endl << std::endl;

	// read int vector from h5 file
	const std::vector<int> intVector_1_r = test_h5r.read_intVector("/", "intVector_1");
	std::cout << "read intVector_1: ";
	print_array_int(intVector_1_r.data(), intVector_1_r.size());
		
	const std::vector<int> intVector_2_r = test_h5r.read_intVector("/dir", "intVector_2");
	std::cout << "read intVector_2 from group: ";
	print_array_int(intVector_2_r.data(), intVector_2_r.size());

	const std::vector<int> intVector_3_r = test_h5r.read_intVector("/", "intVector_3");
	std::cout << "read intVector_3: ";
	print_array_int(intVector_3_r.data(), intVector_3_r.size());
		
	const std::vector<int> intVector_4_r = test_h5r.read_intVector("/dir", "intVector_4");
	std::cout << "read intVector_4 from group: ";
	print_array_int(intVector_4_r.data(), intVector_4_r.size());

	// read double vector from h5 file
	const std::vector<double> doubleVector_1_r = test_h5r.read_doubleVector("/", "doubleVector_1");
	std::cout << "read doubleVector_1: ";
	print_array_double(doubleVector_1_r.data(), doubleVector_1_r.size());
		
	const std::vector<double> doubleVector_2_r = test_h5r.read_doubleVector("/dir", "doubleVector_2");
	std::cout << "read doubleVector_2 from group: ";
	print_array_double(doubleVector_2_r.data(), doubleVector_2_r.size());

	const std::vector<double> doubleVector_3_r = test_h5r.read_doubleVector("/", "doubleVector_3");
	std::cout << "read doubleVector_3: ";
	print_array_double(doubleVector_3_r.data(), doubleVector_3_r.size());
		
	const std::vector<double> doubleVector_4_r = test_h5r.read_doubleVector("/dir", "doubleVector_4");
	std::cout << "read doubleVector_4 from group: ";
	print_array_double(doubleVector_4_r.data(), doubleVector_4_r.size());

	H5Fclose(outfile_id);

	return 0;
}
// EOF
