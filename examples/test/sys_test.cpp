#include <vector>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "SymmTensor4_3D.hpp"
#include "HDF5_Reader.hpp"
#include "HDF5_Writer.hpp"
#include "HDF5_Tools.hpp"

int main(int argc, char *argv[])
{
	SYS_T::execute("rm -rf *.h5");
    
	hid_t infile_id   = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t group_id  = H5Gcreate(infile_id, "/dir", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	HDF5_Writer  test_h5w(infile_id);

	const std::string str_1 { "SUSTech" };
	const std::string str_2 { "MAE" };

	// write strings into h5 file
	test_h5w.write_string("string_1", str_1);
	test_h5w.write_string(group_id, "string_2", str_2);
	std::cout << "write string: " << str_1 << std::endl;
	std::cout << "write string into group " << str_2 << std::endl;

	const int intScalar_1 =  1;
	const int intScalar_2 = -2;

	// write int scalars into h5 file
	test_h5w.write_intScalar("intScalar_1", intScalar_1);
	test_h5w.write_intScalar(group_id, "intScalar_2", intScalar_2);
	std::cout << "write intScalar " << intScalar_1 << std::endl;
	std::cout << "write intScalar into group " << intScalar_2 << std::endl;

	const double doubleScalar_1 =  1.5;
	const double doubleScalar_2 = -9.3;

	// write double scalars into h5 file
	test_h5w.write_doubleScalar("doubleScalar_1", doubleScalar_1);
	test_h5w.write_doubleScalar(group_id, "doubleScalar_2", doubleScalar_2);
	std::cout << "write doubleScalar " << doubleScalar_1 << std::endl;
	std::cout << "write doubleScalar into group " << doubleScalar_2 << std::endl << std::endl;

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
	std::cout << "read doubleScalar from group: " << doubleScalar_2_r << std::endl;

	H5Fclose(outfile_id);

	return 0;
}
// EOF
