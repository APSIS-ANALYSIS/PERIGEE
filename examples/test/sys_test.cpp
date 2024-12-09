#include <iostream>
#include "FEType.hpp"

void print_res(bool result)
{
	if (result)
		std::cout << "Succeed" << std::endl;
	else
		std::cout << "Failed" << std::endl;
}

int main(int argc, char *argv[])
{
	// test FEType to string func()
	std::cout << "Test FEType to string" << std::endl;
	std::cout << "FEType::Tet4: ";
	print_res(std::string("Tet4")==FE_T::to_string(FEType::Tet4));
	
	std::cout << "FEType::Tet10: ";
	print_res(std::string("Tet10")==FE_T::to_string(FEType::Tet10));

	std::cout << "FEType::Hex8: ";
	print_res(std::string("Hex8")==FE_T::to_string(FEType::Hex8));

	std::cout << "FEType::Hex27: ";
	print_res(std::string("Hex27")==FE_T::to_string(FEType::Hex27));

	std::cout << "FEType::Tri3: ";
	print_res(std::string("Tri3")==FE_T::to_string(FEType::Tri3));

	std::cout << "FEType::Tri6: ";
	print_res(std::string("Tri6")==FE_T::to_string(FEType::Tri6));

	std::cout << "FEType::Quad4: ";
	print_res(std::string("Quad4")==FE_T::to_string(FEType::Quad4));

	std::cout << "FEType::Quad9: ";
	print_res(std::string("Quad9")==FE_T::to_string(FEType::Quad9));

	std::cout << "FEType::Tri3_der0: ";
	print_res(std::string("Tri3_der0")==FE_T::to_string(FEType::Tri3_der0));

	std::cout << "FEType::Tri6_der0: ";
	print_res(std::string("Tri6_der0")==FE_T::to_string(FEType::Tri6_der0));

	std::cout << "FEType::Quad4_der0: ";
	print_res(std::string("Quad4_der0")==FE_T::to_string(FEType::Quad4_der0));

	std::cout << "FEType::Quad9_der0: ";
	print_res(std::string("Quad9_der0")==FE_T::to_string(FEType::Quad9_der0));

	std::cout << "FEType::Unknown: ";
	print_res(std::string("Unknown")==FE_T::to_string(FEType::Unknown));

	std::cout << std::endl;

	std::cout << "Test string to FEType" << std::endl;
	std::cout << "Tet4: ";
	print_res(FEType::Tet4==FE_T::to_FEType(std::string("Tet4")));

	std::cout << "Tet10: ";
	print_res(FEType::Tet10==FE_T::to_FEType(std::string("Tet10")));

	std::cout << "Hex8: ";
	print_res(FEType::Hex8==FE_T::to_FEType(std::string("Hex8")));

	std::cout << "Hex27: ";
	print_res(FEType::Hex27==FE_T::to_FEType(std::string("Hex27")));

	std::cout << "Tri3: ";
	print_res(FEType::Tri3==FE_T::to_FEType(std::string("Tri3")));

	std::cout << "Tri6: ";
	print_res(FEType::Tri6==FE_T::to_FEType(std::string("Tri6")));

	std::cout << "Quad4: ";
	print_res(FEType::Quad4==FE_T::to_FEType(std::string("Quad4")));

	std::cout << "Quad9: ";
	print_res(FEType::Quad9==FE_T::to_FEType(std::string("Quad9")));

	std::cout << "Tri3_der0: ";
	print_res(FEType::Tri3_der0==FE_T::to_FEType(std::string("Tri3_der0")));

	std::cout << "Tri6_der0: ";
	print_res(FEType::Tri6_der0==FE_T::to_FEType(std::string("Tri6_der0")));

	std::cout << "Quad4_der0: ";
	print_res(FEType::Quad4_der0==FE_T::to_FEType(std::string("Quad4_der0")));

	std::cout << "Quad9_der0: ";
	print_res(FEType::Quad9_der0==FE_T::to_FEType(std::string("Quad9_der0")));

	std::cout << "Unknown: ";
	print_res(FEType::Unknown==FE_T::to_FEType(std::string("Unknown")));

	std::cout << "Illegal input: ";
	print_res(FEType::Unknown==FE_T::to_FEType(std::string("Illegal")));

	return 0;
}
// EOF
