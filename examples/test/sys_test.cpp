#include <chrono>
#include <thread>
#include <unistd.h>
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "Tensor4_3D.hpp"
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

int main(int argc, char *argv[])
{
	QuadPts_debug qpts{{1,1,1}, {1,2,3}};

	std::cout << "print_info():\n";
	qpts.print_info();

	std::cout << "get_dim():\n";
	std::cout << qpts.get_dim() << '\n';

	std::cout << "get_num_quadPts():\n";
	std::cout << qpts.get_num_quadPts() << '\n';

//	std::cout << "get_num_quadPts_x():\n";
//	std::cout << qpts.get_num_quadPts_x() << '\n';

//	std::cout << "get_num_quadPts_y():\n";
//	std::cout << qpts.get_num_quadPts_y() << '\n';

	std::cout << "get_qp(const int &ii):\n";
	std::cout << qpts.get_qp(1) << '\n';

	std::cout << "get_qp(const int &ii, const int &comp):\n";
	std::cout << qpts.get_qp(0,2) << '\n';

	std::cout << "get_qw(const int &ii):\n";
	std::cout << qpts.get_qw(1) << '\n';

  return EXIT_SUCCESS;
}

// EOF
