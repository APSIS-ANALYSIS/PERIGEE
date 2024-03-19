#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "PLocAssem_VMS_NS_SemiBDF1.hpp"
#include <ctime>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
  // Number of quadrature points for tets and triangles
  // const int nqp_tet = 5, nqp_tri = 4;

  // Number of quadrature points for hexs and quadrangles
  const int nqp_vol_1D = 2, nqp_sur_1D = 2;

  // Numeber of quadrature points
  const int nqp_vol = nqp_vol_1D * nqp_vol_1D * nqp_vol_1D;
  const int nqp_sur = nqp_sur_1D * nqp_sur_1D;

  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;
  IQuadPts   * quadv    = nullptr;
  IQuadPts   * quads    = nullptr;
  IPLocAssem * locAssem_ptr = nullptr;

  elementv = new FEAElement_Hex8( nqp_vol );  
  elements = new FEAElement_Quad4_3D_der0( nqp_sur );
  quadv    = new QuadPts_Gauss_Hex( nqp_vol_1D );
  quads    = new QuadPts_Gauss_Quad( nqp_sur_1D );

  double hex8_ctrl_x[8] = {0, 1, 3, 2, 1, 2, 4, 3};
  double hex8_ctrl_y[8] = {1, 0, 1, 2, 2, 1, 2, 3};
  double hex8_ctrl_z[8] = {0, 0, 0, 0, 1, 1, 1, 1};

  // fluid properties
  const double fluid_density = 1.0;
  const double fluid_mu = 1.0;
  const double bs_beta  = 0.0;
  
  // stablization parameters
  const double c_ct   = 1.0;
  const double c_tauc = 1.0;

  // sol & sol_n
  double sol[32], sol_0[32];
  srand(time(0));
  for(int ii=0; ii<32; ++ii)
  {
	sol[ii]   = ((double)rand()) / RAND_MAX * 2.0 - 1.0;
	sol_0[ii] = ((double)rand()) / RAND_MAX * 2.0 - 1.0;
  }
 
  // Hex8
  locAssem_ptr = new PLocAssem_VMS_NS_SemiBDF1(
      elementv->get_nLocBas(), quadv->get_num_quadPts(), elements->get_nLocBas(),
      fluid_density, fluid_mu, bs_beta, 601, c_ct, c_tauc );

  const double time = 0.0;
  const double dt   = 0.1;
  const double epsilon = 1e-3;
  const int index = 0;

  locAssem_ptr->Assem_Tangent_Residual(time, dt, sol_0, sol, 
  	   elementv, hex8_ctrl_x, hex8_ctrl_y, hex8_ctrl_z, quadv);
  
  double R_0[32] = {0};
  for(int ii=0; ii<32; ++ii)
  {	
  	 R_0[ii] = locAssem_ptr->Residual[ii]; 
  }
   
  sol[index] += epsilon;
  
  locAssem_ptr->Assem_Tangent_Residual(time, dt, sol_0, sol, 
  	  elementv, hex8_ctrl_x, hex8_ctrl_y, hex8_ctrl_z, quadv);
  
 for(int ii=0; ii<32; ++ii)
 {	
 	double  R_1   = locAssem_ptr->Residual[ii];
 	double  dR_de = locAssem_ptr->Tangent[32*ii+index];
	cout << abs(dR_de - (R_1 - R_0[ii])/epsilon) << endl;
 	//cout << (R_1 - R_0[ii])/epsilon << endl;
 	//cout << dR_de << endl;
 }


  delete elementv;
  delete elements;
  delete quadv;
  delete quads;
  delete locAssem_ptr;

  return 0;
}
