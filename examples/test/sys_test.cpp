#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include "FEType.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad4.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "FEAElement_Quad9.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle3.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Triangle6.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"

void writeCoordinatesToFile(const char* filename, double* x, double* y, double* z, size_t size) {
    std::ofstream file(filename, std::ios::app);
	if (!file.is_open()) {
        std::cerr << "Invalid filename: " << filename << std::endl;
        return;
    }
    file << std::fixed << std::setprecision(15);

    for (size_t i = 0; i < size; ++i) {
        file << x[i];
        if (i < size - 1) file << " ";
    }
    file << std::endl;

    for (size_t i = 0; i < size; ++i) {
        file << y[i];
        if (i < size - 1) file << " ";
    }
    file << std::endl;

    for (size_t i = 0; i < size; ++i) {
        file << z[i];
        if (i < size - 1) file << " ";
    }
    file << std::endl;

    file.close();
}

void writeCoordinatesToFile(const char* filename, double* x, double* y, size_t size) {
    std::ofstream file(filename, std::ios::app);
	if (!file.is_open()) {
        std::cerr << "Invalid filename: " << filename << std::endl;
        return;
    }
	file << std::fixed << std::setprecision(15);

    for (size_t i = 0; i < size; ++i) {
        file << x[i];
        if (i < size - 1) file << " ";
    }
    file << std::endl;

    for (size_t i = 0; i < size; ++i) {
        file << y[i];
        if (i < size - 1) file << " ";
    }
    file << std::endl;

    file.close();
}

void readCoordinatesFromFile(const char* filename, size_t startLine, double* x, double* y, double* z, size_t size) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Invalid filename: " << filename << std::endl;
        return;
    }

    std::string line;
    size_t currentLine = 0;

    while (currentLine < startLine && std::getline(file, line)) {
        ++currentLine;
    }

    if (currentLine != startLine) {
		std::cerr << "Invalid startline: " << startLine << std::endl;
        return;
    }

    if (std::getline(file, line)) {
        std::istringstream iss(line);
        for (size_t i = 0; i < size; ++i) {
            if (!(iss >> x[i])) {
				std::cerr << "Error occurs in reading x coordinate" << std::endl;
                return;
            }
        }
    } else {
		std::cerr << "Error: Empty line for x coordinate" << std::endl;
        return;
    }

    if (std::getline(file, line)) {
        std::istringstream iss(line);
        for (size_t i = 0; i < size; ++i) {
            if (!(iss >> y[i])) {
				std::cerr << "Error occurs in reading y coordinate" << std::endl;
                return;
            }
        }
    } else {
		std::cerr << "Error: Empty line for y coordinate" << std::endl;
        return;
    }

    if (std::getline(file, line)) {
        std::istringstream iss(line);
        for (size_t i = 0; i < size; ++i) {
            if (!(iss >> z[i])) {
				std::cerr << "Error occurs in reading z coordinate" << std::endl;
                return;
            }
        }
    } else {
		std::cerr << "Error: Empty line for z coordinate" << std::endl;
        return;
    }

    file.close();
}

void readCoordinatesFromFile(const char* filename, size_t startLine, double* x, double* y, size_t size) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Invalid filename: " << filename << std::endl;
        return;
    }

    std::string line;
    size_t currentLine = 0;

    while (currentLine < startLine && std::getline(file, line)) {
        ++currentLine;
    }

    if (currentLine != startLine) {
		std::cerr << "Invalid startline: " << startLine << std::endl;
        return;
    }

    if (std::getline(file, line)) {
        std::istringstream iss(line);
        for (size_t i = 0; i < size; ++i) {
            if (!(iss >> x[i])) {
				std::cerr << "Error occurs in reading x coordinate" << std::endl;
                return;
            }
        }
    } else {
		std::cerr << "Error: Empty line for x coordinate" << std::endl;
        return;
    }

    if (std::getline(file, line)) {
        std::istringstream iss(line);
        for (size_t i = 0; i < size; ++i) {
            if (!(iss >> y[i])) {
				std::cerr << "Error occurs in reading y coordinate" << std::endl;
                return;
            }
        }
    } else {
		std::cerr << "Error: Empty line for y coordinate" << std::endl;
        return;
    }

    file.close();
}

void generateHex8_Quad4_der0(double * x, double * y, double * z,
	double * x_der0, double * y_der0, double * z_der0) 
{
    double baseX[8] = {0, 1, 1, 0, 0, 1, 1, 0};
    double baseY[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    double baseZ[8] = {0, 0, 0, 0, 1, 1, 1, 1};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_length(0.0, 0.3);
    std::uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    for (size_t i = 0; i < 8; ++i) 
    {
        double length = dist_length(gen);
        double theta = dist_angle(gen);
        double phi = dist_angle(gen);

        double dx = length * std::sin(phi) * std::cos(theta);
        double dy = length * std::sin(phi) * std::sin(theta);
        double dz = length * std::cos(phi);

        x[i] = baseX[i] + dx;
        y[i] = baseY[i] + dy;
        z[i] = baseZ[i] + dz;
    }

    for (size_t i = 0; i < 4; ++i)
    {
	    x_der0[i] = x[i];
	    y_der0[i] = y[i];
	    z_der0[i] = z[i];
    }
}

void generateHex27_Quad9_der0(double * x, double * y, double * z,
	double * x_der0, double * y_der0, double * z_der0) 
{
    double baseX[27] = {0, 1, 1, 0, 
                        0, 1, 1, 0, 
                        0.5, 1, 0.5, 0, 
                        0.5, 1, 0.5, 0, 
                        0, 1, 1, 0, 
                        0, 0.5, 1, 0.5, 
                        0.5, 0.5, 0.5};
    double baseY[27] = {0, 0, 1, 1, 
                        0, 0, 1, 1, 
                        0, 0.5, 1, 0.5, 
                        0, 0.5, 1, 0.5, 
                        0, 0, 1, 1, 
                        0.5, 0.5, 0, 1,
                        0.5, 0.5, 0.5};
    double baseZ[27] = {0, 0, 0, 0, 
                        1, 1, 1, 1,
                        0, 0, 0, 0,
                        1, 1, 1, 1,
                        0.5, 0.5, 0.5, 0.5,
                        0.5, 0.5, 0.5, 0.5,
                        0, 0.5, 1};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_length(0.0, 0.15);
    std::uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    for (size_t i = 0; i < 27; ++i)
	{
        double length = dist_length(gen);
        double theta = dist_angle(gen);
        double phi = dist_angle(gen);

        double dx = length * std::sin(phi) * std::cos(theta);
        double dy = length * std::sin(phi) * std::sin(theta);
        double dz = length * std::cos(phi);

        x[i] = baseX[i] + dx;
        y[i] = baseY[i] + dy;
        z[i] = baseZ[i] + dz;
    }

    for (size_t i = 0; i < 4; ++i)
    {
	    x_der0[i] = x[i];
	    y_der0[i] = y[i];
	    z_der0[i] = z[i];

	    x_der0[4+i] = x[8+i];
	    y_der0[4+i] = y[8+i];
	    z_der0[4+i] = z[8+i];
    }

    x_der0[8] = x[24];
    y_der0[8] = y[24];
    z_der0[8] = z[24];
}

void generateTet4_Tri3_der0(double * x, double * y, double * z,
	double * x_der0, double * y_der0, double * z_der0)
{
    double baseX[4] = {0, 1, 0, 0};
    double baseY[4] = {0, 0, 1, 0};
    double baseZ[4] = {0, 0, 0, 1};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_length(0.0, 0.15);
    std::uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    for (size_t i = 0; i < 4; ++i)
	{
        double length = dist_length(gen);
        double theta = dist_angle(gen);
        double phi = dist_angle(gen);

        double dx = length * std::sin(phi) * std::cos(theta);
        double dy = length * std::sin(phi) * std::sin(theta);
        double dz = length * std::cos(phi);

        x[i] = baseX[i] + dx;
        y[i] = baseY[i] + dy;
        z[i] = baseZ[i] + dz;
    }

    for (size_t i = 0; i < 3; ++i)
    {
        x_der0[i] = x[i];
        y_der0[i] = y[i];
        z_der0[i] = z[i];
    }
}

void generateTet10_Tri6_der0(double * x, double * y, double * z,
	double * x_der0, double * y_der0, double * z_der0)
{
    double baseX[10] = {0, 1, 0, 0, 0.5, 0.5, 0, 0, 0.5, 0};
    double baseY[10] = {0, 0, 1, 0, 0, 0.5, 0.5, 0, 0, 0.5};
    double baseZ[10] = {0, 0, 0, 1, 0, 0, 0, 0.5, 0.5, 0.5};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_length(0.0, 0.1);
    std::uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    for (size_t i = 0; i < 10; ++i)
	{
        double length = dist_length(gen);
        double theta = dist_angle(gen);
        double phi = dist_angle(gen);

        double dx = length * std::sin(phi) * std::cos(theta);
        double dy = length * std::sin(phi) * std::sin(theta);
        double dz = length * std::cos(phi);

        x[i] = baseX[i] + dx;
        y[i] = baseY[i] + dy;
        z[i] = baseZ[i] + dz;
    }

    for (size_t i = 0; i < 3; ++i)
    {
        x_der0[i] = x[i];
        y_der0[i] = y[i];
        z_der0[i] = z[i];

        x_der0[3+i] = x[4+i];
        y_der0[3+i] = y[4+i];
        z_der0[3+i] = z[4+i];
    }
}

void generateQuad4(double * x, double * y) {
    double baseX[4] = {0, 1, 1, 0};
    double baseY[4] = {0, 0, 1, 1};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_length(0.0, 0.2);
    std::uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    for (size_t i = 0; i < 4; ++i)
	{
        double length = dist_length(gen);
        double theta = dist_angle(gen);

        double dx = length * std::cos(theta);
        double dy = length * std::sin(theta);

        x[i] = baseX[i] + dx;
        y[i] = baseY[i] + dy;
    }
}

void generateQuad9(double * x, double * y) {
    double baseX[9] = {0, 1, 1, 0, 0.5, 1, 0.5, 0, 0.5};
    double baseY[9] = {0, 0, 1, 1, 0, 0.5, 1, 0.5, 0.5};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_length(0.0, 0.1);
    std::uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    for (size_t i = 0; i < 9; ++i)
	{
        double length = dist_length(gen);
        double theta = dist_angle(gen);

        double dx = length * std::cos(theta);
        double dy = length * std::sin(theta);

        x[i] = baseX[i] + dx;
        y[i] = baseY[i] + dy;
    }
}

void generateTri3(double * x, double * y) {
    double baseX[3] = {0, 1, 0};
    double baseY[3] = {0, 0, 1};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_length(0.0, 0.2);
    std::uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    for (size_t i = 0; i < 3; ++i)
	{
        double length = dist_length(gen);
        double theta = dist_angle(gen);

        double dx = length * std::cos(theta);
        double dy = length * std::sin(theta);

        x[i] = baseX[i] + dx;
        y[i] = baseY[i] + dy;
    }
}

void generateTri6(double * x, double * y) {
    double baseX[6] = {0, 1, 0, 0.5, 0.5, 0};
    double baseY[6] = {0, 0, 1, 0, 0.5, 0.5};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist_length(0.0, 0.1);
    std::uniform_real_distribution<double> dist_angle(0.0, 2 * M_PI);

    for (size_t i = 0; i < 6; ++i)
	{
        double length = dist_length(gen);
        double theta = dist_angle(gen);

        double dx = length * std::cos(theta);
        double dy = length * std::sin(theta);

        x[i] = baseX[i] + dx;
        y[i] = baseY[i] + dy;
    }
}

void outputVector(double * const vec, size_t size)
{
	for (size_t i = 0; i < size; ++i)
	{
        std::cout << vec[i];
        if (i < size - 1) std::cout << " ";
    }
    std::cout << std::endl;
}

void outputVector(const std::vector<double> vec, size_t size)
{
	for (size_t i = 0; i < size; ++i)
	{
        std::cout << vec[i];
        if (i < size - 1) std::cout << " ";
    }
    std::cout << std::endl;
}

void outputArray(const std::array<double, 9> arr)
{
	for (size_t i = 0; i < 9; ++i)
	{
        std::cout << arr[i];
        if (i < 8) std::cout << " ";
    }
    std::cout << std::endl;
}

void outputArray(const std::array<double, 4> arr)
{
	for (size_t i = 0; i < 4; ++i)
	{
        std::cout << arr[i];
        if (i < 3) std::cout << " ";
    }
    std::cout << std::endl;
}

void testMemberFunction(FEAElement * const &elementv, FEAElement * const &elements,
    const double * const &ctrl_x, const double * const &ctrl_y, const double * const &ctrl_z,
	const double * const &ctrl_x_der0, const double * const &ctrl_y_der0, const double * const &ctrl_z_der0)
{
	int nqp_vol_1D = 0, nqp_sur_1D = 0;
	int nqp_vol = 0, nqp_sur = 0;

	IQuadPts * quadv = nullptr;
	IQuadPts * quads = nullptr;

	if( elementv->get_Type() == FEType::Tet4 && elements->get_Type() == FEType::Tri3_der0 )
    {
    	nqp_vol = 5;
		nqp_sur = 4;

    	quadv = new QuadPts_Gauss_Tet( nqp_vol );
		quads = new QuadPts_Gauss_Triangle( nqp_sur );
  	}
 	else if( elementv->get_Type() == FEType::Tet10 && elements->get_Type() == FEType::Tri6_der0 )
  	{
		nqp_vol = 29;
		nqp_sur = 13;

    	quadv = new QuadPts_Gauss_Tet( nqp_vol );
		quads = new QuadPts_Gauss_Triangle( nqp_sur );
  	}
  	else if( elementv->get_Type() == FEType::Hex8 && elements->get_Type() == FEType::Quad4_der0 )
  	{
		nqp_vol_1D = 3;
		nqp_sur_1D = 2;

    	quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
		quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  	}
  	else if( elementv->get_Type() == FEType::Hex27 && elements->get_Type() == FEType::Quad9_der0 )
  	{
		nqp_vol_1D = 4;
		nqp_sur_1D = 3;

    	quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
		quads = new QuadPts_Gauss_Quad( nqp_sur_1D );
  	}
  	else SYS_T::print_fatal("Error: Element type not supported.\n");

	const int dim = elementv->get_elemDim();
	std::cout << dim << std::endl;

	const int numQuapts = elementv->get_numQuapts();
	std::cout << numQuapts << std::endl;

	const int nLocBas = elementv->get_nLocBas();
	std::cout << nLocBas << std::endl;

	elementv->print_info();
	elementv->buildBasis(quadv, ctrl_x, ctrl_y, ctrl_z);
	
	const double h = elementv->get_h(ctrl_x, ctrl_y, ctrl_z);
	std::cout << h << std::endl;
	
	double * R = new double[nLocBas];
	double * R_dx = new double[nLocBas];
	double * R_dy = new double[nLocBas];
	double * R_dz = new double[nLocBas];
	double * R_dxx = new double[nLocBas];
	double * R_dyy = new double[nLocBas];
	double * R_dzz = new double[nLocBas];
	double * R_dxy = new double[nLocBas];
	double * R_dxz = new double[nLocBas];
	double * R_dyz = new double[nLocBas];

	double * dx_dr = new double[9];
	double * dr_dx = new double[9];

	double detJac = 0.0;
	
	for(int qua=0; qua<numQuapts; ++qua)
	{	
		elementv->get_R(qua, R);
		outputVector(R, nLocBas);

		std::vector<double> RR = elementv->get_R(qua);
		outputVector(RR, nLocBas);

		elementv->get_gradR(qua, R_dx, R_dy, R_dz);
		outputVector(R_dx, nLocBas);
		outputVector(R_dy, nLocBas);
		outputVector(R_dz, nLocBas);

		elementv->get_3D_R_dR_d2R(qua, R, R_dx, R_dy, R_dz,
						R_dxx, R_dyy, R_dzz,
						R_dxy, R_dxz, R_dyz);
		outputVector(R, nLocBas);
		outputVector(R_dx, nLocBas);
		outputVector(R_dy, nLocBas);
		outputVector(R_dz, nLocBas);
		outputVector(R_dxx, nLocBas);
		outputVector(R_dyy, nLocBas);
		outputVector(R_dzz, nLocBas);
		outputVector(R_dxy, nLocBas);
		outputVector(R_dxz, nLocBas);
		outputVector(R_dyz, nLocBas);

		elementv->get_3D_R_gradR_LaplacianR(qua, R, R_dx, R_dy, R_dz,
								R_dxx, R_dyy, R_dzz);
		outputVector(R, nLocBas);
		outputVector(R_dx, nLocBas);
		outputVector(R_dy, nLocBas);
		outputVector(R_dz, nLocBas);
		outputVector(R_dxx, nLocBas);
		outputVector(R_dyy, nLocBas);
		outputVector(R_dzz, nLocBas);

		elementv->get_Jacobian(qua, dx_dr);

		outputVector(dx_dr, 9);
		
		const std::array<double, 9> Jacobian = elementv->get_Jacobian(qua);
		outputArray(Jacobian);

		elementv->get_invJacobian(qua, dr_dx);
		outputVector(dr_dx, 9);

		const std::array<double, 9> invJacobian = elementv->get_invJacobian(qua);
		outputArray(invJacobian);

		detJac = elementv->get_detJac(qua);
		std::cout << detJac << std::endl;

		double area = 0.0;
		const Vector_3 outwardnormal = elementv->get_2d_normal_out(qua, area);

		std::cout << area << std::endl;
		std::cout << outwardnormal.x() << ' ' << outwardnormal.y() << ' ' << outwardnormal.z() << std::endl;
	}

	const int dim_der0 = elements->get_elemDim();
	std::cout << dim_der0 << std::endl;

	const int numQuapts_der0 = elements->get_numQuapts();
	std::cout << numQuapts_der0 << std::endl;

	const int nLocBas_der0 = elements->get_nLocBas();
	std::cout << nLocBas_der0 << std::endl;

	elements->print_info();
	elements->buildBasis(quads, ctrl_x_der0, ctrl_y_der0, ctrl_z_der0);
	
	double * R_der0 = new double[nLocBas_der0];

	double detJac_der0 = 0.0;

	for(int qua=0; qua<numQuapts_der0; ++qua)
	{
		elements->get_R(qua, R);
		outputVector(R, nLocBas_der0);

		std::vector<double> RR = elements->get_R(qua);
		outputVector(RR, nLocBas_der0);

		double area = 0.0;
		const Vector_3 outwardnormal = elements->get_2d_normal_out(qua, area);

		std::cout << area << std::endl;
		std::cout << outwardnormal.x() << ' ' << outwardnormal.y() << ' ' << outwardnormal.z() << std::endl;

		detJac_der0 = elements->get_detJac(qua);
		std::cout << detJac_der0 << std::endl;

		const Vector_3 dx_dr_der0 = elements->get_dx_dr(qua, ctrl_x_der0, ctrl_y_der0, ctrl_z_der0);
		std::cout << dx_dr_der0.x() << dx_dr_der0.y() << dx_dr_der0.z() << std::endl;

		const Vector_3 dx_ds_der0 = elements->get_dx_ds(qua, ctrl_x_der0, ctrl_y_der0, ctrl_z_der0);
		std::cout << dx_ds_der0.x() << dx_ds_der0.y() << dx_ds_der0.z() << std::endl;

		const Vector_3 d2x_drr_der0 = elements->get_d2x_drr(qua, ctrl_x_der0, ctrl_y_der0, ctrl_z_der0);
		std::cout << d2x_drr_der0.x() << d2x_drr_der0.y() << d2x_drr_der0.z() << std::endl;

		const Vector_3 d2x_dss_der0 = elements->get_d2x_dss(qua, ctrl_x_der0, ctrl_y_der0, ctrl_z_der0);
		std::cout << d2x_dss_der0.x() << d2x_dss_der0.y() << d2x_dss_der0.z() << std::endl;

		const Vector_3 d2x_drs_der0 = elements->get_d2x_drs(qua, ctrl_x_der0, ctrl_y_der0, ctrl_z_der0);
		std::cout << d2x_drs_der0.x() << d2x_drs_der0.y() << d2x_drs_der0.z() << std::endl;
	}

	for(int ii=0; ii<nLocBas_der0; ++ii)
	{
		const Vector_3 int_pt(0.5, 0.5, 1.0);
		const Vector_3 sur_pt(ctrl_x_der0[ii], ctrl_y_der0[ii], ctrl_z_der0[ii]);

		double area = 0.0;
		const Vector_3 outwardnormal = elements->get_normal_out(ii, sur_pt, int_pt, area);

		std::cout << area << std::endl;
		std::cout << outwardnormal.x() << ' ' << outwardnormal.y() << ' ' << outwardnormal.z() << std::endl;
	}

	delete quadv; quadv = nullptr;
	delete quads; quads = nullptr;

	delete [] R; R = nullptr;
	delete [] R_dx; R_dx = nullptr;
	delete [] R_dy; R_dy = nullptr;
	delete [] R_dz; R_dz = nullptr;
	delete [] R_dxx; R_dxx = nullptr;
	delete [] R_dyy; R_dyy = nullptr;
	delete [] R_dzz; R_dzz = nullptr;
	delete [] R_dxy; R_dxy = nullptr;
	delete [] R_dxz; R_dxz = nullptr;
	delete [] R_dyz; R_dyz = nullptr;

	delete [] dx_dr; dx_dr = nullptr;
	delete [] dr_dx; dr_dx = nullptr;

	delete [] R_der0; R_der0 = nullptr;
}

void testMemberFunction(FEAElement * const &element, const double * const &ctrl_x,
	const double * const &ctrl_y)
{
	int nqp_sur_1D = 0;
	int nqp_sur = 0;

    IQuadPts * quad = nullptr;

	if( element->get_Type() == FEType::Tri3 )
    {
		nqp_sur = 4;
    	quad = new QuadPts_Gauss_Triangle( nqp_sur );
  	}
 	else if( element->get_Type() == FEType::Tri6 )
  	{
		nqp_sur = 13;
    	quad = new QuadPts_Gauss_Triangle( nqp_sur );
  	}
  	else if( element->get_Type() == FEType::Quad4 )
  	{
		nqp_sur_1D = 2;
    	quad = new QuadPts_Gauss_Quad( nqp_sur_1D );
  	}
  	else if( element->get_Type() == FEType::Quad9 )
  	{
		nqp_sur_1D = 3;
    	quad = new QuadPts_Gauss_Quad( nqp_sur_1D );
  	}
  	else SYS_T::print_fatal("Error: Element type not supported.\n");

	const int dim = element->get_elemDim();
	std::cout << dim << std::endl;

	const int numQuapts = element->get_numQuapts();
	std::cout << numQuapts << std::endl;

	const int nLocBas = element->get_nLocBas();
	std::cout << nLocBas << std::endl;

	element->print_info();
	element->buildBasis(quad, ctrl_x, ctrl_y);
	
	const double h = element->get_h(ctrl_x, ctrl_y);
	std::cout << h << std::endl;
	
	double * R = new double[nLocBas];
	double * R_dx = new double[nLocBas];
	double * R_dy = new double[nLocBas];
	double * R_dxx = new double[nLocBas];
	double * R_dyy = new double[nLocBas];
	double * R_dxy = new double[nLocBas];

	double * dx_dr = new double[4];
	double * dr_dx = new double[4];

	double detJac = 0.0;

	for(int qua=0; qua<numQuapts; ++qua)
	{
		element->get_R(qua, R);
		outputVector(R, nLocBas);

		std::vector<double> RR = element->get_R(qua);
		outputVector(RR, nLocBas);

		element->get_gradR(qua, R_dx, R_dy);
		outputVector(R_dx, nLocBas);
		outputVector(R_dy, nLocBas);

		element->get_2D_R_dR_d2R(qua, R, R_dx, R_dy,
						R_dxx, R_dyy, R_dxy);
		outputVector(R, nLocBas);
		outputVector(R_dx, nLocBas);
		outputVector(R_dy, nLocBas);
		outputVector(R_dxx, nLocBas);
		outputVector(R_dyy, nLocBas);
		outputVector(R_dxy, nLocBas);

		element->get_Jacobian(qua, dx_dr);
		outputVector(dx_dr, 4);
		
		element->get_invJacobian(qua, dr_dx);
		outputVector(dr_dx, 4);

		detJac = element->get_detJac(qua);
		std::cout << detJac << std::endl;
	}

	delete quad; quad = nullptr;

	delete [] R; R = nullptr;
	delete [] R_dx; R_dx = nullptr;
	delete [] R_dy; R_dy = nullptr;
	delete [] R_dxx; R_dxx = nullptr;
	delete [] R_dyy; R_dyy = nullptr;
	delete [] R_dxy; R_dxy = nullptr;

	delete [] dx_dr; dx_dr = nullptr;
	delete [] dr_dx; dr_dx = nullptr;
}

int main(int argc, char *argv[])
{
	std::cout << std::fixed << std::setprecision(15);

	const char* outfile = "coordinates.txt";

	double * ctrl_x_hex8 = new double[8];
	double * ctrl_y_hex8 = new double[8];
	double * ctrl_z_hex8 = new double[8];
	double * ctrl_x_quad4_3d = new double[4];
	double * ctrl_y_quad4_3d = new double[4];
	double * ctrl_z_quad4_3d = new double[4];

	generateHex8_Quad4_der0(ctrl_x_hex8, ctrl_y_hex8, ctrl_z_hex8,
		ctrl_x_quad4_3d, ctrl_y_quad4_3d, ctrl_z_quad4_3d);
	writeCoordinatesToFile(outfile, ctrl_x_hex8, ctrl_y_hex8, ctrl_z_hex8, 8);
	writeCoordinatesToFile(outfile, ctrl_x_quad4_3d, ctrl_y_quad4_3d, ctrl_z_quad4_3d, 4);
	
	double * ctrl_x_hex27 = new double[27];
	double * ctrl_y_hex27 = new double[27];
	double * ctrl_z_hex27 = new double[27];
	double * ctrl_x_quad9_3d = new double[9];
	double * ctrl_y_quad9_3d = new double[9];
	double * ctrl_z_quad9_3d = new double[9];

	generateHex27_Quad9_der0(ctrl_x_hex27, ctrl_y_hex27, ctrl_z_hex27,
		ctrl_x_quad9_3d, ctrl_y_quad9_3d, ctrl_z_quad9_3d);
	writeCoordinatesToFile(outfile, ctrl_x_hex27, ctrl_y_hex27, ctrl_z_hex27, 27);
	writeCoordinatesToFile(outfile, ctrl_x_quad9_3d, ctrl_y_quad9_3d, ctrl_z_quad9_3d, 9);

	double * ctrl_x_tet4 = new double[4];
	double * ctrl_y_tet4 = new double[4];
	double * ctrl_z_tet4 = new double[4];
	double * ctrl_x_tri3_3d = new double[3];
	double * ctrl_y_tri3_3d = new double[3];
	double * ctrl_z_tri3_3d = new double[3];

	generateTet4_Tri3_der0(ctrl_x_tet4, ctrl_y_tet4, ctrl_z_tet4,
		ctrl_x_tri3_3d, ctrl_y_tri3_3d, ctrl_z_tri3_3d);
	writeCoordinatesToFile(outfile, ctrl_x_tet4, ctrl_y_tet4, ctrl_z_tet4, 4);
	writeCoordinatesToFile(outfile, ctrl_x_tri3_3d, ctrl_y_tri3_3d, ctrl_z_tri3_3d, 3);

	double * ctrl_x_tet10 = new double[10];
	double * ctrl_y_tet10 = new double[10];
	double * ctrl_z_tet10 = new double[10];
	double * ctrl_x_tri6_3d = new double[6];
	double * ctrl_y_tri6_3d = new double[6];
	double * ctrl_z_tri6_3d = new double[6];

	generateTet10_Tri6_der0(ctrl_x_tet10, ctrl_y_tet10, ctrl_z_tet10,
		ctrl_x_tri6_3d, ctrl_y_tri6_3d, ctrl_z_tri6_3d);
	writeCoordinatesToFile(outfile, ctrl_x_tet10, ctrl_y_tet10, ctrl_z_tet10, 10);
	writeCoordinatesToFile(outfile, ctrl_x_tri6_3d, ctrl_y_tri6_3d, ctrl_z_tri6_3d, 6);

	double * ctrl_x_quad4 = new double[4];
	double * ctrl_y_quad4 = new double[4];

	generateQuad4(ctrl_x_quad4, ctrl_y_quad4);
	writeCoordinatesToFile(outfile, ctrl_x_quad4, ctrl_y_quad4, 4);

	double * ctrl_x_quad9 = new double[9];
	double * ctrl_y_quad9 = new double[9];

	generateQuad9(ctrl_x_quad9, ctrl_y_quad9);
	writeCoordinatesToFile(outfile, ctrl_x_quad9, ctrl_y_quad9, 9);

	double * ctrl_x_tri3 = new double[3];
	double * ctrl_y_tri3 = new double[3];

	generateTri3(ctrl_x_tri3, ctrl_y_tri3);
	writeCoordinatesToFile(outfile, ctrl_x_tri3, ctrl_y_tri3, 3);

	double * ctrl_x_tri6 = new double[6];
	double * ctrl_y_tri6 = new double[6];

	generateTri6(ctrl_x_tri6, ctrl_y_tri6);
	writeCoordinatesToFile(outfile, ctrl_x_tri6, ctrl_y_tri6, 6);

	// Read coordinates of the control points
	readCoordinatesFromFile(outfile, 0, ctrl_x_hex8, ctrl_y_hex8, ctrl_z_hex8, 8);
	readCoordinatesFromFile(outfile, 3, ctrl_x_quad4_3d, ctrl_y_quad4_3d, ctrl_z_quad4_3d, 4);
	readCoordinatesFromFile(outfile, 6, ctrl_x_hex27, ctrl_y_hex27, ctrl_z_hex27, 27);
	readCoordinatesFromFile(outfile, 9, ctrl_x_quad9_3d, ctrl_y_quad9_3d, ctrl_z_quad9_3d, 9);
	readCoordinatesFromFile(outfile, 12, ctrl_x_tet4, ctrl_y_tet4, ctrl_z_tet4, 4);
	readCoordinatesFromFile(outfile, 15, ctrl_x_tri3_3d, ctrl_y_tri3_3d, ctrl_z_tri3_3d, 3);
	readCoordinatesFromFile(outfile, 18, ctrl_x_tet10, ctrl_y_tet10, ctrl_z_tet10, 10);
	readCoordinatesFromFile(outfile, 21, ctrl_x_tri6_3d, ctrl_y_tri6_3d, ctrl_z_tri6_3d, 6);
	readCoordinatesFromFile(outfile, 24, ctrl_x_quad4, ctrl_y_quad4, 4);
	readCoordinatesFromFile(outfile, 26, ctrl_x_quad9, ctrl_y_quad9, 9);
	readCoordinatesFromFile(outfile, 28, ctrl_x_tri3, ctrl_y_tri3, 3);
	readCoordinatesFromFile(outfile, 30, ctrl_x_tri6, ctrl_y_tri6, 6);

	std::cout << "hex8&quad4_der0" << std::endl;
	FEAElement * hex8 = new FEAElement_Hex8(3*3*3);
	FEAElement * quad4_der0 = new FEAElement_Quad4_3D_der0(2*2);
	testMemberFunction( hex8, quad4_der0, ctrl_x_hex8, ctrl_y_hex8, ctrl_z_hex8,
		ctrl_x_quad4_3d, ctrl_y_quad4_3d, ctrl_z_quad4_3d);
	
	std::cout << "hex27&quad9_der0" << std::endl;
	FEAElement * hex27 = new FEAElement_Hex27(4*4*4);
	FEAElement * quad9_der0 = new FEAElement_Quad9_3D_der0(3*3);
	testMemberFunction( hex27, quad9_der0, ctrl_x_hex27, ctrl_y_hex27, ctrl_z_hex27,
		ctrl_x_quad9_3d, ctrl_y_quad9_3d, ctrl_z_quad9_3d);
	
	std::cout << "tet4&tri3_der0" << std::endl;
	FEAElement * tet4 = new FEAElement_Tet4(5);
	FEAElement * tri3_der0 = new FEAElement_Triangle3_3D_der0(4);
	testMemberFunction( tet4, tri3_der0, ctrl_x_tet4, ctrl_y_tet4, ctrl_z_tet4,
		ctrl_x_tri3_3d, ctrl_y_tri3_3d, ctrl_z_tri3_3d);
	
	std::cout << "tet10&tri6_der0" << std::endl;
	FEAElement * tet10 = new FEAElement_Tet10(29);
	FEAElement * tri6_der0 = new FEAElement_Triangle6_3D_der0(13);
	testMemberFunction( tet10, tri6_der0, ctrl_x_tet10, ctrl_y_tet10, ctrl_z_tet10,
		ctrl_x_tri6_3d, ctrl_y_tri6_3d, ctrl_z_tri6_3d);
	
	std::cout << "quad4" << std::endl;
	FEAElement * quad4 = new FEAElement_Quad4(2*2);
	testMemberFunction( quad4, ctrl_x_quad4, ctrl_y_quad4);

	std::cout << "quad9" << std::endl;
	FEAElement * quad9 = new FEAElement_Quad9(3*3);
	testMemberFunction( quad9, ctrl_x_quad9, ctrl_y_quad9);

	std::cout << "tri3" << std::endl;
	FEAElement * tri3 = new FEAElement_Triangle3(4);
	testMemberFunction( tri3, ctrl_x_tri3, ctrl_y_tri3);

	std::cout << "tri6" << std::endl;
	FEAElement * tri6 = new FEAElement_Triangle6(13);
	testMemberFunction( tri6, ctrl_x_tri6, ctrl_y_tri6);

	delete hex8; hex8 = nullptr;
	delete quad4_der0; quad4_der0 = nullptr;

	delete hex27; hex27 = nullptr;
	delete quad9_der0; quad9_der0 = nullptr;

	delete tet4; tet4 = nullptr;
	delete tri3_der0; tri3_der0 = nullptr;

	delete tet10; tet10 = nullptr;
	delete tri6_der0; tri6_der0 = nullptr;

	delete quad4; quad4 = nullptr;
	delete quad9; quad9 = nullptr;
	delete tri3; tri3 = nullptr;
	delete tri6; tri6 = nullptr;

	delete [] ctrl_x_hex8; ctrl_x_hex8 = nullptr;
	delete [] ctrl_y_hex8; ctrl_y_hex8 = nullptr;
	delete [] ctrl_z_hex8; ctrl_z_hex8 = nullptr;
	delete [] ctrl_x_quad4_3d; ctrl_x_quad4_3d = nullptr;
	delete [] ctrl_y_quad4_3d; ctrl_y_quad4_3d = nullptr;
	delete [] ctrl_z_quad4_3d; ctrl_z_quad4_3d = nullptr;

	delete [] ctrl_x_hex27; ctrl_x_hex27 = nullptr;
	delete [] ctrl_y_hex27; ctrl_y_hex27 = nullptr;
	delete [] ctrl_z_hex27; ctrl_z_hex27 = nullptr;
	delete [] ctrl_x_quad9_3d; ctrl_x_quad9_3d = nullptr;
	delete [] ctrl_y_quad9_3d; ctrl_y_quad9_3d = nullptr;
	delete [] ctrl_z_quad9_3d; ctrl_z_quad9_3d = nullptr;

	delete [] ctrl_x_tet4; ctrl_x_tet4 = nullptr;
	delete [] ctrl_y_tet4; ctrl_y_tet4 = nullptr;
	delete [] ctrl_z_tet4; ctrl_z_tet4 = nullptr;
	delete [] ctrl_x_tri3_3d; ctrl_x_tri3_3d = nullptr; 
	delete [] ctrl_y_tri3_3d; ctrl_y_tri3_3d = nullptr; 
	delete [] ctrl_z_tri3_3d; ctrl_z_tri3_3d = nullptr;

	delete [] ctrl_x_tet10; ctrl_x_tet10 = nullptr;
	delete [] ctrl_y_tet10; ctrl_y_tet10 = nullptr;
	delete [] ctrl_z_tet10; ctrl_z_tet10 = nullptr;
	delete [] ctrl_x_tri6_3d; ctrl_x_tri6_3d = nullptr;
	delete [] ctrl_y_tri6_3d; ctrl_y_tri6_3d = nullptr;
	delete [] ctrl_z_tri6_3d; ctrl_z_tri6_3d = nullptr;

	delete [] ctrl_x_quad4; ctrl_x_quad4 = nullptr; 
	delete [] ctrl_y_quad4; ctrl_y_quad4 = nullptr;

	delete [] ctrl_x_quad9; ctrl_x_quad9 = nullptr;
	delete [] ctrl_y_quad9; ctrl_y_quad9 = nullptr;

	delete [] ctrl_x_tri3; ctrl_x_tri3 = nullptr; 
	delete [] ctrl_y_tri3; ctrl_y_tri3 = nullptr;

	delete [] ctrl_x_tri6; ctrl_x_tri6 = nullptr;
	delete [] ctrl_y_tri6; ctrl_y_tri6 = nullptr;

	delete [] outfile; outfile = nullptr;

	return 0;
}
// EOF
