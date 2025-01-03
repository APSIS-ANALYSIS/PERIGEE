// ============================================================================
// driver.cpp
//
// Finite element code for stress recovery.
//
// Date: Dec. 20 2023
// ============================================================================
#include "AGlobal_Mesh_Info.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_Elem.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "PLocAssem_Stress_Recovery.hpp"
#include "PGAssem_Stress_Recovery.hpp"
#include "PLinear_Solver_PETSc.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points
  int nqp_vol = 5;
  int nqp_vol_1D = 4;

  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  // Part file location
  std::string part_file("./apart/part");

  int time_start = 0, time_step = 1, time_end = 1;
  
  // Input and output files' names
  std::string isol_bname = "SOL_disp_";
  std::string osol_bname = "SOL_Cauchy_";

  // Yaml options
  bool is_loadYaml = true;
  std::string yaml_file("./smooth.yml");

#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_perigee_art();

  // Yaml arguments
  SYS_T::GetOptionBool("-is_loadYaml", is_loadYaml);
  SYS_T::GetOptionString("-yaml_file", yaml_file);

  if(is_loadYaml) SYS_T::InsertFileYAML( yaml_file,  false );
  
  // Read command line arguments
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-nqp_vol", nqp_vol);
  SYS_T::GetOptionInt("-nqp_vol_1d", nqp_vol_1D);
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionString("-isol_bname", isol_bname);
  SYS_T::GetOptionString("-osol_bname", osol_bname);

  // Print arguments
  SYS_T::cmdPrint("-nqp_vol:", nqp_vol);
  SYS_T::cmdPrint("-nqp_vol_1d", nqp_vol_1D);
  SYS_T::cmdPrint("-nz_estimate:", nz_estimate);
  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-time_start", time_start);
  SYS_T::cmdPrint("-time_step", time_step);
  SYS_T::cmdPrint("-time_end", time_end);
  SYS_T::cmdPrint("-isol_bname", isol_bname);
  SYS_T::cmdPrint("-osol_bname", osol_bname);

  MPI_Barrier(PETSC_COMM_WORLD);

  FEANode * fNode = new FEANode(part_file, rank);

  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  AGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info(part_file, rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  APart_Node * pNode = new APart_Node(part_file, rank);

  SYS_T::commPrint("===> Data from HDF5 files are read from disk.\n");

  SYS_T::print_fatal_if( size != PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");
  
  SYS_T::commPrint("===> %d processor(s) are assigned for solution smoother.\n", size);

  // Finite element container & Quadrature rules
  SYS_T::commPrint("===> Setup element container. \n");
  SYS_T::commPrint("===> Build quadrature rules. \n");
  FEAElement * elementv = nullptr;
  IQuadPts * quadv = nullptr;

  if( GMIptr->get_elemType() == FEType::Tet4 )
  {
    if( nqp_vol > 5 ) SYS_T::commPrint("Warning: the tet element is linear and you are using more than 5 quadrature points.\n");
    
    elementv = new FEAElement_Tet4( nqp_vol ); // elem type Tet4
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
  }
  else if( GMIptr->get_elemType() == FEType::Tet10 )
  {
    SYS_T::print_fatal_if( nqp_vol < 29, "Error: not enough quadrature points for tets.\n" );

    elementv = new FEAElement_Tet10( nqp_vol ); // elem type Tet10
    quadv = new QuadPts_Gauss_Tet( nqp_vol );
  }
  else if( GMIptr->get_elemType() == FEType::Hex8 )
  {
    SYS_T::print_fatal_if( nqp_vol_1D < 2, "Error: not enough quadrature points for hex.\n" );
  
    elementv = new FEAElement_Hex8( nqp_vol_1D * nqp_vol_1D * nqp_vol_1D ); // elem type Hex8
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
  }
  else if( GMIptr->get_elemType() == FEType::Hex27 )
  {
    SYS_T::print_fatal_if( nqp_vol_1D < 4, "Error: not enough quadrature points for hex.\n" );
  
    elementv = new FEAElement_Hex27( nqp_vol_1D * nqp_vol_1D * nqp_vol_1D ); // elem type Hex27
    quadv = new QuadPts_Gauss_Hex( nqp_vol_1D );
  }
  else SYS_T::print_fatal("Error: Element type not supported.\n");

  // print the information of element and quadrature rule
  elementv->print_info();
  quadv->print_info();

  // Solution vector
  PDNSolution * disp = new PDNSolution( pNode, 3 );
  PDNSolution * stress = new PDNSolution( pNode, 6 );

  disp -> ScaleValue( 0.0 );
  stress -> ScaleValue( 0.0 );

  // Material model
  const double in_module = 1.0e+5;
  const double in_nu = 0.3;
  
  // Local assembly routine
  IPLocAssem * locAssem_ptr = new PLocAssem_Stress_Recovery(elementv->get_nLocBas(), in_nu, in_module);

  // Global assembly
  SYS_T::commPrint("===> Initializing Mat K and Vec G ... \n");
  IPGAssem * gloAssem_ptr = new PGAssem_Stress_Recovery( locAssem_ptr,
      GMIptr, locElem, locIEN, pNode, nz_estimate );
  
  SYS_T::commPrint("===> Assembly nonzero estimate matrix ... \n");
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, locIEN, pNode );
  // MatView(gloAssem_ptr->K, PETSC_VIEWER_STDOUT_WORLD);

  SYS_T::commPrint("===> Matrix nonzero structure fixed. \n");
  gloAssem_ptr->Fix_nonzero_err_str();
  gloAssem_ptr->Clear_KG();

  // Linear solver
  PLinear_Solver_PETSc * lsolver_ptr = new PLinear_Solver_PETSc();

  // Smooth the solutions
  // Assemble mass matrix
  gloAssem_ptr->Assem_mass_residual(disp, locElem, locAssem_ptr, elementv, quadv, locIEN, fNode, pNode);

  lsolver_ptr->SetOperator( gloAssem_ptr->K );

  // Smooth 
  for(int time = time_start; time<=time_end; time+=time_step)
  {
    std::ostringstream time_index;
    std::string name_to_read(isol_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_read.append(time_index.str());

    std::string name_to_write(osol_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %d: Read %s and Write %s \n",
        time, name_to_read.c_str(), name_to_write.c_str() );

    disp->ReadBinary(name_to_read);

    gloAssem_ptr->Clear_G();

    gloAssem_ptr->Assem_residual(disp, locElem, locAssem_ptr, elementv, quadv, locIEN, fNode, pNode);

    lsolver_ptr->Solve( gloAssem_ptr->G, stress );

    stress->WriteBinary(name_to_write);

    SYS_T::commPrint("\n" );
  }
  // Print complete solver info
  lsolver_ptr -> print_info();

  delete lsolver_ptr; delete gloAssem_ptr; delete locAssem_ptr; delete disp; delete stress;
  delete fNode; delete locIEN; delete GMIptr; delete locElem; delete pNode; delete PartBasic;
  delete quadv; delete elementv;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
