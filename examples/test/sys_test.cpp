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
#include "ElemBC_3D_wall_turbulence.hpp"
#include "EBC_Partition_outflow.hpp"
#include "EBC_Partition_wall_turbulence.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_EBC.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"

int main(int argc, char *argv[])
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const std::string part_file("part");
  int elemType = 501;
  SYS_T::GetOptionInt("-elem_type", elemType);

  // // -------------------------------------------------------------------------------------
  // // test preprocess
  // constexpr int dofNum = 4; // degree-of-freedom for the physical problem
  // constexpr int dofMat = 4; // degree-of-freedom in the matrix problem

  // int weakBC_type = 1;

  // int cpu_size = 1;
  // int in_ncommon = 2;
  // bool isDualGraph = true;

  // std::string geo_file("./whole_vol.vtu");
  // std::string test_sur("./testface_vol");

  // if(elemType == 501 || elemType == 601)
  //   test_sur += ".vtp";
  // else if(elemType == 502 || elemType == 602)
  //   test_sur += ".vtu";
  // else
  //   SYS_T::print_fatal("Wrong -ele_type\n");

  // // volumetric mesh
  // int nFunc, nElem;
  // std::vector<int> vecIEN;
  // std::vector<double> ctrlPts;
  
  // VTK_T::read_vtu_grid(geo_file, nFunc, nElem, ctrlPts, vecIEN);
  
  // IIEN * IEN = new IEN_FEM(nElem, vecIEN);
  // VEC_T::clean( vecIEN );

  //   IMesh * mesh = nullptr;

  // switch( elemType )
  // {
  //   case 501:
  //     mesh = new Mesh_Tet(nFunc, nElem, 1);
  //     break;
  //   case 502:
  //     mesh = new Mesh_Tet(nFunc, nElem, 2);
  //     break;
  //   // case 601:
  //   //   mesh = new Mesh_FEM(nFunc, nElem, 8, 1);
  //   //   break;
  //   // case 602:
  //   //   mesh = new Mesh_FEM(nFunc, nElem, 27, 2);
  //   //   break;
  //   default:
  //     SYS_T::print_fatal("Error: elemType %d is not supported.\n", elemType);
  //     break;
  // }
  
  // SYS_T::print_fatal_if( IEN->get_nLocBas() != mesh->get_nLocBas(), "Error: the nLocBas from the Mesh %d and the IEN %d classes do not match. \n", mesh->get_nLocBas(), IEN->get_nLocBas());

  // mesh -> print_info();

  // std::vector<std::string> dir_list {}; // empty list
  // std::vector<INodalBC *> NBC_list( dofMat, nullptr );
  // NBC_list[0] = new NodalBC( nFunc );
  // NBC_list[1] = new NodalBC( dir_list, nFunc );
  // NBC_list[2] = new NodalBC( dir_list, nFunc );
  // NBC_list[3] = new NodalBC( dir_list, nFunc );

  // std::vector< std::string > sur_file_out {test_sur};

  // // previous
  // std::vector< Vector_3 > outlet_outvec( sur_file_out.size() );

  // if(elemType == 501 || elemType == 502)
  // {
  //     for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
  //   outlet_outvec[ii] = TET_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN );
  // }
  // // else if(elemType == 601 || elemType == 602)
  // // {
  // //     for(unsigned int ii=0; ii<sur_file_out.size(); ++ii)
  // //   outlet_outvec[ii] = HEX_T::get_out_normal( sur_file_out[ii], ctrlPts, IEN );
  // // }
  // else
  //   SYS_T::print_fatal("Wrong -ele_type\n");

  // ElemBC * ebc = new ElemBC_3D_outflow( sur_file_out, outlet_outvec, elemType );

  // ebc -> resetSurIEN_outwardnormal( IEN );

  // // new
  // ElemBC * wbc = new ElemBC_3D_wall_turbulence( sur_file_out, weakBC_type, IEN, elemType );

  // // partition
  // IGlobal_Part * global_part = new Global_Part_Serial( mesh, "epart", "npart" );
  // // new Global_Part_METIS( cpu_size, in_ncommon,
  // //   isDualGraph, mesh, IEN, "epart", "npart" );

  // Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  //   mnindex->write_hdf5("node_mapping");

  // SYS_T::Timer * mytimer = new SYS_T::Timer();
  // for(int proc_rank = 0; proc_rank < cpu_size; ++proc_rank)
  // {
  //   mytimer->Reset();
  //   mytimer->Start();
  //   IPart * part = new Part_FEM( mesh, global_part, mnindex, IEN,
  //       ctrlPts, proc_rank, cpu_size, dofNum, dofMat, elemType );
    
  //   mytimer->Stop();
  //   cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

  //   part -> write( part_file );
  //   part -> print_part_loadbalance_edgecut();

  //   // Partition Elemental BC and write to h5 file
  //   EBC_Partition * ebcpart = new EBC_Partition_outflow(part, mnindex, ebc, NBC_list);
  //   ebcpart -> write_hdf5( part_file );

  //   // Partition Weak BC and write to h5 file
  //   EBC_Partition * wbcpart = new EBC_Partition_wall_turbulence(part, mnindex, wbc);
  //   wbcpart -> write_hdf5( part_file );

  //   delete wbcpart; delete ebcpart; delete part;
  // }

  // // Finalize the test preprocess
  // delete mytimer;
  // delete ebc; delete wbc;
  // for(auto it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;
  // delete mnindex;
  // delete global_part;
  // delete mesh; delete IEN;

  // -------------------------------------------------------------------------------------
  // test analysis
  const PetscMPIInt rank = 0; // SYS_T::get_MPI_rank();

  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;
  IQuadPts * quads = nullptr;
  switch(elemType)
  {
    case 501:
      elementv = new FEAElement_Tet4( 4 ); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elements = new FEAElement_Triangle3_3D_der0( 4 );
      quads = new QuadPts_Gauss_Triangle( 4 );
      break;
    case 502:
      elementv = new FEAElement_Tet10_v2( 13 ); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elements = new FEAElement_Triangle6_3D_der0( 13 );
      quads = new QuadPts_Gauss_Triangle( 13 );
      break;
    // case 601:
    //   elementv = new FEAElement_Hex8( 2 ); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //   elements = new FEAElement_Quad4_3D_der0( 2 );
    //   quads = new QuadPts_Gauss_Quad( 2 );
    //   break;
    // case 602:
    //   elementv = new FEAElement_Hex27( 4 ); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //   elements = new FEAElement_Quad9_3D_der0( 4 );
    //   quads = new QuadPts_Gauss_Quad( 4 );
    //   break;
    default:
      SYS_T::print_fatal("Error: elemType %d is not supported.\n", elemType);
      break;
  }

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  // previous
  ALocal_EBC * locebc = new ALocal_EBC_outflow(part_file, rank);

  // new
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank); // volumetric IEN
  FEANode * fNode = new FEANode(part_file, rank); // volumetric nodes
  ALocal_WeakBC * locwbc = new ALocal_WeakBC(part_file, rank, 1.0);

  int s_nElem = locwbc -> get_num_ele();
  int s_nLocBas = locebc -> get_cell_nLocBas(0);
  int v_nLocBas = GMIptr -> get_nLocBas();

  // mapped node id in each element
  std::vector<int> s_node_global (s_nLocBas, -1), s_node_local (s_nLocBas, -1);
  
  std::vector<int> v_node_global (v_nLocBas, -1);

  // intNA in EBC_outflow
  std::vector<double> s_intNA = locebc -> get_intNA(0);

  // node xyz in each volume element
  double vctrl_x[v_nLocBas], vctrl_y[v_nLocBas], vctrl_z[v_nLocBas];

  // new surface integral of each NA
  std::vector<double> new_intNA (GMIptr->get_nFunc(), 0.0);
  for(int ee{0}; ee < s_nElem; ++ee)
  {
    int vol_ee = locwbc -> get_part_vol_ele_id(ee);
    for(int jj{0}; jj < v_nLocBas; ++jj)
      v_node_global[jj] = locIEN -> get_LIEN(vol_ee, jj);

    fNode -> get_ctrlPts_xyz(v_nLocBas, &v_node_global[0], vctrl_x, vctrl_y, vctrl_z);

    int face_id = locwbc -> get_ele_face_id(ee);
    elementv -> buildBasis(quads, face_id, vctrl_x, vctrl_y, vctrl_z);

    for(int qua{0}; qua < quads->get_num_quadPts(); ++qua)
    {
      std::vector<double> NA = elementv -> get_R(qua);

      double surface_area {0.0};
      Vector_3 n_out = elementv -> get_2d_normal_out(qua, surface_area);
      double gwts = surface_area * quads->get_qw(qua);

      for(int jj{0}; jj < v_nLocBas; ++jj)
        new_intNA[v_node_global[jj]] += gwts * NA[jj];
    }
  }

  // compare by element
  std::vector<double> old_int (s_nLocBas, 0.0), new_int (v_nLocBas, 0.0);

  for(int ee{0}; ee < s_nElem; ++ee)
  {
    // previous
    for(int ii{0}; ii < s_nLocBas; ++ii)
    {
      s_node_local[ii] = locebc -> get_local_cell_ien(0, s_nLocBas * ee + ii);
      s_node_global[ii] = locebc -> get_local_cell_node_pos(0, s_node_local[ii]);

      old_int[ii] = s_intNA[s_node_local[ii]];
    }

    // new
    int vol_ee = locwbc -> get_part_vol_ele_id(ee);
    for(int jj{0}; jj < v_nLocBas; ++jj)
    {
      v_node_global[jj] = locIEN -> get_LIEN(vol_ee, jj);

      new_int[jj] = new_intNA[v_node_global[jj]];
    }

    SYS_T::commPrint("ee = %d \n", ee);
    SYS_T::commPrint("  Old:\n");
    SYS_T::commPrint("    Node:\n\t");
    VEC_T::print(s_node_global);
    SYS_T::commPrint("    intNA:\n\t");
    VEC_T::print(old_int);

    SYS_T::commPrint("  New:\n");
    SYS_T::commPrint("    Node:\n\t");
    VEC_T::print(v_node_global);
    SYS_T::commPrint("    intNA:\n\t");
    VEC_T::print(new_int);
    SYS_T::commPrint("\n");
  }

  delete quads;
  delete elements; delete elementv;

  // -------------------------------------------------------------------------------------

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
