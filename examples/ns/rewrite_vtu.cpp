#include "VTK_Tools.hpp"
#include "Tet_Tools.hpp"
#include "HDF5_Writer.hpp"

std::string gen_file_name(std::string &me, std::string &co,
  std::string &su, int &id);

int main( int argc, char * argv[] )
{
  #if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif

  // SYS_T::execute("rm -rf inherited.h5");
  int id_start = 34280;
  int id_step = 450;
  int id_end = 43280;
  std::string cond = "c1";
  std::string mesh_used = "22M";

  SYS_T::GetOptionString("-mesh_used", mesh_used);
  SYS_T::GetOptionString("-cond", cond);
  SYS_T::GetOptionInt("-id_start", id_start);
  SYS_T::GetOptionInt("-id_step", id_step);
  SYS_T::GetOptionInt("-id_end", id_end);

  std::string f_suffix = "0_0.vtu";

  std::string empty = "";

  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;

  std::vector<double> elem_vol;

  for(int index = id_start; index <= id_end; index += id_step)
  {
    std::string input_file = gen_file_name(empty, cond, f_suffix, index);

    std::cout<< "Reading file: "<< input_file <<std::endl;
    
    if(index == id_start)
      VTK_T::read_vtu_grid(input_file, nFunc, nElem, ctrlPts, vecIEN, elem_vol);
    else
      VTK_T::read_vtu_grid(input_file, nFunc, ctrlPts);

    // nElem = VEC_T::get_size(elem_vol);
    // std::cout << nElem << " " << VEC_T::get_size(vecIEN) << std::endl;

    std::vector<double> sol_p = VTK_T::read_double_PointData(input_file, "PressureStagnation");
    std::vector<Vector_3> sol_velo = VTK_T::read_Vector3_PointData(input_file, "Velocity");

    std::cout<< "Read Successfully\n";

    // write whole domain
    std::vector<DataVecStr<int>> input_vtk_data {};

    std::vector<DataVecStr<double>> double_vtk_data {};

    std::vector<DataVecStr<Vector_3>> vec_vtk_data {};
    
    std::vector<int> temp_nid(nFunc, 0);
    for(int ii=0; ii<nFunc; ++ii) temp_nid[ii] = ii;
    input_vtk_data.push_back({temp_nid, "GlobalNodeID", AssociateObject::Node});

    std::vector<int> temp_eid(nElem, 0);
    for(int ii=0; ii<nElem; ++ii) temp_eid[ii] = ii;
    input_vtk_data.push_back({temp_eid, "GlobalElementID", AssociateObject::Cell});

    double_vtk_data.push_back({sol_p, "Pressure", AssociateObject::Node});
    double_vtk_data.push_back({elem_vol, "CellVolume", AssociateObject::Cell});

    vec_vtk_data.push_back({sol_velo, "Velocity", AssociateObject::Node});

    std::string output_file_name = gen_file_name(mesh_used, cond, empty, index);

    std::cout<< "Writing file: "<< output_file_name <<std::endl;
    TET_T::write_tet_grid( output_file_name, nFunc, nElem, ctrlPts,
        vecIEN, input_vtk_data, double_vtk_data, vec_vtk_data, 1 );
  }

  std::cout<< "Rewriting finished" <<std::endl;

  PetscFinalize();

  return EXIT_SUCCESS;
}

std::string gen_file_name(std::string &me, std::string &co,
  std::string &su, int &id)
{
  std::ostringstream ss;
  
  if(me != "")
    ss << me << "_";

  ss << co << "_" << id;

  if(su != "")
    ss << "_" << su;

  return ss.str();
}