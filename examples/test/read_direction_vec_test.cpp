#include <unistd.h>
#include <vector>
#include "Sys_Tools.hpp"
#include "Tet_Tools.hpp"
#include "Vec_Tools.hpp"
#include "IEN_Tetra_P1.hpp"

int main(int argc, char *argv[])
{
  // Input file
  std::string geo_s_file("./outsolid.vtu");
    
  //read solid node data
  std::vector<int> global_node_index = TET_T::read_int_PointData(geo_s_file, "GlobalNodeID");
  std::vector<double> radial_normal_vec = TET_T::read_double_vec_3_PointData(geo_s_file, "radial_normal");
  std::vector<double> longitudinal_normal_vec = TET_T::read_double_vec_3_PointData(geo_s_file, "longitudinal_normal");
  std::vector<double> circumferential_normal_vec = TET_T::read_double_vec_3_PointData(geo_s_file, "circumferential_normal");

  cout<<"global_node_index size: "<<global_node_index.size()<<endl;
  cout<<"radial_normal_vec size: "<<radial_normal_vec.size()<<endl;
  cout<<"longitudinal_normal_vec size: "<<longitudinal_normal_vec.size()<<endl;
  cout<<"circumferential_normal_vec size: "<<circumferential_normal_vec.size()<<endl;

  for(unsigned int ii=0; ii<global_node_index.size(); ii++)
  {
    cout<<"nodeID: "<<global_node_index[ii]<<endl;
    cout<<"radial normal: "<<radial_normal_vec[ii*3]<<"\t"<<radial_normal_vec[ii*3+1]<<"\t"<<radial_normal_vec[ii*3+2]<<endl;
    cout<<"longitudinal normal: "<<longitudinal_normal_vec[ii*3]<<"\t"<<longitudinal_normal_vec[ii*3+1]<<"\t"<<longitudinal_normal_vec[ii*3+2]<<endl;
    cout<<"circumferential normal: "<<circumferential_normal_vec[ii*3]<<"\t"<<circumferential_normal_vec[ii*3+1]<<"\t"<<circumferential_normal_vec[ii*3+2]<<endl;
  }
  
  return EXIT_SUCCESS;
}

// EOF
