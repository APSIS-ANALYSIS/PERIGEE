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
  std::string geo_file("./whole_vol.vtu");
  
  //read sold data
  std::vector<int> global_node_index = TET_T::read_int_PointData(geo_s_file, "GlobalNodeID");
  std::vector<double> radial_normal_vec = TET_T::read_double_vec_3_PointData(geo_s_file, "radial_normal");
  std::vector<double> longitudinal_normal_vec = TET_T::read_double_vec_3_PointData(geo_s_file, "longitudinal_normal");
  std::vector<double> circumferential_normal_vec = TET_T::read_double_vec_3_PointData(geo_s_file, "circumferential_normal");

  //read whole geo data
  int nFunc, nElem;
  std::vector<int> vecIEN, phy_tag;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid( geo_file, nFunc, nElem, ctrlPts, vecIEN, phy_tag );
  
  for(unsigned int ii=0; ii<phy_tag.size(); ++ii)
  {
    if(phy_tag[ii] != 0 && phy_tag[ii] != 1) SYS_T::print_fatal("Error: FSI problem, the physical tag for element should be 0 (fluid domain) or 1 (solid domain).\n");
  }

  // Generate IEN
  IIEN * IEN_v = new IEN_Tetra_P1( nElem, vecIEN );

  // Generate the list of nodes for fluid and solid
  std::vector<int> v_node_f, v_node_s; v_node_f.clear(); v_node_s.clear();

  for(int ee=0; ee<nElem; ++ee)
  {
    if( phy_tag[ee] == 0 )
    {
      for(int ii=0; ii<4; ++ii) v_node_f.push_back( IEN_v->get_IEN(ee, ii) );
    }
    else
    {
      for(int ii=0; ii<4; ++ii) v_node_s.push_back( IEN_v->get_IEN(ee, ii) );
    }
  }

  VEC_T::sort_unique_resize( v_node_f ); VEC_T::sort_unique_resize( v_node_s );

  //compare global_node_index vs. v_node_s
  if (v_node_s == global_node_index) cout<<"global_node_index is equal to v_node_s"<<endl;

  //for(auto it = radial_normal_vec.begin(); it != radial_normal_vec.end(); it++) cout<<*it<<endl;

  cout<<"v_node_s size: "<<v_node_s.size()<<endl;
  cout<<"global_node_index size: "<<global_node_index.size()<<endl;
  cout<<"radial_normal_vec size: "<<radial_normal_vec.size()<<endl;
  cout<<"longitudinal_normal_vec size: "<<longitudinal_normal_vec.size()<<endl;
  cout<<"circumferential_normal_vec size: "<<circumferential_normal_vec.size()<<endl;

  /*
  for(unsigned int ii=0; ii<global_node_index.size(); ii++)
  {
    cout<<"nodeID: "<<global_node_index[ii]<<endl;
    cout<<"radial normal: "<<radial_normal_vec[ii*3]<<"\t"<<radial_normal_vec[ii*3+1]<<"\t"<<radial_normal_vec[ii*3+2]<<endl;
    cout<<"longitudinal normal: "<<longitudinal_normal_vec[ii*3]<<"\t"<<longitudinal_normal_vec[ii*3+1]<<"\t"<<longitudinal_normal_vec[ii*3+2]<<endl;
    cout<<"circumferential normal: "<<circumferential_normal_vec[ii*3]<<"\t"<<circumferential_normal_vec[ii*3+1]<<"\t"<<circumferential_normal_vec[ii*3+2]<<endl;
  }
  */

  return EXIT_SUCCESS;
}

// EOF
