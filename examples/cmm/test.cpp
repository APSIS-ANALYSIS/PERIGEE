// Test the new TET_T functions

#include "Tet_Tools.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
  vector<double> ppt;
  vector<int> ien;
  int numpt, numcl;
  vector<int> eindex, nindex;

  ppt.clear(); ien.clear();
  numpt = 6;
  numcl = 1;
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(1.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(1.0); ppt.push_back(0.0);
  //ppt.push_back(0.5); ppt.push_back(0.0); ppt.push_back(0.0);
  //ppt.push_back(0.8); ppt.push_back(0.5); ppt.push_back(0.0);
  //ppt.push_back(0.1); ppt.push_back(0.5); ppt.push_back(0.0);

  ien.push_back(0); ien.push_back(1); ien.push_back(2);
  //ien.push_back(3); ien.push_back(4); ien.push_back(5);

  nindex.push_back(1);
  nindex.push_back(2);
  nindex.push_back(3);
  //nindex.push_back(4);
  //nindex.push_back(5);
  //nindex.push_back(-3);

  eindex.push_back(2);

  //TET_T::write_quadratic_triangle_grid("sample_tri", numpt, numcl, ppt, ien, nindex, eindex);

  vtkPolyData * grid_w = vtkPolyData::New();

  TET_T::gen_triangle_grid( grid_w, 3, 1, ppt, ien );
  TET_T::add_int_PointData( grid_w, nindex, "gnid");
  TET_T::add_int_CellData( grid_w, eindex, "geid");

  std::vector<double> n2idx;
  n2idx.push_back(2.5);
  n2idx.push_back(-1.5);
  n2idx.push_back(32.3333);
  
  TET_T::add_double_PointData( grid_w, n2idx, "testid" );

  TET_T::write_vtkXMLPolyData("sample_tri", grid_w);

  grid_w -> Delete();

  return 0;
}

//EOF
