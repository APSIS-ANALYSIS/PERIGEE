// Test the new TET_T functions

#include "Tet_Tools.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
  vector<double> ppt;
  vector<int> ien;
  int numpt, numcl;
  vector<int> eindex, e2index, nindex, phytag;

  ppt.clear(); ien.clear();
  numpt = 4;
  numcl = 1;
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(1.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(1.0); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(3.0);
  //ppt.push_back(0.8); ppt.push_back(0.5); ppt.push_back(0.0);
  //ppt.push_back(0.1); ppt.push_back(0.5); ppt.push_back(0.0);

  ien.push_back(0); ien.push_back(1); ien.push_back(2);
  ien.push_back(3); //ien.push_back(4); ien.push_back(5);

  nindex.clear(); eindex.clear(); phytag.clear();
  nindex.push_back(1);
  nindex.push_back(2);
  nindex.push_back(3);
  nindex.push_back(4);
  //nindex.push_back(5);
  //nindex.push_back(-3);

  eindex.push_back(2);
  //e2index.push_back(-3);
  
  phytag.push_back(1);

  // TET_T::write_tet_grid("sample_tet", numpt, numcl, ppt, ien);
  // TET_T::write_tet_grid_node_elem_index("sample_tet_idx", numpt, numcl, ppt, ien);
  // TET_T::write_tet_grid_node_elem_index("sample_tet_idx_2", numpt, numcl, ppt, ien, nindex, eindex);
  // TET_T::write_tet_grid("sample_tet_phytag", numpt, numcl, ppt, ien, phytag, true, 1);
  TET_T::write_tet_grid("sample_tet_phytag_2", true, numpt, numcl, ppt, ien, phytag, nindex, eindex);

  return 0;
}

//EOF
