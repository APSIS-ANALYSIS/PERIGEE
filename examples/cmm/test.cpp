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
  numpt = 10;
  numcl = 1;
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(1.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(1.0); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(1.0);
  ppt.push_back(0.5); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(0.5); ppt.push_back(0.5); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(0.5); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(0.5);
  ppt.push_back(0.5); ppt.push_back(0.3); ppt.push_back(0.5);
  ppt.push_back(0.1); ppt.push_back(0.5); ppt.push_back(0.5);

  ien.push_back(0); ien.push_back(1); ien.push_back(2);
  ien.push_back(3); ien.push_back(4); ien.push_back(5);
  ien.push_back(6); ien.push_back(7); ien.push_back(8);
  ien.push_back(9);

  nindex.clear(); eindex.clear(); phytag.clear();
  nindex.push_back(1);
  nindex.push_back(2);
  nindex.push_back(3);
  nindex.push_back(4);
  nindex.push_back(5);
  nindex.push_back(-3);
  nindex.push_back(-2);
  nindex.push_back(-1);
  nindex.push_back(-5);
  nindex.push_back(5);

  eindex.push_back(2);
  
  phytag.push_back(1);

  TET_T::write_tet_grid("sample_tet", numpt, numcl, ppt, ien);
  TET_T::write_tet_grid("sample_tet_phytag", numpt, numcl, ppt, ien, phytag, true, 1);
  TET_T::write_tet_grid("sample_tet_phytag_2", numpt, numcl, ppt, ien, phytag, false, 22 );
  TET_T::write_tet_grid("sample_tet_phytag_3", numpt, numcl, ppt, ien, nindex, eindex, phytag, false );
  TET_T::write_tet_grid("sample_tet_phytag_4", numpt, numcl, ppt, ien, nindex, eindex, phytag, true );

  ppt.clear();
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(1.1); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(0.5); ppt.push_back(3.0); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(1.0);

  TET_T::Tet4 * a = new TET_T::Tet4(ppt);
  a -> print_info();
  a -> write_vtu("hello");
  delete a;
  return 0;
}

//EOF
