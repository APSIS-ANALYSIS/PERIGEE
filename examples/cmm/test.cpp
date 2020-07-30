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
  numpt = 4;
  numcl = 2;
  ppt.push_back(0.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(15.0); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(0.0); ppt.push_back(1.0); ppt.push_back(0.0);
  ppt.push_back(-1.0); ppt.push_back(0.0); ppt.push_back(1.0);

  ien.push_back(0); ien.push_back(1); ien.push_back(2);
  ien.push_back(0); ien.push_back(3); ien.push_back(2);

  nindex.push_back(1);
  nindex.push_back(2);
  nindex.push_back(3);
  nindex.push_back(4);

  eindex.push_back(2);
  eindex.push_back(1);

  TET_T::write_triangle_grid("sample_tri", numpt, numcl, ppt, ien, nindex, eindex);

  return 0;
}

//EOF
