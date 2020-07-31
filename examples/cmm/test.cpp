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
  ppt.push_back(0.5); ppt.push_back(0.0); ppt.push_back(0.0);
  ppt.push_back(0.8); ppt.push_back(0.5); ppt.push_back(0.0);
  ppt.push_back(0.1); ppt.push_back(0.5); ppt.push_back(0.0);

  ien.push_back(0); ien.push_back(1); ien.push_back(2);
  ien.push_back(3); ien.push_back(4); ien.push_back(5);

  nindex.push_back(1);
  nindex.push_back(2);
  nindex.push_back(3);
  nindex.push_back(4);
  nindex.push_back(5);
  nindex.push_back(-3);

  eindex.push_back(2);

  TET_T::write_quadratic_triangle_grid("sample_tri", numpt, numcl, ppt, ien, nindex, eindex);

  return 0;
}

//EOF
