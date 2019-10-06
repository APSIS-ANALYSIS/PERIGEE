#include "AExtractor_2D_NURBS_xy.hpp"

AExtractor_2D_NURBS_xy::AExtractor_2D_NURBS_xy(const std::string &part_file, const int &rank)
{
  HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);
  
  h5reader->get_EXT_x(extractor_x);
  h5reader->get_EXT_y(extractor_y);

  h5reader->get_GMI_degree(sdegree, tdegree);

  delete h5reader;
}


AExtractor_2D_NURBS_xy::AExtractor_2D_NURBS_xy( const HDF5_PartReader * const &h5reader )
{
  h5reader->get_EXT_x(extractor_x);
  h5reader->get_EXT_y(extractor_y);

  h5reader->get_GMI_degree(sdegree, tdegree);
}


AExtractor_2D_NURBS_xy::~AExtractor_2D_NURBS_xy()
{}


void AExtractor_2D_NURBS_xy::get_EXT_x(const int &e, std::vector<double> &ext_x) const
{
  ext_x.clear();
  int ssize = (sdegree+1) * (sdegree+1);
  ext_x.resize(ssize);

  for(int ii=0; ii<ssize; ++ii)
    ext_x[ii] = extractor_x[ssize*e+ii];

  VEC_T::shrink2fit(ext_x);
}


void AExtractor_2D_NURBS_xy::get_EXT_x(const int &e, double * &ext_x) const
{
  int ssize = (sdegree+1) * (sdegree+1);
  ext_x = new double [ssize];
  for(int ii=0; ii<ssize; ++ii)
    ext_x[ii] = extractor_x[ssize*e+ii];
}


void AExtractor_2D_NURBS_xy::get_EXT_y(const int &e, std::vector<double> &ext_y) const
{
  ext_y.clear();
  int tsize = (tdegree+1)*(tdegree+1);
  ext_y.resize(tsize);

  for(int ii=0; ii<tsize; ++ii)
    ext_y[ii] = extractor_y[tsize*e + ii];

  VEC_T::shrink2fit(ext_y);
}

void AExtractor_2D_NURBS_xy::get_EXT_y(const int &e, double * &ext_y) const
{
  int tsize = (tdegree+1) * (tdegree+1);
  ext_y = new double [tsize];
  for(int ii=0; ii<tsize; ++ii)
    ext_y[ii] = extractor_y[tsize*e+ii];
}


void AExtractor_2D_NURBS_xy::get_EXT_z(const int &e, std::vector<double> &ext_z) const
{
  ext_z.clear();
}



void AExtractor_2D_NURBS_xy::get_EXT_z(const int &e, double * &ext_z) const
{
  ext_z = NULL;
}







// EOF
