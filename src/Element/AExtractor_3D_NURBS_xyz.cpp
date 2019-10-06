#include "AExtractor_3D_NURBS_xyz.hpp"

AExtractor_3D_NURBS_xyz::AExtractor_3D_NURBS_xyz(const std::string &part_file, const int &rank)
{
  HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);
  
  h5reader->get_EXT_x(extractor_x);
  h5reader->get_EXT_y(extractor_y);
  h5reader->get_EXT_z(extractor_z);

  h5reader->get_GMI_degree(sdegree, tdegree, udegree);

  delete h5reader;
}

AExtractor_3D_NURBS_xyz::AExtractor_3D_NURBS_xyz( const HDF5_PartReader * const &h5reader )
{
  h5reader->get_EXT_x(extractor_x);
  h5reader->get_EXT_y(extractor_y);
  h5reader->get_EXT_z(extractor_z);

  h5reader->get_GMI_degree(sdegree, tdegree, udegree);
}

AExtractor_3D_NURBS_xyz::~AExtractor_3D_NURBS_xyz()
{}

void AExtractor_3D_NURBS_xyz::get_EXT_x(const int &e, std::vector<double> &ext_x) const
{
  ext_x.clear();
  int ssize = (sdegree+1) * (sdegree+1);
  ext_x.resize(ssize);

  for(int ii=0; ii<ssize; ++ii)
    ext_x[ii] = extractor_x[ssize*e+ii];

  VEC_T::shrink2fit(ext_x);
}

void AExtractor_3D_NURBS_xyz::get_EXT_x(const int &e, double * &ext_x) const
{
  int ssize = (sdegree+1) * (sdegree+1);
  ext_x = new double [ssize];
  for(int ii=0; ii<ssize; ++ii)
    ext_x[ii] = extractor_x[ssize*e+ii];
}

void AExtractor_3D_NURBS_xyz::get_EXT_y(const int &e, std::vector<double> &ext_y) const
{
  ext_y.clear();
  int tsize = (tdegree+1)*(tdegree+1);
  ext_y.resize(tsize);

  for(int ii=0; ii<tsize; ++ii)
    ext_y[ii] = extractor_y[tsize*e + ii];

  VEC_T::shrink2fit(ext_y);
}

void AExtractor_3D_NURBS_xyz::get_EXT_y(const int &e, double * &ext_y) const
{
  int tsize = (tdegree+1) * (tdegree+1);
  ext_y = new double [tsize];
  for(int ii=0; ii<tsize; ++ii)
    ext_y[ii] = extractor_y[tsize*e+ii];
}

void AExtractor_3D_NURBS_xyz::get_EXT_z(const int &e, std::vector<double> &ext_z) const
{
  ext_z.clear();
  int usize = (udegree+1) * (udegree+1);
  ext_z.resize(usize);

  for(int ii=0; ii<usize; ++ii)
    ext_z[ii] = extractor_z[usize*e+ii];

  VEC_T::shrink2fit(ext_z);
}

void AExtractor_3D_NURBS_xyz::get_EXT_z(const int &e, double * &ext_z) const
{
  int usize = (udegree+1) * (udegree+1);
  ext_z = new double [usize];
  for(int ii=0; ii<usize; ++ii)
    ext_z[ii] = extractor_z[usize*e+ii];
}
