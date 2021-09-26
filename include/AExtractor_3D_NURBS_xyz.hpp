#ifndef AEXTRACTOR_3D_NURBS_XYZ_HPP
#define AEXTRACTOR_3D_NURBS_XYZ_HPP
// ==================================================================
// AExtractor_3D_NURBS_xyz.hpp
// This is a collector of local elements' extractors. It is used for
// 3D NURBS, and the extractor are given in each direction.
//
// Date:
// Nov. 8 2013
// ==================================================================
#include <string>
#include <iostream>
#include <vector>
#include "IAExtractor.hpp"

class AExtractor_3D_NURBS_xyz : public IAExtractor
{
  public:
    AExtractor_3D_NURBS_xyz( const std::string &part_file, const int &rank );
    
    AExtractor_3D_NURBS_xyz( const HDF5_PartReader * const &h5reader );
    
    virtual ~AExtractor_3D_NURBS_xyz();

    virtual void get_EXT_x(const int &e, std::vector<double> &ext_x) const;
    virtual void get_EXT_y(const int &e, std::vector<double> &ext_y) const;
    virtual void get_EXT_z(const int &e, std::vector<double> &ext_z) const;

    virtual void get_EXT_x(const int &e, double * &ext_x) const;
    virtual void get_EXT_y(const int &e, double * &ext_y) const;
    virtual void get_EXT_z(const int &e, double * &ext_z) const;
  private:
    int sdegree, tdegree, udegree;
    std::vector<double> extractor_x, extractor_y, extractor_z; 
};
#endif
