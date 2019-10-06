#ifndef AEXTRACTOR_2D_NURBS_XY_HPP
#define AEXTRACTOR_2D_NURBS_XY_HPP
// ==================================================================
// AExtractor_2D_NURBS_xy.hpp
// This is a collector of local elements' extractors. It is used for
// 2D NURBS, and the extractor are given in each direction (x and y).
//
// Date: April 13 2014
// ==================================================================
#include "IAExtractor.hpp"
#include "HDF5_PartReader.hpp"

class AExtractor_2D_NURBS_xy : public IAExtractor
{
  public:
    AExtractor_2D_NURBS_xy( const std::string &part_file, const int &rank );
    
    AExtractor_2D_NURBS_xy( const HDF5_PartReader * const &h5reader );
    
    virtual ~AExtractor_2D_NURBS_xy();

    virtual void get_EXT_x(const int &e, std::vector<double> &ext_x) const;
    virtual void get_EXT_y(const int &e, std::vector<double> &ext_y) const;
    virtual void get_EXT_z(const int &e, std::vector<double> &ext_z) const;

    virtual void get_EXT_x(const int &e, double * &ext_x) const;
    virtual void get_EXT_y(const int &e, double * &ext_y) const;
    virtual void get_EXT_z(const int &e, double * &ext_z) const;
  private:
    int sdegree, tdegree;
    std::vector<double> extractor_x, extractor_y; 
};
#endif
