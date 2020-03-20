#ifndef FEANODE_HPP
#define FEANODE_HPP
// ==================================================================
// FEANode.hpp
// Interface for Finite Element Node
// Object:
// This class should 
// 1. record the [x,y,z,w] corrdinates of all the control point for 
//    the subdomain, if the geometry is described by NURBS/T-splines,
//    construct by using the HDF5_PartReader tool.
//
// 2. record the [x,y,z] coordinates of all the control points for 
//    the domain, if the geometry is described by classical FEM mesh. 
//    Construct this type of data structure by input the filebasename,
//    and cpu_rank, and read by using the HDF5_Reader tool.
//
// Date Created: Nov. 5th 2013
// Date Modified: Jan. 20 2017 
// ==================================================================
#include "HDF5_Reader.hpp"

class FEANode
{
  public:
    // -------------------------------------------------------------- 
    // Read ctrl_x/y/z directly from part file with base name and 
    // given cpu index. If the weight exists in the h5 file, the 
    // constructor will read its value as well; otherwise, the ctrlPts_w
    // vector will be empty.
    // -------------------------------------------------------------- 
    FEANode( const std::string &fileBaseName, const int &cpu_rank );

    virtual ~FEANode();

    virtual double get_ctrlPts_x(const int &index) const 
    {return ctrlPts_x[index];}
    
    virtual double get_ctrlPts_y(const int &index) const 
    {return ctrlPts_y[index];}
    
    virtual double get_ctrlPts_z(const int &index) const 
    {return ctrlPts_z[index];}
    
    virtual double get_ctrlPts_w(const int &index) const 
    {return ctrlPts_w[index];}

    // -------------------------------------------------------------- 
    // Get n control points' x-y-z index in a batch
    // \para num: the number of control points to be passed out
    // \para index: index[0] - index[num-1] ctrlPts_* will be passed 
    //              out to ctrl_x, ctrl_y, ctrl_z
    // usually, num = nLocBas, and index = LIEN[e][]
    // Note: Users are responsible for allocating num solts for both
    //       index and ctrl_..; also delete these pointers after usage.
    // --------------------------------------------------------------
    virtual void get_ctrlPts_xyz( const int &num, const int * const &index, 
        double * const &ctrl_x, double * const &ctrl_y, 
        double * const &ctrl_z ) const;
    
    // --------------------------------------------------------------
    // Get n control points' x-y-z-w index in a batch
    // \para num: the number of control points to be passed out
    // \para index: index[0] - index[num-1] ctrlPts_* will be passed 
    //              out to ctrl_x, ctrl_y, ctrl_z, ctrl_w
    // usually, num = nLocBas, and index = LIEN[e][]
    // Note: Users are responsible for allocating num solts for both 
    //       index and ctrl_..; also delete these pointers after usage.
    // --------------------------------------------------------------
    virtual void get_ctrlPts_xyzw( const int &num, const int * const &index, 
        double * const &ctrl_x, double * const &ctrl_y, 
        double * const &ctrl_z, double * const &ctrl_w ) const;
   
    // --------------------------------------------------------------
    // get n contrl points' x-y-w in a batch
    // Users are responsible for allocating num solts for index, 
    // ctrl_... and free the memory for these pointers after usage.
    // --------------------------------------------------------------
    virtual void get_ctrlPts_xyw( const int &num, const int * const &index,
        double * const &ctrl_x, double * const &ctrl_y,
        double * const &ctrl_w ) const;

    // --------------------------------------------------------------
    // get n contrl points' x-y in a batch
    // Users are responsible for allocating num solts for index, 
    // ctrl_... and free the memory for these pointers after usage.
    // --------------------------------------------------------------
    virtual void get_ctrlPts_xy( const int &num, const int * const &index,
        double * const &ctrl_x, double * const &ctrl_y ) const;

    virtual void print_info() const;

    // --------------------------------------------------------------
    // get_memory_usage()
    // returns the memory usage of the private data in Byte
    // --------------------------------------------------------------
    virtual double get_memory_usage() const;

  private:
    std::vector<double> ctrlPts_x, ctrlPts_y, ctrlPts_z, ctrlPts_w;
};

#endif
