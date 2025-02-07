#ifndef FEANODE_HPP
#define FEANODE_HPP
// ============================================================================
// FEANode.hpp
// Interface for Finite Element Node
// Object:
// This class should 
// 1. record the [x,y,z,w] corrdinates of all control points for 
//    the local partition, if the geometry is described by NURBS/T-splines.
//
// 2. record the [x,y,z] coordinates of all control points for 
//    the local partition, if the geometry is described by classical FEM mesh. 
//
// Constructed from inputs of the file basename and cpu_rank, and read using
// the  HDF5_Reader tool.
//
// Author: Ju Liu
// Date Created:  Nov.  5 2013
// Date Modified: Jan. 20 2017 
// ============================================================================
#include "HDF5_Reader.hpp"
#include "Vector_3.hpp"

class FEANode
{
  public:
    // ------------------------------------------------------------------------
    // ! Constructor
    //   Read ctrl_x/y/z directly from the part file with the given base name 
    //   and cpu index. If weights exist in the h5 file, they will also be read; 
    //   otherwise, the ctrlPts_w vector will be empty.
    // ------------------------------------------------------------------------
    FEANode( const std::string &fileBaseName, const int &cpu_rank );

    // ------------------------------------------------------------------------
    // ! Destructor
    // ------------------------------------------------------------------------
    virtual ~FEANode() = default;

    // ------------------------------------------------------------------------
    // ! Functions that give access to the coordinates (and weights).
    //   Input: index ranges in [ 0 , nlocghonode )
    // ------------------------------------------------------------------------
    virtual Vector_3 get_ctrlPts_xyz(const int &index) const
    { return Vector_3( ctrlPts_x[index], ctrlPts_y[index], ctrlPts_z[index] ); }
    
    virtual double get_ctrlPts_x(const int &index) const {return ctrlPts_x[index];}
    
    virtual double get_ctrlPts_y(const int &index) const {return ctrlPts_y[index];}
    
    virtual double get_ctrlPts_z(const int &index) const {return ctrlPts_z[index];}
    
    virtual double get_ctrlPts_w(const int &index) const {return ctrlPts_w[index];}

    // ------------------------------------------------------------------------
    // ! Get n control points' x-y-z in a batch
    //   Input : \para num  : the number of control points (n) to be passed out
    //           \para index: index[0] - index[num-1] ctrlPts_* will be passed 
    //                        to ctrl_x, ctrl_y, ctrl_z
    //   Usually, num = nLocBas, and index = LIEN[e][]
    //   Note: Users are responsible for allocating and deallocating memory
    //         for index and ctrl_(x/y/z).
    // ------------------------------------------------------------------------
    virtual void get_ctrlPts_xyz( const int &num, const int * const &index, 
        double * const &ctrl_x, double * const &ctrl_y, 
        double * const &ctrl_z ) const;

    virtual std::array<std::vector<double>, 3> get_ctrlPts_xyz( 
        const std::vector<int> &index ) const;

    // ------------------------------------------------------------------------
    // ! Get n control points' x-y-z-w in a batch
    //   Note: Users are responsible for allocating and deallocating memory
    //         for index and ctrl_(x/y/z/w).
    // ------------------------------------------------------------------------
    virtual void get_ctrlPts_xyzw( const int &num, const int * const &index, 
        double * const &ctrl_x, double * const &ctrl_y, 
        double * const &ctrl_z, double * const &ctrl_w ) const;
   
    // ------------------------------------------------------------------------
    // ! Get n control points' x-y-w in a batch
    //   NOTE: Users are responsible for allocating and deallocating memory
    //         for index and ctrl_(x/y/w).
    // ------------------------------------------------------------------------
    virtual void get_ctrlPts_xyw( const int &num, const int * const &index,
        double * const &ctrl_x, double * const &ctrl_y,
        double * const &ctrl_w ) const;

    // ------------------------------------------------------------------------
    // ! Get n contrl points' x-y in a batch
    //   NOTE: Users are responsible for allocating and deallocating memory
    //         for index and ctrl_(x/y).
    // ------------------------------------------------------------------------
    virtual void get_ctrlPts_xy( const int &num, const int * const &index,
        double * const &ctrl_x, double * const &ctrl_y ) const;

    // ------------------------------------------------------------------------
    // ! Print the info for this class.
    // ------------------------------------------------------------------------
    virtual void print_info() const;

    // ------------------------------------------------------------------------
    // ! Returns the memory usage of the private data in Bytes
    // ------------------------------------------------------------------------
    virtual double get_memory_usage() const;

  private:
    // The control points' coordinates and weights if used for NURBS
    std::vector<double> ctrlPts_x {}, ctrlPts_y {}, ctrlPts_z {}, ctrlPts_w {};

    FEANode() = delete;
};

#endif
