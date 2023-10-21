#ifndef ELEMBC_3D_WEAK_HPP
#define ELEMBC_3D_WEAK_HPP
// ============================================================================
// ElemBC_3D_weak.hpp
//
// This is an instantation of ElemBC_3D for weakly enforced Dirichlet boundary 
// condition on the wall of a fluid domain or on a moving interface of two
// fluid domains.
//
// Since the surface integral involves the first order derivative of basis function,
// we will use the local ien of volume element that has faces on the surface.
// Hence we should write the following data of each boundary into h5 file:
//
// 1. Number of surface elements;
// 2. Face2elem id, i.e. the global id of the volume element where the surface
//    element is attached. We will use it to find local ien of volume element
//    in local assembly.
// 3. Face id of the volume element, we will use the face id and the local ien of
//    volume element to build the quadrature rule and the basis;
// For type 2 weak BC:
// 4. Number of surface nodes;
// 5. Global id of surface nodes;
// 6. Rotation matrices at each surface node.
//
// To be distinguished from other ElemBC_3D, it uses EBC_Partition_weak to 
// write h5 file and serves a unique ALocal_WeakBC object.
//
// Author: Xuanming Huang
// Date Created: Oct. 20th  2023
// ============================================================================

#include "ElemBC_3D.hpp"

class ElemBC_3D_weak : public ElemBC_3D {
  public:
    ElemBC_3D_weak( const std::vector<std::string> &vtkfileList,
                    const int &in_weak_bc_type,
                    const double &in_C_bI,
                    const IIEN * const &VIEN,
                    const int &elemtype );

    virtual ~ElemBC_3D_weak();

  private:
    // Weak bc type
    // type = 0 : strongly no-slip bc in all direction (do nothing)
    // type = 1 : weakly no-slip bc in all direction
    // type = 2 : strongly no-slip bc in wall-normal direction,
    //            and weakly no-slip bc in wall-tangent direction
    const int weak_bc_type;

    // value of coefficient C_bI
    const double C_bI;

    // the face id of the volume element
    std::vector< std::vector<int> > face_id;

    // ------------------------------------------------------------------------
    // Disallow default constructor
    ElemBC_3D_weak() = delete;
     
};

#endif