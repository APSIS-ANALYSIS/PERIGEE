#ifndef WBC_TOOLS_HPP
#define WBC_TOOLS_HPP
// ============================================================================
// WBC_Tools.hpp
// ----------------------------------------------------------------------------
// This is a temporary file for the implementation of weakly enforced Dirichlet
// boundary condition.
// WBC_T namespace contains the temporary functions which are necessary for the 
// formulation or the algorithm. Perhaps many of them will be sent to specific
// classes as member functions, and this file will be removed in future.
// ============================================================================

#include "Sys_Tools.hpp"
#include "Vector_3.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "FEAElement.hpp"

namespace WBC_T
{
  
}

// ----------------------------------------------------------------
// class QuadPts_on_face : express a triangular quadrature rule on a tet's face
//                         or express a quadrilateral quadrature rule on a hex's face
//                         with volume coordinate
//  
// Input: \para qp_surface : the quadrature rule of a surface element
//        \para face_id    : on which face of the volume element
// ----------------------------------------------------------------
class QuadPts_on_face : public IQuadPts
{ 
  public:
    // The input IQuadPts should be QuadPts_Gauss_Quad or QuadPts_Gauss_Triangle
    QuadPts_on_face(const IQuadPts const * qp_surface, const int &face_id)
    { 
      num_pts = qp_surface->get_num_quadPts();
      dim = qp_surface->get_dim() + 1;

      if(dim == 4) // tet
      {
        qp.assign(4 * num_pts, 0.0);
        switch (face_id)
        {
          case 0: // u = 0 : node3 = node0', node1 = node1', node2 = node2'
          {
            for(unsigned int ii{0}; ii < num_pts; ++ii)
            {
              qp[4*ii + 0] = qp_surface->get_qp(3*ii + 0);  // r = r'
              qp[4*ii + 1] = qp_surface->get_qp(3*ii + 1);  // s = s'
              qp[4*ii + 2] = qp_surface->get_qp(3*ii + 2);  // t = t'
            }
          } break;
          case 1: // r = 0 : node0 = node0', node3 = node1', node2 = node2'
          {
            for(unsigned int ii{0}; ii < num_pts; ++ii)
            {
              qp[4*ii + 1] = qp_surface->get_qp(3*ii + 1);  // s = s'
              qp[4*ii + 2] = qp_surface->get_qp(3*ii + 0);  // t = r'
              qp[4*ii + 3] = qp_surface->get_qp(3*ii + 2);  // u = t'
            }
          } break;
          case 2: // s = 0 : node0 = node0', node1 = node1', node3 = node2'
          {
            for(unsigned int ii{0}; ii < num_pts; ++ii)
            {
              qp[4*ii + 0] = qp_surface->get_qp(3*ii + 0);  // r = r'
              qp[4*ii + 2] = qp_surface->get_qp(3*ii + 1);  // t = s'
              qp[4*ii + 3] = qp_surface->get_qp(3*ii + 2);  // u = t'
            }
          } break;
          case 3: // t = 0 : node0 = node0', node2 = node1', node1 = node2'
          {
            for(unsigned int ii{0}; ii < num_pts; ++ii)
            {
              qp[4*ii + 0] = qp_surface->get_qp(3*ii + 1);  // r = s'
              qp[4*ii + 1] = qp_surface->get_qp(3*ii + 0);  // s = r'
              qp[4*ii + 3] = qp_surface->get_qp(3*ii + 2);  // u = t'
            }
          } break;
          default:
            SYS_T::print_fatal("Error: QuadPts_on_face, wrong face id for tet.");
            break;
        }
      }
      else if(dim == 3) // hex
      {
        ; // Unimplemented
      }
      else
        SYS_T::print_fatal("Error: QuadPts_on_face, unsupported IQuadPts.\n");

      qw.resize(num_pts);
      for(unsigned int ii{0}; ii < num_pts; ++ii)
        qw[ii] = qp_surface->get_qw(ii);
    }

    ~QuadPts_on_face(){}

    virtual int get_dim() const {return dim;};

    virtual int get_num_quadPts() const {return num_pts;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[dim*ii + comp];}

    virtual double get_qw(unsigned int ii) const
    {return qw[ii];}

    virtual void print_info() const { ;}

  private:
    int num_pts;
    std::vector<double> qp {};
    std::vector<double> qw {};

    int dim;
};

#endif