#ifndef QUADPTSFACTORY_HPP
#define QUADPTSFACTORY_HPP
// ============================================================================
// FEAElementFactory.hpp
//
// This is the factory that construct an element class based on the FEType and
// number of quadrature points.
//
// Author: Ju Liu, liujuy@gmail.com
// Date: Jan. 1st 2025
// ============================================================================
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"

class QuadPtsFactory
{
  public:
    static std::unique_ptr<IQuadPts> createVolumeQuadrature(const FEType &elemType, 
        const int &nqp_vol)
    {
      switch (elemType)
      {
        case FEType::Tet4:
          return SYS_T::make_unique<QuadPts_Gauss_Tet>(nqp_vol);
        case FEType::Tet10:
          return SYS_T::make_unique<QuadPts_Gauss_Tet>(nqp_vol);
        case FEType::Hex8:
          return SYS_T::make_unique<QuadPts_Gauss_Hex>(nqp_vol_1D);
        case FEType::Hex27:
          return SYS_T::make_unique<QuadPts_Gauss_Hex>(nqp_vol_1D);
        default:
          SYS_T::print_fatal("Error: Volume quadrature type not supported for this element.\n");
          return nullptr;
      }
    }

    static std::unique_ptr<IQuadPts> createSurfaceQuadrature(const FEType &elemType, 
        const int &nqp_sur)
    {
      switch (elemType)
      {
        case FEType::Tet4:
          return SYS_T::make_unique<QuadPts_Gauss_Triangle>(nqp_sur);
        case FEType::Tet10:
          return SYS_T::make_unique<QuadPts_Gauss_Triangle>(nqp_sur);
        case FEType::Hex8:
          return SYS_T::make_unique<QuadPts_Gauss_Quad>(nqp_sur_1D);
        case FEType::Hex27:
          return SYS_T::make_unique<QuadPts_Gauss_Quad>(nqp_sur_1D);
        default:
          SYS_T::print_fatal("Error: Surface quadrature type not supported for this element.\n");
          return nullptr;
      }
    }

  private:
    // Check if a number is a perfect cube
    static bool isPerfectCube(const int &nn, int &root)
    {
      root = static_cast<int>(std::round(std::cbrt(static_cast<double>(nn))));
      return (root * root * root == nn);
    }

    // Check if a number is a perfect square
    static bool isPerfectSquare(const int &nn, int &root)
    {
      root = static_cast<int>(std::round(std::sqrt(static_cast<double>(nn))));
      return (root * root == nn);
    }
};

#endif
