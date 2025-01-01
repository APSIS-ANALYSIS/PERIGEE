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
#include "QuadPts_vis_tet4.hpp"
#include "QuadPts_vis_tet10.hpp"
#include "QuadPts_vis_hex8.hpp"
#include "QuadPts_vis_hex27.hpp"

class QuadPtsFactory
{
  public:
    static std::unique_ptr<IQuadPts> createVolumeQuadrature(const FEType &elemType, 
        const int &nqp_vol)
    {
      if (nqp_vol <= 0)
      {
        SYS_T::print_fatal("Error: Invalid input. nqp_vol must be positive.\n");
        return nullptr;
      }

      if (elemType == FEType::Tet4 || elemType == FEType::Tet10)
        return SYS_T::make_unique<QuadPts_Gauss_Tet>(nqp_vol);
      else if(elemType == FEType::Hex8 || elemType == FEType::Hex27)
      {
        // For hexahedral elements, check if the volume quadrature points is a perfect 
        // cube
        bool flg = false; 
        const int nqp_vol_1D = isPerfectCube(nqp_vol, flg)
        if(flg == false)
        {
          SYS_T::print_fatal("Error: Invalid nqp_vol for Hex element (cube check).\n");
          return nullptr;
        }
        return SYS_T::make_unique<QuadPts_Gauss_Hex>(nqp_vol_1D);
      }
      else
      {
        SYS_T::print_fatal("Error: FEType not supported.\n");
        return nullptr;
      }
    }

    static std::unique_ptr<IQuadPts> createSurfaceQuadrature(const FEType &elemType, 
        const int &nqp_sur)
    {
      if (nqp_sur <= 0)
      {
        SYS_T::print_fatal("Error: Invalid input. nqp_sur must be positive.\n");
        return nullptr;
      }

      if (elemType == FEType::Tet4 || elemType == FEType::Tet10)
        return SYS_T::make_unique<QuadPts_Gauss_Triangle>(nqp_sur);
      else if(elemType == FEType::Hex8 || elemType == FEType::Hex27)
      {
        // For hexahedral elements, check if the surface quadrature points is a perfect 
        // square
        bool flg = false; 
        const int nqp_vol_1D = isPerfectSquare(nqp_sur, flg)
        if(flg == false)
        {
          SYS_T::print_fatal("Error: Invalid nqp_sur for Hex element (square check).\n");
          return nullptr;
        }
        return SYS_T::make_unique<QuadPts_Gauss_Quad>(nqp_vol_1D);
      }
      else
      {
        SYS_T::print_fatal("Error: FEType not supported.\n");
        return nullptr;
      }
    }
    
    static std::unique_ptr<IQuadPts> createVisQuadrature(const FEType &elemType)
    {
      switch (elemType)
      {
        case FEType::Tet4:
          return std::make_unique<QuadPts_vis_Tet4>();
        case FEType::Tet10:
          return std::make_unique<QuadPts_vis_Tet10>();
        case FEType::Hex8:
          return std::make_unique<QuadPts_vis_Hex8>();
        case FEType::Hex27:
          return std::make_unique<QuadPts_vis_Hex27>();
        default:
          SYS_T::print_fatal("Error: FEType not supported.\n");
          return nullptr;
      }
    }

  private:
    // Check if a number is a perfect cube
    static int isPerfectCube(const int &nn, bool &flg)
    {
      const int root = static_cast<int>(std::round(std::cbrt(static_cast<double>(nn))));
      flg = (root * root * root == nn);
      return root;
    }

    // Check if a number is a perfect square
    static int isPerfectSquare(const int &nn, bool &flg)
    {
      const int root = static_cast<int>(std::round(std::sqrt(static_cast<double>(nn))));
      flg = (root * root == nn);
      return root;
    }
};

#endif
