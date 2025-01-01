#ifndef FEAELEMENTFACTORY_HPP
#define FEAELEMENTFACTORY_HPP

#include "FEAElement.hpp"

class ElementFactory
{
  public:
    static std::unique_ptr<FEAElement> createVolumeElement(FEType elemType, const int &nqp)
    {
      switch(elemType)
      {
        case FEType::Tet4:
          return SYS_T::make_unique<FEAElement_Tet4>(nqp);
        case FEType::Tet10:
          return SYS_T::make_unique<FEAElement_Tet10>(nqp);
        case FEType::Hex8:
          return SYS_T::make_unique<FEAElement_Hex8>(nqp);
        case FEType::Hex27:
          return SYS_T::make_unique<FEAElement_Hex27>(nqp);
        default:
          SYS_T::print_fatal("Error: Volume element type not supported.\n");
          return nullptr;
      }
    }

    static std::unique_ptr<FEAElement> createSurfaceElement(FEType elemType, const int &nqp)
    {
      switch (elemType)
      {
        case FEType::Tet4:
          return SYS_T::make_unique<FEAElement_Triangle3_3D_der0>(nqp);
        case FEType::Tet10:
          return SYS_T::make_unique<FEAElement_Triangle6_3D_der0>(nqp);
        case FEType::Hex8:
          return SYS_T::make_unique<FEAElement_Quad4_3D_der0>(nqp);
        case FEType::Hex27:
          return SYS_T::make_unique<FEAElement_Quad9_3D_der0>(nqp);
        default:
          SYS_T::print_fatal("Error: Surface element type not supported.\n");
          return nullptr;
      }
    }
};

#endif
