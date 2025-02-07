#ifndef FETYPE_HPP
#define FETYPE_HPP
// ============================================================================
// FEType.hpp
// Object: This is an enumerate class for Finite Element classes.
//
// Author: Chi Ding
// Date:   Nov. 19, 2024
// ============================================================================
#include <string>
#include <unordered_map>

enum class FEType
{
  Unknown,
  Tet4,
  Tet10,
  Hex8,
  Hex27,
  Tri3,
  Tri6,
  Quad4,
  Quad9,
  Tri3_der0,
  Tri6_der0,
  Quad4_der0,
  Quad9_der0
};

namespace FE_T 
{
  // Convert FEType to std::string
  inline std::string to_string(const FEType &type) 
  {
    switch (type) 
    {
      case FEType::Unknown: return "Unknown";
      case FEType::Tet4: return "Tet4";
      case FEType::Tet10: return "Tet10";
      case FEType::Hex8: return "Hex8";
      case FEType::Hex27: return "Hex27";
      case FEType::Tri3: return "Tri3";
      case FEType::Tri6: return "Tri6";
      case FEType::Quad4: return "Quad4";
      case FEType::Quad9: return "Quad9";
      case FEType::Tri3_der0: return "Tri3_der0";
      case FEType::Tri6_der0: return "Tri6_der0";
      case FEType::Quad4_der0: return "Quad4_der0";
      case FEType::Quad9_der0: return "Quad9_der0";
      default: return "Unknown";
    }
  }

  // Convert std::string to FEType
  inline FEType to_FEType(const std::string &str) 
  {
    static const std::unordered_map<std::string, FEType> string_to_fetype = 
    {
      {"Unknown", FEType::Unknown},
      {"Tet4", FEType::Tet4},
      {"Tet10", FEType::Tet10},
      {"Hex8", FEType::Hex8},
      {"Hex27", FEType::Hex27},
      {"Tri3", FEType::Tri3},
      {"Tri6", FEType::Tri6},
      {"Quad4", FEType::Quad4},
      {"Quad9", FEType::Quad9},
      {"Tri3_der0", FEType::Tri3_der0},
      {"Tri6_der0", FEType::Tri6_der0},
      {"Quad4_der0", FEType::Quad4_der0},
      {"Quad9_der0", FEType::Quad9_der0}
    };

    auto it = string_to_fetype.find(str);
    if (it != string_to_fetype.end())
      return it->second;
    else
      return FEType::Unknown;
  }
  
  // Convert FEType to its number of local basis functions
  inline int to_nLocBas(const FEType &type) 
  {
    switch (type) 
    {
      case FEType::Unknown: return 0;
      case FEType::Tet4: return 4;
      case FEType::Tet10: return 10;
      case FEType::Hex8: return 8;
      case FEType::Hex27: return 27;
      case FEType::Tri3: return 3;
      case FEType::Tri6: return 6;
      case FEType::Quad4: return 4;
      case FEType::Quad9: return 9;
      case FEType::Tri3_der0: return 3;
      case FEType::Tri6_der0: return 6;
      case FEType::Quad4_der0: return 4;
      case FEType::Quad9_der0: return 9;
      default: return 0;
    }
  }
}

#endif
