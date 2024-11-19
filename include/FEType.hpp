#ifndef FETYPE_HPP
#define FETYPE_HPP
// ============================================================================
// FEType.hpp
// Object: This is an enumerate class for Finite Element classes.
//
// Author: Chi Ding
// Date:   Nov. 19, 2024
// ============================================================================

enum class FEType
{
  Unknown,
  Tet4,
  Tet10,
  Tet10_v2,
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

#endif
