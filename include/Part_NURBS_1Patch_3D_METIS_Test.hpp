#ifndef PART_NURBS_1PATCH_3D_METIS_TEST_HPP
#define PART_NURBS_1PATCH_3D_METIS_TEST_HPP
// ==================================================================
// Part_NURBS_1Patch_3D_METIS_Test.hpp
// Object:
// Perform tests to verify the functions in Part_NURBS_1Patch_3D_METIS
// class.
//
// Date:
// Oct. 4th 2013
// ==================================================================
#include <cassert>

#include "IGlobal_Part.hpp"
#include "Part_NURBS_1Patch_3D_METIS.hpp"
#include "IIEN.hpp"
#include "IMesh.hpp"
#include "Map_Node_Index.hpp"

using std::cout;
using std::endl;
using std::cerr;

void Part_NURBS_1Patch_3D_METIS_LIEN_Test(
    const class IPart * const &part,
    const class Map_Node_Index * const &mnindex,
    const class IIEN * const &IEN );

void Part_NURBS_1Patch_3D_METIS_Node_Test(
    const class IPart * const &part,
    const class Map_Node_Index * const &mnindex);

void Part_NURBS_1Patch_3D_METIS_Mesh_Test(
    const class IPart * const &part,
    const class IMesh * const &mesh);

void Part_NURBS_1Patch_3D_METIS_CtrlPts_Test(
    const class IPart * const &part,
    const class Map_Node_Index * const &mnindex,
    const std::vector<double> &ctrlPts );

void Part_NURBS_1Patch_3D_METIS_Print_Node(
    const class IMesh * const &mesh,
    const class Map_Node_Index * const &mnindex,
    const class IGlobal_Part * const &gpart );

void Part_NURBS_1Patch_3D_METIS_Print_Elem(
    const class IMesh * const &mesh,
    const class IGlobal_Part * const &gpart );

#endif
