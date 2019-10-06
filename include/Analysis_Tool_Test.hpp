#ifndef ANALYSIS_TOOL_TEST_HPP
#define ANALYSIS_TOOL_TEST_HPP
// ==================================================================
// Analysis_Tool_Test.hpp
// This is a collection of test functions for testing the data structures
// in Analysis Tools.
//
// Date: Nov. 8th 2013
// ==================================================================
#include <iostream>
#include <cassert>
#include <cstdlib>

#include "Sys_Tools.hpp"
#include "HDF5_PartReader.hpp"
#include "IAExtractor.hpp"
#include "AExtractor_3D_NURBS_xyz.hpp"
#include "ALocal_IEN.hpp"
#include "FEANode.hpp"
#include "APart_Node.hpp"

// Make sure the Extraction operator are read correctly
void AnalysisTool_EXT_check(const std::string &part_file, const int &rank);

// This routine checkes the correctness of ALocal_IEN class
void AnalysisTool_LIEN_check(const HDF5_PartReader * const &h5reader);

#endif
