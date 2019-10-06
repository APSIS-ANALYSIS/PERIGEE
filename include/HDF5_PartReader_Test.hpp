#ifndef HDF5_PART_READER_TEST_HPP
#define HDF5_PART_READER_TEST_HPP
// ==================================================================
// HDF5_PartReader_Test.hpp
// Description:
// Perform tests / print job for HDF5_PartReader class.
//
// Date:
// Oct. 20 2013
// ==================================================================
#include <iostream>
#include <cstdlib>

#include "HDF5_PartReader.hpp"
#include "Vec_Tools.hpp"
#include "IPart.hpp"
#include "BC_Partition.hpp"

using std::endl;
using std::cout;
using std::cerr;

void HDF5_PartReader_Print_GMI( const class HDF5_PartReader * const &preader );
void HDF5_PartReader_Print_EXT( const class HDF5_PartReader * const &preader,
   const int &e );
void HDF5_PartReader_Print_LIEN( const class HDF5_PartReader * const &preader,
    const int &e );
void HDF5_PartReader_Print_LE( const class HDF5_PartReader * const &preader );
void HDF5_PartReader_Print_LN( const class HDF5_PartReader * const &preader );
void HDF5_PartReader_Print_PI( const class HDF5_PartReader * const &preader );

void HDF5_PartReader_Print_CPL( const class HDF5_PartReader * const &preader );

void HDF5_PartReader_Print_BC_LID_dof( const class HDF5_PartReader * const &preader );
void HDF5_PartReader_Print_BC_LDN( const class HDF5_PartReader * const &preader );
void HDF5_PartReader_Print_BC_LP( const class HDF5_PartReader * const &preader );
void HDF5_PartReader_Print_BC_BCE( const class HDF5_PartReader * const &preader );

// CHECK READ USING IPART
void HDF5_PartReader_Check( const class HDF5_PartReader * const &preader,
    const class IPart * const &part );

void HDF5_PartReader_Check_CtrlPts( const class HDF5_PartReader * const &preader,
    const class IPart * const &part );

void HDF5_PartReader_Check_LIEN( const class HDF5_PartReader * const &preader,
    const class IPart * const &part );

void HDF5_PartReader_Check_LID( const class HDF5_PartReader * const &preader,
    const class BC_Partition * const &bcpart );
#endif
