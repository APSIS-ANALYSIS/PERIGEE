#ifndef GLOBAL_PART_RELOAD_HPP
#define GLOBAL_PART_RELOAD_HPP
// ============================================================================
// Global_Part_Reload.hpp
// Load the mesh partitioning from HDF5 files.
//
// Author: Ju Liu
// Date: Dec. 9 2021
// ============================================================================
#include "IGlobal_Part.hpp"
#include "HDF5_Reader.hpp"

class Global_Part_Reload : public IGlobal_Part
{
  public:
    Global_Part_Reload(  const int &cpu_size, const int &in_ncommon, 
        const bool &isDualGraph, const std::string &element_part_name = "epart",
        const std::string &node_part_name = "npart" );

    virtual ~Global_Part_Reload();

    virtual idx_t get_epart( const int &ee ) const {return static_cast<idx_t>(epart[ee]);}

    virtual idx_t get_npart( const int &nn, const int &field ) const 
    {return static_cast<idx_t>(npart[nn + field_offset[field]]);}

    virtual bool get_isMETIS() const {return isMETIS;};

    virtual bool get_isDual() const {return isDual;};

    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}

    virtual bool is_serial() const {return isSerial;}

  private:
    bool isMETIS, isDual, isSerial;
    int dual_edge_ncommon;

    std::vector<int> field_offset;
    std::vector<int> epart, npart;
};

#endif
