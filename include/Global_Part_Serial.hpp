#ifndef GLOBAL_PART_SERIAL_HPP
#define GLOBAL_PART_SERIAL_HPP
// ==================================================================
// Global_Part_Serial.hpp
// Object:
// Trival mesh file for serial run
// element partition information: epart
// node partition information: npart.
//
// Date: Oct 2nd 2013
// ==================================================================
#include "Sys_Tools.hpp"
#include "IGlobal_Part.hpp"
#include "HDF5_Writer.hpp"

class Global_Part_Serial : public IGlobal_Part
{
  public:
    Global_Part_Serial( const int &in_nelem, const int &in_nfunc,
        const std::string &element_part_name = "epart",
        const std::string &node_part_name = "npart" );

    Global_Part_Serial( const int &num_fields,
        const std::vector<int> &in_nelem_list,
        const std::vector<int> &in_nfunc_list,
        const std::string &element_part_name = "epart",
        const std::string &node_part_name = "npart" );

    ~Global_Part_Serial() override;

    idx_t get_epart( int ee ) const override {return epart[ee];}

    idx_t get_npart( int nn, int field = 0 ) const override
    {return npart[nn + field_offset[field]];}

    bool get_isMETIS() const override {return false;};

    bool get_isDual() const override {return false;};

    int get_dual_edge_ncommon() const override {return 0;}

    bool is_serial() const override {return true;}

  private:
    idx_t * epart, * npart;

    std::vector<int> field_offset;

    virtual void write_part_hdf5( const std::string &fileName, 
        const idx_t * const &part_in, const int &part_size ) const;
};

#endif
