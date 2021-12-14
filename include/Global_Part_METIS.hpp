#ifndef GLOBAL_PART_METIS_HPP
#define GLOBAL_PART_METIS_HPP
// ==================================================================
// Global_Part_METIS.hpp
// Object:
// Use METIS routien to get the global mesh partition stored as:
// element partition information: epart
// node partition information: npart.
//
// Date Created: Oct 2nd 2013
// ==================================================================
#include "IMesh.hpp"
#include "IIEN.hpp"
#include "IGlobal_Part.hpp"
#include "HDF5_Writer.hpp"

class Global_Part_METIS : public IGlobal_Part
{
  public:
    // Constructor:
    // It will create eptr and eind arrays and call METIS_PartMeshDual or
    // METIS_PartMeshNodal for mesh partition
    Global_Part_METIS( const int &cpu_size,
        const int &in_ncommon, const bool &isDualGraph,
        const IMesh * const &mesh,
        const IIEN * const &IEN,
        const std::string &element_part_name = "epart",
        const std::string &node_part_name = "npart" );

    // It will call METIS to partition a multi-field mesh. An assumption is that
    // the nElem in all meshes are the same.
    Global_Part_METIS( const int &num_fields, const int &cpu_size,
        const int &in_ncommon, const bool &isDualGraph,
        const std::vector<IMesh const *> &mesh_list,
        const std::vector<IIEN const *> &IEN_list,
        const std::string &element_part_name = "epart",
        const std::string &node_part_name = "npart" );

    virtual ~Global_Part_METIS();

    virtual idx_t get_epart( const int &ee ) const {return epart[ee];}

    virtual idx_t get_npart( const int &nn ) const {return npart[nn];}

    virtual bool get_isMETIS() const {return isMETIS;};

    virtual bool get_isDual() const {return isDual;};

    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}

    virtual bool is_serial() const {return false;}

  private:
    const bool isMETIS, isDual;
    const int dual_edge_ncommon;

    idx_t * epart;
    idx_t * npart;

    virtual void write_part_hdf5( const std::string &fileName, 
        const idx_t * const &part_in,
        const int &part_size, const int &cpu_size,
        const bool &part_isdual, const int &in_ncommon,
        const bool &isMETIS ) const;

    // --------------------------------------------------------------
    // This function will write the data of part_in in 64bit HDF5 format. 
    // This function should be called if idx_t is the int64_t.
    // --------------------------------------------------------------
    virtual void write_part_hdf5_64bit( const std::string &fileName, 
        const int64_t * const &part_in,
        const int64_t &part_size, const int &cpu_size,
        const bool &part_isdual, const int &in_ncommon,
        const bool &isMETIS ) const;
};

#endif
