#ifndef GLOBAL_PART_METIS_HPP
#define GLOBAL_PART_METIS_HPP
// ============================================================================
// Global_Part_METIS.hpp
// Object:
// Use METIS routien to get the global mesh partition stored as:
// element partition information: epart
// node partition information: npart.
// 
// Author: Ju Liu
// Date Created: Oct 2nd 2013
// ============================================================================
#include "Sys_Tools.hpp"
#include "IIEN.hpp"
#include "IGlobal_Part.hpp"
#include "HDF5_Writer.hpp"

class Global_Part_METIS : public IGlobal_Part
{
  public:
    // ------------------------------------------------------------------------
    // Constructor:
    // It will create eptr and eind arrays and call METIS_PartMeshDual or
    // METIS_PartMeshNodal for mesh partition
    // ------------------------------------------------------------------------
    Global_Part_METIS( const int &cpu_size,
        const int &in_ncommon, const bool &isDualGraph,
        const int &in_nelem, const int &in_nfunc, const int &in_nlocbas, 
        const IIEN * const &IEN,
        const std::string &element_part_name = "epart",
        const std::string &node_part_name = "npart" );

    // ------------------------------------------------------------------------
    // It will call METIS to partition a multi-field mesh. An assumption is that
    // the nElem in all meshes are the same.
    // An example to think in mind could be the Taylor-Hood discretization for
    // Stokes problems, where we use two field (pressure and velocity). The IEN
    // for pressure is for 4 node tet, while the IEN for velocity is for 10 node
    // tet. Therefore, we create a 14 node connectivity array for the mesh
    // partitioning. This ensures that the mesh for pressure and the mesh for
    // velocity have the same partition results in terms of element partition.
    // ------------------------------------------------------------------------
    Global_Part_METIS( const int &num_fields, const int &cpu_size,
        const int &in_ncommon, const bool &isDualGraph,
        const std::vector<int> &in_nelem_list,
        const std::vector<int> &in_nfunc_list,
        const std::vector<int> &in_nlocbas_list,
        const std::vector<IIEN const *> &IEN_list,
        const std::string &element_part_name = "epart",
        const std::string &node_part_name = "npart" );

    virtual ~Global_Part_METIS();

    virtual idx_t get_epart( const int &ee ) const {return epart[ee];}

    // ------------------------------------------------------------------------
    // For multifield partition, we allow the access of partition results by
    // specifiying the field index, and nn ranges in 
    // [ 0 , mesh[field]->get_nFunc )
    // ------------------------------------------------------------------------
    virtual idx_t get_npart( const int &nn, const int &field = 0 ) const 
    {return npart[nn + field_offset[field]];}

    virtual bool get_isMETIS() const {return true;};

    virtual bool get_isDual() const {return isDual;};

    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}

    virtual bool is_serial() const {return false;}

  private:
    const bool isDual;
    const int dual_edge_ncommon;

    idx_t * epart, * npart;

    // ------------------------------------------------------------------------
    // For multi-field partition, we may list different field's mesh nodal index
    // in a single list of long nodal indices, starting from the first field's
    // indices. The IEN for METIS partitioning will be augmented by thinking of
    // an IEN that connect different fields' indices. The partitioined npart
    // will have the length equaling the sum of nFunc of all fields (or meshes).
    // To obtain each separate fields npart value (aka the CPU it belongs), we
    // store the field offset as the starting index in that field.
    // ------------------------------------------------------------------------
    std::vector<int> field_offset;

    virtual void write_part_hdf5( const std::string &fileName, 
        const idx_t * const &part_in,
        const int &part_size, const int &cpu_size ) const;

    // ------------------------------------------------------------------------
    // This function will write the data of part_in in 64bit HDF5 format. 
    // This function should be called if idx_t is the int64_t.
    // ------------------------------------------------------------------------
    virtual void write_part_hdf5_64bit( const std::string &fileName, 
        const int64_t * const &part_in,
        const int64_t &part_size, const int &cpu_size ) const;
};

#endif
