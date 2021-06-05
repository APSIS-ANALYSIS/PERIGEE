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
#include "IGlobal_Part.hpp"
#include "IMesh.hpp"
#include "hdf5.h"

class Global_Part_Serial : public IGlobal_Part
{
  public:
    Global_Part_Serial( const IMesh * const &mesh,
       const char * const &element_part_name,
       const char * const &node_part_name );

    virtual ~Global_Part_Serial();

    virtual idx_t get_epart( const int &ee ) const {return epart[ee];}
    
    virtual idx_t get_npart( const int &nn ) const {return npart[nn];}

    virtual bool get_isMETIS() const {return isMETIS;};
    
    virtual bool get_isDual() const {return isDual;};
    
    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}
    
    virtual bool is_serial() const {true;}

  private:
    const bool isMETIS, isDual;
    const int dual_edge_ncommon;
    
    idx_t * epart;
    idx_t * npart;

    virtual void write_part_hdf5( const char * const &fileName, 
        const idx_t * const &part_in,
        const int &part_size, const int &cpu_size,
        const bool &part_isdual, const int &in_ncommon,
        const bool &isMETIS ) const;
};

#endif
