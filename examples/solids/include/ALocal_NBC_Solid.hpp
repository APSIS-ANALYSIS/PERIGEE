#ifndef ALOCAL_NBC_SOLID_HPP
#define ALOCAL_NBC_SOLID_HPP
// ============================================================================
// ALocal_NBC_Solid.hpp
//
// Extension of ALocal_NBC with displacement-driven tags.
//
// Date: Feb 02 2026
// ============================================================================
#include "ALocal_NBC.hpp"

class ALocal_NBC_Solid : public ALocal_NBC
{
  public:
    ALocal_NBC_Solid( const std::string &fileBaseName,
        const int &cpu_rank, const std::string &gname="/nbc" );

    ALocal_NBC_Solid( const HDF5_Reader * const &h5r,
        const std::string &gname="/nbc" );

    virtual ~ALocal_NBC_Solid() = default;

    virtual int get_LDN_is_disp_driven( const int &dof_index, const int &node ) const
    { return LDN_is_disp.empty() ? 0 : LDN_is_disp[LD_offset[dof_index] + node]; }

  private:
    std::vector<int> LDN_is_disp;

    void read_disp_flag( const HDF5_Reader * const &h5r,
        const std::string &gname );

    ALocal_NBC_Solid() = delete;
};

#endif
