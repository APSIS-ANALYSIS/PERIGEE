#ifndef ALOCAL_EBC_WINTPTS_HPP
#define ALOCAL_EBC_WINTPTS_HPP
// ==================================================================
// ALocal_EBC_wIntPts.hpp
//
// Analysis use local subdomain's elemental boundary condition. This
// is derived from ALocal_EBC.hpp to add the interior point coordinates
// to the class.
// 
// The data is the vector local_intpts with variable length 
// num_ebc x num_local_cell[ii]
//
// Author: Ju Liu
// Date: Nov. 22 2017
// ==================================================================
#include "ALocal_EBC.hpp"

class ALocal_EBC_wIntPts : public ALocal_EBC
{
  public:
    ALocal_EBC_wIntPts( const std::string &fileBaseName,
        const int &cpu_rank,
        const std::string &gname="/ebc" );

    ~ALocal_EBC_wIntPts();

    virtual void get_intPts_xyz(const int &ii, const int &eindex, 
        double &coor_x, double &coor_y, double &coor_z ) const
    {
      coor_x = local_intpts[ii][3*eindex];
      coor_y = local_intpts[ii][3*eindex+1];
      coor_z = local_intpts[ii][3*eindex+2];
    }

  private:
    std::vector< std::vector<double> > local_intpts;
};

#endif
