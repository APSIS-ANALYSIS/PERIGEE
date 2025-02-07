#ifndef ALOCAL_EBC_OUTFLOW_HPP
#define ALOCAL_EBC_OUTFLOW_HPP
// ============================================================================
// ALocal_EBC_outflow.hpp
//
// Analysis-use: local subdomain's elemental boundary condition. This is a 
// derived class from ALocal_EBC to add int_NA, the integral of the basis 
// functions over the outlet face, and the outlet's outward normal vector.
//
// Author: Ju Liu
// Date: Mar. 19 2019
// ============================================================================
#include "ALocal_EBC.hpp"

class ALocal_EBC_outflow : public ALocal_EBC
{
  public:
    ALocal_EBC_outflow( const std::string &fileBaseName,
        const int &cpu_rank, const std::string &gname="/ebc" );

    virtual ~ALocal_EBC_outflow() = default;

    virtual void print_info() const;

    // ------------------------------------------------------------------------
    // The following get functions will give meaningful data if the
    // num_local_cell[ii] > 0, that is, the partition owns the ii-th outlet
    // surface
    // Note, the get functions do NOT check the range of the ii index
    // ------------------------------------------------------------------------
    virtual int get_num_face_nodes(const int &ii) const {return num_face_nodes[ii];}

    virtual std::vector<double> get_intNA(const int &ii) const {return intNA[ii];}

    virtual std::vector<int> get_LID(const int &ii) const {return LID[ii];}

    virtual Vector_3 get_outvec( const int &ii ) const {return outvec[ii];}

  protected:
    // Length num_ebc.
    // If this partition owns any part of face xx, num_face_nodes[xx] records
    // the *total* number of face nodes, which is also the length of intNA[xx].
    // This times 3 gives the length of LID[xx]. If this partition does not
    // own any part of face xx, then num_face_nodes[xx] will be zero.
    std::vector<int> num_face_nodes;

    // Length num_ebc.
    // If this partition owns any part of face xx, intNA[xx] is of length
    // num_face_nodes[xx] and records the face integral of each basis function NA. 
    // Otherwise, intNA[xx] is of length 0. 
    std::vector< std::vector<double> > intNA;

    // Length num_ebc.
    // If this partition owns any part of face xx, LID[xx] is of length
    // 3 * num_face_nodes[xx] and stores each face node's LID values
    // (new numbering) corresponding to the 3 velocity dofs. Otherwise,
    // LID[xx] is of length 0.
    std::vector< std::vector<int> > LID;

    // Length num_ebc
    // It has x-y-z components of the unit outward normal vector
    std::vector< Vector_3 > outvec;
};

#endif
