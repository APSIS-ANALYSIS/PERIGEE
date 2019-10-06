#ifndef ALOCAL_EBC_OUTFLOW_HPP
#define ALOCAL_EBC_OUTFLOW_HPP
// ==================================================================
// ALocal_EBC_outflow.hpp
//
// Analysis use local subdomain's elemental boundary condition. This
// is a derived class from ALocal_EBC to add int_NA, the face integral
// of the basis functions, and outward normal vector of the face.
//
// Author: Ju Liu
// Date: Mar. 19 2019
// ==================================================================
#include "ALocal_EBC.hpp"

class ALocal_EBC_outflow : public ALocal_EBC
{
  public:
    // Constructor will read the input lumped parameters from a file
    // lpn_file.txt is the default name.
    // It is allowed to use # to comment lines in the file.
    ALocal_EBC_outflow( const std::string &fileBaseName,
        const int &cpu_rank, const std::string &gname="/ebc" );

    virtual ~ALocal_EBC_outflow();

    virtual void print_info() const;

    virtual int get_num_face_nodes(const int &ii) const
    {return num_face_nodes[ii];}

    virtual void get_intNA(const int &ii, std::vector<double> &out) const
    {out = intNA[ii];}

    virtual void get_LID(const int &ii, std::vector<int> &out) const
    {out = LID[ii];}

    virtual void get_outvec( const int &ii, double &nx, double &ny,
        double &nz ) const
    {
      nx = outvec[ii][0];
      ny = outvec[ii][1];
      nz = outvec[ii][2];
    }

  protected:
    // If this partition owns the face, it will record the number of
    // face nodes, which is also the length of intNA[xx], xx is the
    // face id. It times 3 gives the length of LID[xx].
    // If this partition does not own the face, it will be zero.
    std::vector<int> num_face_nodes;

    // local copy of the integral of NA basis on the face
    // length is num_ebc, and intNA[ii] length is 0 if this partition
    // does not own the face; or the number of nodes on the whole 
    // face if it owns any cell on the face.
    std::vector< std::vector<double> > intNA;

    // The corresponding LID value for the nodes associated with intNA.
    // Length is num_ebc and LID[ii] length is 0 if this partition does
    // not own the face; or the number of nodes on the whole face x 3
    // if it owns any cell on the face.
    std::vector< std::vector<int> > LID;

    std::vector< std::vector<double> > outvec;
};

#endif
