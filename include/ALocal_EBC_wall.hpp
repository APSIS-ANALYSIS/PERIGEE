#ifndef ALOCAL_EBC_WALL_HPP
#define ALOCAL_EBC_WALL_HPP
// ============================================================================
// ALocal_EBC_wall.hpp
//
// Analysis use: local subdomain's elemental boundary condition. This is a 
// derived class from ALocal_EBC to add the wall thickness and young's modulus. 
//
// Author: Ju Liu
// Date: Aug. 10 2020
// ============================================================================
#include "ALocal_EBC.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "HDF5_Writer.hpp"

class ALocal_EBC_wall : public ALocal_EBC
{
  public:
    ALocal_EBC_wall( const std::string &fileBaseName,
        const int &in_cpu_rank, const int &in_face_nqp,
        const std::string &gname = "ebc_wall",
        const std::string &in_ps_filebasename = "./prestress/prestress" );

    virtual ~ALocal_EBC_wall() = default;

    virtual void print_info() const;

    virtual int get_num_local_node_on_sur() const {return num_local_node_on_sur;}

    // ------------------------------------------------------------------------
    // ee is the local element index, ranges in [ 0, num_local_cell[0] )
    // elem_thickness is the output with length cell_nLocBas[0], which
    // represents the thickness at each local node
    // ------------------------------------------------------------------------
    virtual void get_thickness(const int &ee, double * const &elem_thickness) const;

    // ------------------------------------------------------------------------
    // elem_youngsmod is the output with length cell_nLocBas[0], which
    // represents the Young's modulus at each local node
    // ------------------------------------------------------------------------
    virtual void get_youngsmod(const int &ee, double * const &elem_youngsmod) const;

    // ------------------------------------------------------------------------
    // elem_springconst is the output with length cell_nLocBas[0], which
    // represents the spring constant at each local node
    // ------------------------------------------------------------------------
    virtual void get_springconst(const int &ee, double * const &elem_springconst) const;

    // ------------------------------------------------------------------------
    // elem_dampingconst is the output with length cell_nLocBas[0], which
    // represents the damping constant at each local node
    // ------------------------------------------------------------------------
    virtual void get_dampingconst(const int &ee, double * const &elem_dampingconst) const;

    virtual int get_local_node_on_sur_pos( const int &ii ) const 
    {return local_node_on_sur_pos[ii];}

    // ------------------------------------------------------------------------
    // elem_quaprestress is the output with length 6 x num_local_cell[0] x
    // num_face_nqp, which represents the prestress values at each quadrature
    // points. The data is arranged the following way,
    // elem_quaprestress[6*ii+jj] is the stress tensor jj-th component in Voigt
    // notation at the ii-th quadrature point.
    // ------------------------------------------------------------------------
    virtual void get_prestress(const int &ee, double * const &elem_quaprestress) const;

    // ------------------------------------------------------------------------
    // elem_quaprestress is the input with length 6 x num_face_nqp.
    // ------------------------------------------------------------------------
    virtual void set_prestress(const int &ee, const double * const &elem_quaprestress);

    virtual void write_prestress_hdf5() const;

  protected:
    // ------------------------------------------------------------------------
    // The rank or id of the subdomain
    // ------------------------------------------------------------------------
    const int cpu_rank;

    // ------------------------------------------------------------------------
    // num quadrature points for surface element
    // ------------------------------------------------------------------------
    const int face_nqp;

    // ------------------------------------------------------------------------
    // file name for storing the prestress values.
    // It is used in the write_prestress_hdf5 function
    // ------------------------------------------------------------------------
    const std::string ps_fileBaseName;

    // ------------------------------------------------------------------------
    // The number of local nodes belonging to this subdomain that is on the
    // surface
    // ------------------------------------------------------------------------
    int num_local_node_on_sur;

    // ------------------------------------------------------------------------
    // If this partition owns any part of the wall, these wall property
    // vectors are each of length num_local_node[0].
    // Otherwise, these vectors are of length 0.
    // ------------------------------------------------------------------------
    std::vector<double> thickness, youngsmod, springconst, dampingconst;

    // ------------------------------------------------------------------------
    // If num_local_node_on_sur > 0, this partition owns wall surface nodes, and
    // their position in the local_to_global array is stored in
    // local_node_on_sur_pos
    // Length is num_local_node_on_sur
    // ------------------------------------------------------------------------
    std::vector<int> local_node_on_sur_pos;

    // ------------------------------------------------------------------------
    // Prestress components (6 in Voigt notation) per quadrature point
    // Length 6 x face_nqp x num_local_cell[0]
    // ------------------------------------------------------------------------
    std::vector<double> qua_prestress;
};

#endif
