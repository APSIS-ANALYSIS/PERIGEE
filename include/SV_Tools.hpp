#ifndef SV_TOOLS_HPP
#define SV_TOOLS_HPP
// ==================================================================
// SV_Tools.hpp
//
// This is a suite of tools that handle the IO files between this
// code and the SimVascular files.
//
// Author: Ju Liu
// Date Created: Aug. 6 2017
// ==================================================================
#include "Tet_Tools.hpp"

namespace SV_T
{
  // ================================================================
  // ===> 1. The first suite of tools are designed to read volumetric
  //         and surface mesh files and find the starting index for
  //         nodes and elements, then correct the index to make them
  //         start from zero. 
  // ================================================================
  // ----------------------------------------------------------------  
  // ! update_sv_vtu : read the volumetric mesh file and find the 
  //                   smallest index for nodes and elements. If they
  //                   are nonzero, update them to make the starting
  //                   index as zero. Output the simvascular's node
  //                   and element's starting index.
  //                   The new file will be write in the writename.
  //   note: in simvascular, the starting index is often 1.
  // ----------------------------------------------------------------  
  void update_sv_vtu( const std::string &filename,
      const std::string &writename,
      int &sv_node_start, int &sv_elem_start );

  // ----------------------------------------------------------------
  // ! gen_sv_fsi_vtus : merge fluid and solid vtu files into one 
  //                     single file for the whole domain, and write 
  //                     the SOLID mesh with the updated nodal and 
  //                     element indices. Since SV will generate fluid
  //                     and solid domains separately. In the merging,
  //                     we take fluid element and nodal indices first,
  //                     meaning their values will not be changed. Then
  //                     the solid domain element and nodal indices 
  //                     follow those of the fluid domain, meaning we
  //                     will adjust the nodes in the solid domain, adjust
  //                     the solid IEN array, and append the solid IEN
  //                     array after the fluid IEN array to generate
  //                     the new whole domain IEN array.
  //  \para filename_f      : the original fluid mesh vtu
  //  \para filename_s      : the original solid mesh vtu
  //  \para filename_f_wall : the UPDATED fluid wall vtp file
  //  \para writename       : the file name for the FSI whole domain
  //  \para wirtename_solid : the file name for the solid subdomain
  //  
  //  \output map_s_node : the mapping for the solid node to the FSI
  //                       indices.
  //  \output map_s_elem : the mapping for the solid element to the
  //                       FSI indices.
  // ----------------------------------------------------------------
  void gen_sv_fsi_vtus( const std::string &filename_f,
      const std::string &filename_s,
      const std::string &filename_f_wall,
      const std::string &writename, 
      const std::string &writename_solid,
      std::vector<int> &map_s_node, std::vector<int> &map_s_elem ); 

  // ---------------------------------------------------------------- 
  // ! update_sv_vtp : read the surface mesh file and update the nodal
  //                   and elemental indices using the starting index
  //                   obtained from uodate_sv_vtu, the volumetric
  //                   update function.
  //                   The new file will be written in writename.
  // ---------------------------------------------------------------- 
  void update_sv_vtp( const std::string &filename,
      const std::string &writename,
      const int &nstart, const int &estart );

  // ---------------------------------------------------------------- 
  // ! update_sv_sur_vtu : read the surface mesh written in a vtu file
  //                       and update the nodal & elemental indices.
  //                       This function is designed for quadratic
  //                       triangle mesh primarily.
  // ---------------------------------------------------------------- 
  void update_sv_sur_vtu( const std::string &filename,
      const std::string &writename,
      const int &nstart, const int &estart );

  // ---------------------------------------------------------------- 
  // update the vtp file with a node map and element map to map
  // the nodal and elemental indices to new values.
  // \para filename : the original SV generated solid surface vtp file.
  // \para writename : the name to be written
  // \para nmap : the nodal mapping from the solid alone index to FSI 
  //              index
  // \para emap : the elemental mapping from the solid alone index to 
  //              FSI index
  // ---------------------------------------------------------------- 
  void update_sv_vtp( const std::string &filename,
      const std::string &writename,
      const int &nstart, const int &estart,
      const std::vector<int> &nmap, const std::vector<int> &emap );

  // ----------------------------------------------------------------
  // ! compare_sv_vtp : compare two vtp and make sure they are matched
  //                    in the sense that they contain the same number
  //                    of nodes and the nodal xyz coordinates match,
  //                    and the element IEN match.
  //                    It is used to make sure in the FSI mesh the
  //                    lumen wall and the tissue inner surface are
  //                    from the same geometry.
  // ----------------------------------------------------------------
  void compare_sv_vtp( const std::string &filename_1,
      const std::string &filename_2 );

  // ----------------------------------------------------------------
  // ! find_idx : locate the index of the coordinates in the point
  //              array; otherwise, return -1.
  //   \para len : the number of points (pt.size() / 3)
  // ----------------------------------------------------------------
  int find_idx( const std::vector<double> &pt, const int &len, 
      const double &x, const double &y, const double &z,
      const double tol = 1.0e-15 );
}

#endif
