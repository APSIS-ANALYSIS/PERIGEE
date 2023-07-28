#ifndef GMSH_FILEIO_HPP
#define GMSH_FILEIO_HPP
// ==================================================================
// Gmsh_FileIO.hpp
//
// This is a class of tools that read and write Gmsh's .msh files.
//
// Date Created: July 1 2017
// Author: Ju Liu
// ==================================================================
#include "Tet_Tools.hpp"
#include "HDF5_Writer.hpp"

class Gmsh_FileIO
{
  public:
    Gmsh_FileIO( const std::string &in_file_name );

    ~Gmsh_FileIO();

    void print_info() const;

    int get_num_phy_domain() const {return num_phy_domain;}

    int get_num_phy_domain_3d() const {return num_phy_domain_3d;}

    int get_num_phy_domain_2d() const {return num_phy_domain_2d;}

    int get_num_phy_domain_1d() const {return num_phy_domain_1d;}

    int get_phy_id_3d( const int &ii ) const {return phy_3d_index[ii];}
    
    int get_phy_id_2d( const int &ii ) const {return phy_2d_index[ii];}
    
    int get_phy_id_1d( const int &ii ) const {return phy_1d_index[ii];}

    std::string get_phy_name_3d( const int &ii ) const {return phy_3d_name[ii];} 
    
    std::string get_phy_name_2d( const int &ii ) const {return phy_2d_name[ii];} 
    
    std::string get_phy_name_1d( const int &ii ) const {return phy_1d_name[ii];} 

    // --------------------------------------------------------------
    // In FSI problems, we require that the 3d physical domain index
    // 0 be the fluid domain and 3d physical domain index 1 be the 
    // solid domain. This function will check the ordering once the
    // Gmsh is constructed and throw a fatal error message if the 
    // ordering is wrong.
    // Sometimes, we use "lumen" as the name for fluid sub-domain,
    // and "tissue" as the name for the solid sub-domain, in the gmsh
    // file, so the user may put proper name in the argument.
    // --------------------------------------------------------------
    void check_FSI_ordering( const std::string &phy1="fluid",
       const std::string &phy2="solid" ) const;

    // --------------------------------------------------------------
    // In FSI problems, in addition to make the elements start from
    // the fluid sub-domain, we may also need to have the nodal
    // indices start from the fluid domain. This funciton will update
    // the nodal indices such that the fluid domain owns the first
    // a few nodal points, and the solid domain owns the next.
    // --------------------------------------------------------------
    void update_FSI_nodal_ordering();

    // --------------------------------------------------------------
    // update the IEN array to accomodate for the VTK ordering for
    // quadratic tetrahedral elements. For quadratic tetrahedral
    // element, the 8-th and 9-th nodes in Gmsh correspond to the 
    // 9-th and 8-th nodes in VTK format.
    // \para index_3d : the 3D domain index. so the 
    //                  eIEN[ phy_3d_index[index_3d]][ ]
    //                  will be modified.
    // --------------------------------------------------------------
    void update_quadratic_tet_IEN( const int &index_3d );

    // --------------------------------------------------------------
    // write a vtp file for an interior surface between two physical
    // subdomains, with name surfaceName_vol1Name_vol2Name.vtp
    // The user is responsible to give the correct indices for the
    // surface and the volume such that the surface is in between the
    // two volumetric domains.
    // \para index_sur : the index for the surface in phy_2d_index array
    // \para index_vol1 : the index for the volume in phy_3d_index array
    // \para index_vol2 : the index for the volume in phy_3d_index array
    // in the output vtp file, two face2elem mapping will be associated
    // with the face triangle element,
    // face2elem_1 is associated with the element in index_vol1 domain;
    // face2elem_2 is associated with the element in index_vol2 domain.
    // All element indices are global, which means phy_3d_start_index
    // is used to give the global domain element index for all output
    // elemental quantities.
    // --------------------------------------------------------------
    void write_interior_vtp( const int &index_sur,
       const int &index_vol1, const int &index_vol2 ) const;

    // --------------------------------------------------------------
    // write a vtp file for a surface associated with a volume mesh
    // with name surfacename_volumename.vtp
    // \para index_sur : the index for the surface in phy_2d_index array
    // \para index_vol : the index for the associated volume
    //                   in phy_3d_index array
    // The function will check the surface's triangle's associated
    // volume mesh index and write as its element index.
    // This function is similar to TET_T::tetgenio2vtp.
    // The face2elem writes the global element number that the faces
    // belong to, this means for many sub-volumetric domains, the
    // phy_3d_start_index is used to give the global whole domain element
    // index.
    // If isf2e is set to false, the face2elem will not be calcualted to
    // save time. In practice, we only need this for surfaces that has
    // boundary integral calculations.
    // In principle, we suggest that the surface belong to only one
    // volumetric physical domain. In case that a surface spans over
    // many physical volumetric domain, the face2elem mapping is -1.
    // --------------------------------------------------------------
    void write_vtp(const int &index_sur, const int &index_vol,
        const bool &isf2e = false) const;
  
    // --------------------------------------------------------------
    // write a vtu file for a surface associated with a volume mesh
    // with name surfacename_volumename.vtu
    // This function is used specifically for quadratic triangle mesh
    // and will write the surface mesh into a vtu file. It functionality
    // is quite close to write_vtp.
    // --------------------------------------------------------------
    void write_quadratic_sur_vtu( const int &index_sur, 
        const int &index_vol, const bool &isf2e = false) const;

    // --------------------------------------------------------------
    // Write a vtu file for all volumetric physical domain together.
    // a physical tag will be added to distinguish the physical problem
    // domain.
    // --------------------------------------------------------------
    void write_vtu( const std::string &in_fname, const bool &isXML ) const;

    // --------------------------------------------------------------
    // write a separate vtu file for each physical volumetric domain.
    // E.G. for a FSI problem a fluid and a solid domain will be written;
    //      for a fluid problem, a single fluid domain will be written.
    // --------------------------------------------------------------
    void write_each_vtu() const;
    
    // --------------------------------------------------------------
    // write a h5 file for the 2D simplex domain with given index : 
    // index_2d, and write the 2d domain's line boundaries specified by 
    // the index-list index_1d. Within this write function, we locate 
    // the associated 2d element with each 1d element and write is as an
    // attribute for the 1d element entry.
    // index_2d in [ 0, num_phy_domain_2d ),
    // index_1d in [ 0, num_phy_domain_1d ).
    // --------------------------------------------------------------
    void write_tri_h5(const int &index_2d, 
        const std::vector<int> &index_1d ) const;

    // --------------------------------------------------------------
    // write a h5 file for the 3D simplex domain with given index :
    // index_3d, and write the 3d domain's surface boundaries specified
    // by the index-list index_2d.
    // index_3d in [ 0, num_phy_domain_3d ),
    // index_2d in [ 0, num_phy_domain_2d ).
    // --------------------------------------------------------------
    void write_tet_h5( const int &index_3d,
        const std::vector<int> &index_2d ) const;

    // --------------------------------------------------------------
    // write a h5 file for the 3D simplex domain with given index,
    // and the index list of the faces that needs a face2element map.
    // index_3d in [ 0 , num_phy_domain_3d ),
    // index_2d in [ 0 , num_phy_domain_2d ),
    // index_2d_need_facemap is a subset of index_2d
    // The index in the index_2d_need_facemap will generate the
    // face2elem array;
    // otherwise, the face2elem array will be empty. 
    // --------------------------------------------------------------
    void write_tet_h5( const int &index_3d,
        const std::vector<int> &index_2d,
        const std::vector<int> &index_2d_need_facemap ) const;

  private:
    Gmsh_FileIO(){}; // disallow default constructor

    const std::string filename; // file = xxx.msh

    // --------------------------------------------------------------
    // physical info 
    // number of physical domains.
    int num_phy_domain;

    // details of the physical domain
    // the number of 3d, 2d, and 1d domains.
    // num_phy_domain = num_phy_domain_3d + num_phy_domain_2d
    //                 + num_phy_domain_1d
    int num_phy_domain_3d, num_phy_domain_2d, num_phy_domain_1d;

    // the indices (physical tag) of the 1d/2d/3d domains respectively
    // {phy_index} = {phy_3d_index} + {phy_2d_index} + {phy_1d_index} 
    std::vector<int> phy_3d_index, phy_2d_index, phy_1d_index;

    // the physicla subdomain names of the corresponding phy_xd_index.
    std::vector<std::string> phy_3d_name, phy_2d_name, phy_1d_name;

    // stores the number of 3d/2d/1d element respectively
    // vector lengths are num_phy_domain_3d/2d/1d respectively
    std::vector<int> phy_3d_nElem, phy_2d_nElem, phy_1d_nElem;

    // Stores the starting index for the 3d/2d volume mesh, with length
    // num_phy_domain_3d/num_phy_domain_2d.
    // This is a data generated within this code. It is utilized such
    // that the first phy domain starts with the first a few elements,
    // the second phy domain element indices follows, etc. This is because
    // Gmsh has an element id for 1d, 2d, and 3d elements all together,
    // we need to manage our own element id ourselves. 
    std::vector<int> phy_3d_start_index;
    std::vector<int> phy_2d_start_index;

    // --------------------------------------------------------------
    // Geometry info    
    int num_node; // number of nodal points
    
    std::vector<double> node; // 3 x num_node: x-y-z coordinates

    int num_elem; // number of total elem (2D and 3D together)

    // size is {num_phy_domain} x { ele_nlocbas[ii] times phy_domain_nElem[ii] }
    // note: the eIEN values are with nodal indices starting from 0.       
    //       the first argument 0 <= ii < num_phy_domain is the physical tag
    std::vector< std::vector<int> > eIEN;

    // stores the number of basis functions in each physical domain. 
    // This implicitly implies that the we use the same type of element 
    // in each physical domain. 
    // ele_nlocbas stores the number of local basis functions in element.
    // the argument 0<= ii < num_phy_domain is the physical tag
    std::vector<int> ele_nlocbas;

    // stores the element type for each physical subdomain.
    // the argument 0<= ii < num_phy_domain is the physical tag
    // its value and the corresponding ele_nlocbas satisfies the ele_nlocbas
    // mapping, which is defined by Gmsh.
    std::vector<int> ele_type;
    // --------------------------------------------------------------
};

#endif
