#ifndef ELEMBC_HPP
#define ELEMBC_HPP
// ==================================================================
// ElemBC.hpp
// 
// This is the abstract template header definition for the data structure
// recording the surface integral geometry and necessary mapping from the
// 2D/1D geometry back to the volumetric 3D/2D global data.
//
// The following data structure are necessary.
//
// 0: The number of different boundary surface domains;
// 1. The number of nodes on the Neumann boundary surface;
// 2. The number of cell on this surface (here, cell is typically 2d);
// 3. The x-y-z coordinates of the nodal points;
// 4. The connectivity array for the cell (here, the ien is for the 2d
//    cell, rather than the 3d cell.)
// 5. The mapping for the local nodes to the volume mesh node index;
// 6. The mapping from the 2d cell to the 3d cell.
//
// Note: 
// 1. we presumbly may have multiple different surface domains for
//    boundary treatment. So we we have the 0-th data in this class.
//    Alternatively, one can have all surface bc glued together and
//    differentiate them by the prescribed traction functions.
// 
// 2. The number of ebc type has to match the number of boundary 
//    prescribed data. (For example, in the previous NURBS code,
//    the number of ebc type is always 6, corresponding to the 6 faces
//    of the structured cube mesh. And we have 6 boundary integral
//    routine in the assembly code.)
//
// Date: Jan. 10 2017
// ==================================================================
#include "Hex_Tools.hpp"
#include "IIEN.hpp"
#include "FEType.hpp"

class ElemBC
{
  public:
    ElemBC() = default;

    virtual ~ElemBC() = default;

    // This returns the number of surface domains that requires boundary
    // integral.
    virtual int get_num_ebc() const
    {SYS_T::commPrint("Warning: get_num_ebc() is not implemented. \n"); return 0;}
    
    // This returns the number of nodes on the surface with ebc_id, wherein
    // 0 <= ebc_id < get_num_ebc() 
    virtual int get_num_node(const int &ebc_id) const
    {SYS_T::commPrint("Warning: get_num_node is not implemented. \n"); return 0;}
    
    // This returns the number of cells on the surface with ebc_id, wherein
    // 0 <= ebc_id < get_num_ebc()
    virtual int get_num_cell(const int &ebc_id) const
    {SYS_T::commPrint("Warning: get_num_cell is not implemented. \n"); return 0;}

    // This returns the number of local nodes in the boundary cell. For example,
    // linear triangles return 3;
    // bilinear quad returns 4;
    // quadratic triangles return 6.
    virtual int get_cell_nLocBas(const int &ebc_id) const
    {SYS_T::commPrint("Warning: get_cell_nLocBas is not implemented. \n"); return 0;}
    
    // This returns the coordinates of the nodal points, wherein
    // 0 <= ebc_id < get_num_ebc();
    // 0 <= node   < get_num_node(ebc_id);
    // 0 <= dir < 3.
    // dir = 0 => x direction
    // dir = 1 => y direction
    // dir = 2 => z direction
    virtual double get_pt_xyz(const int &ebc_id, const int &node, 
        const int &dir) const
    {SYS_T::commPrint("Warning: get_pt_xyz is not implemented. \n"); return 0.0;}

    // This returns the connectivity for the surface cell with the node list
    // It can be regarded as a IEN array for lower-dimensional elements.
    // 0 <= ebc_id < get_num_ebc();
    // 0 <= cell   < get_num_cell();
    // 0 <= lnode  < get_cell_nLocBas().
    virtual int get_ien(const int &ebc_id, const int &cell, const int &lnode) const
    {SYS_T::commPrint("Warning: get_ien is not implemented. \n"); return 0;}

    // This returns the global node index (i.e., the node indices in the
    // volumetric mesh) for nodes.
    // 0 <= ebc_id < get_num_ebc();
    // 0 <= node   < get_num_node(ebc_id);
    virtual int get_global_node(const int &ebc_id, const int &node_index) const
    {SYS_T::commPrint("Warning: get_global_node is not implemented. \n"); return 0;}

    virtual void get_global_node(const int &ebc_id,std::vector<int> &out) const
    {SYS_T::commPrint("Warning: get_global_node is not implemented. \n");}

    virtual std::vector<int> get_global_node(const int &ebc_id) const
    {SYS_T::commPrint("Warning: get_global_node is not implemented. \n"); return {};}

    // This returns the global volumetric element index for the surface cells.
    virtual int get_global_cell(const int &ebc_id, const int &cell_index) const
    {SYS_T::commPrint("Warning: get_global_cell is not implemented. \n"); return 0;}

    // This returns the cell's interior point's xyz coordinates.
    virtual double get_intpt_xyz(const int &ebc_id, const int &cell_index,
        const int &dir) const
    {SYS_T::commPrint("Warning: get_intpt_xyz is not implemented. \n"); return 0.0;}

    // This function returns the outward normal vector for faces
    virtual Vector_3 get_normal_vec( const int &ebc_id ) const
    {SYS_T::commPrint("Warning: get_normal_vec is not implemented. \n"); return Vector_3();}

    // This function returns the integral of basis NA on the faces
    virtual std::vector<double> get_intNA( const int &ebc_id ) const
    {SYS_T::commPrint("Warning: get_intNA is not implemented.\n"); return {};}

    // Access the data in ElemBC_3D_wall_turbulence, wall model type
    virtual int get_wall_model_type() const
    {SYS_T::commPrint("Warning: get_wall_model_type is not implemented. \n"); return -1;}

    // Access the data in ElemBC_3D_wall_turbulence, face id of volume element
    virtual int get_faceID( const int &cell_index ) const
    {SYS_T::commPrint("Warning: get_faceID is not implemented. \n"); return {};}
 
    // Overwrite ElemBC_3D_wall properties from a vtp/vtu file
    virtual void overwrite_from_vtk( const std::string &wallprop_vtk,
        const int &type, const std::string &vtk_fieldname )
    {SYS_T::commPrint("Warning: overwrite_from_vtk is not implemented. \n");}

    // write the surface to a vtk/vtu file
    virtual void write_vtk( const int &ebc_id,
       const std::string &filename="elembc_surface" ) const
    {SYS_T::commPrint("Warning: write_vtk is not implemented. \n");}

    // print the information of the ebc object on screen.    
    virtual void print_info() const
    {SYS_T::commPrint("Warning: print_info is not implemented. \n");}

    virtual void resetSurIEN_outwardnormal( const IIEN * const &VIEN )
    {SYS_T::print_fatal("Warning: resetSurIEN_outwardnormal is not implemented. \n");}
};

#endif
