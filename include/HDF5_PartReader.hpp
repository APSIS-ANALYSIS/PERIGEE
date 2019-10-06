#ifndef HDF5_PARTREADER_HPP
#define HDF5_PARTREADER_HPP
// ==================================================================
// HDF5_PartReader.hpp
// This is an interface for various reading function for the partition
// files stored in hdf5 format.
// ==================================================================
#include <cassert>
#include "Vec_Tools.hpp"

class HDF5_PartReader
{
  public:
    HDF5_PartReader( std::string partFileBaseName, const int &rank );
    virtual ~HDF5_PartReader();

    // EXT: extraction operator
    void get_EXT_x( std::vector<double> &ext_x ) const;
    void get_EXT_y( std::vector<double> &ext_y ) const;
    void get_EXT_z( std::vector<double> &ext_z ) const;
   
    // ------------------------------------------------------------------------
    // get_EXT_TS_full:
    // use the h5 part file reader to get the extraction operator for
    // T-Splines. The ext in part h5 file is a one-dimensional array for the 
    // rows of each elements. 
    // ------------------------------------------------------------------------
    void get_EXT_TS_full( std::vector<double> &ext ) const;
    
    void get_EXT_elemType( int &elemType ) const; 

    // GMI : group name: Global-Mesh-Info
    void get_GMI_degree( int &sDegree, int &tDegree, int &uDegree ) const;
    void get_GMI_degree( int &sDegree, int &tDegree ) const;
    void get_GMI_h_max( double &hx_max, double &hy_max, double &hz_max) const;
    void get_GMI_h_max( double &hx_max, double &hy_max) const;
    void get_GMI_h_min( double &hx_min, double &hy_min, double &hz_min) const;
    void get_GMI_h_min( double &hx_min, double &hy_min) const;
    void get_GMI_nElem( int &nElem, int &nElem_x, int &nElem_y, int &nElem_z ) const;
    void get_GMI_nElem( int &nElem, int &nElem_x, int &nElem_y ) const;
    void get_GMI_nElem( int &nElem ) const;
    void get_GMI_nFunc( int &nFunc, int &nFunc_x, int &nFunc_y, int &nFunc_z ) const;
    void get_GMI_nFunc( int &nFunc, int &nFunc_x, int &nFunc_y ) const;
    void get_GMI_nFunc( int &nFunc ) const;
    void get_GMI_nLocBas( int &nLocBas ) const;
    void get_GMI_nLocBas( std::vector<int> &nLocBas ) const;
    void get_GMI_probDim( int &probDim ) const;
    void get_GMI_dofNum( int &dofNum ) const;
    void get_GMI_dofMat( int &dofMat ) const;

    // LIEN
    void get_LIEN( const int &e, std::vector<int> &LIEN ) const;
    void get_LIEN( std::vector<int> &LIEN ) const;
    void get_1D_LIEN( std::vector<int> &LIEN ) const;

    // Local Element 
    void get_LE( std::vector<int> &elem_loc, int &nlocalele ) const;
    
    // hx hy hz
    void get_hxyz( std::vector<double> &hx, std::vector<double> &hy, std::vector<double> &hz ) const;
    void get_hxy( std::vector<double> &hx, std::vector<double> &hy ) const;

    // Local Node
    void get_LN( int &nlocalnode, int &nghostnode, int &nbadnode, 
        int &nlocghonode, int &ntotalnode, std::vector<int> &local_to_global,
        std::vector<int> &node_ghost, std::vector<int> &node_loc,
        std::vector<int> &node_loc_original ) const;

    void get_L2GN( int &nlocghonode, std::vector<int> &local_to_global ) const;
    
    // PI: group_name: Part_Info
    void get_PI( int &cpu_rank, int &cpu_size, int &dual_edge_ncommon ) const;

    // CPL: ControlPointsLocal
    void get_CPL( std::vector<double> &ctrlPts_x_loc, std::vector<double> &ctrlPts_y_loc,
        std::vector<double> &ctrlPts_z_loc, std::vector<double> &ctrlPts_w_loc ) const;

    
    // BC: group_name: bc
    //     return the local bc info.
    void get_BC_LID_dof( std::vector<int> &LID, int &dof ) const;
    void get_BC_LD(std::vector<int> &LDN, std::vector<int> &Num_LD) const;
    void get_BC_LP(std::vector<int> &LPSN, std::vector<int> &LPMN, std::vector<int> &Num_LP) const;
    void get_BC_LBCE( std::vector<int> &Num_LBCElem, std::vector<int> &LFront_Elem,
        std::vector<int> &LBack_Elem, std::vector<int> &LLeft_Elem,
        std::vector<int> &LRight_Elem, std::vector<int> &LTop_Elem,
        std::vector<int> &LBottom_Elem ) const;


    // NBC: group name: nbc
    //      get the nodal bc indices.
    void get_NBC_LID_dof( std::vector<int> &LID, int &dof ) const;
    void get_NBC_LD(std::vector<int> &LDN, std::vector<int> &Num_LD) const;
    void get_NBC_LPS(std::vector<int> &LPSN, std::vector<int> &LPMN, std::vector<int> &Num_LPS) const;
    void get_NBC_LPM(std::vector<int> &LocalMaster, std::vector<int> &LocalMasterSlave,
        std::vector<int> &Num_LPM) const;


    // EBCL group name: ebc
    //      get the elemental bc indices
    void get_EBC( std::vector<int> &num_lbcelem, std::vector<int> &lfront_elem,
      std::vector<int> &lback_elem, std::vector<int> &lleft_elem,
      std::vector<int> &lright_elem, std::vector<int> &ltop_elem,
      std::vector<int> &lbottom_elem  ) const;


  private:
    hid_t file_id;

    // Failed hdf5 functions returns a negative value, while successful
    // functions return a non-negative integer.
    void check_error(const herr_t &status, const char * const &funname ) const
    {
      if(status < 0)
      {
        std::cerr<<"Error: HDF5_PartReader::"<<funname<<std::endl;
        exit(EXIT_FAILURE);
      }
    }

    // Note: The following four functions can read in integer, double
    // datas as a whole from the disk, or can read in integer, double
    // row array in a two dim array in disk.
    //
    // The user is responsible to provide a int / double pointer, and 
    // delete the pointer to free memory after usage.
    int get_intScalar( const char * group_name, const char * data_name ) const;

    void get_intData( const char * group_name, const char * data_name, 
        hid_t &data_rank, hsize_t * &data_dims, int * &data ) const;

    void get_doubleData( const char * group_name, const char * data_name, 
        hid_t &data_rank, hsize_t * &data_dims, double * &data ) const;

    void get_intRowData( const char * group_name, const char * data_name,
        hid_t &data_rank, hsize_t * &data_dims, const int &row,
        int * &row_data ) const;

    void get_doubleRowData( const char * group_name, const char * data_name,
        hid_t &data_rank, hsize_t * &data_dims, const int &row,
        double * &row_data ) const;

    // CPL: group name: ctrlPts_loc
    //      return the local control points:
    //      ctrlPts_w_loc, ctrlPts_x_loc, ctrlPts_y_loc, ctrlPts_z_loc
    //      user is responsible for deleting the memory allocations.
    void get_CPL( double * &ctrlPts_x_loc, double * &ctrlPts_y_loc,
        double * &ctrlPts_z_loc, double * &ctrlPts_w_loc, int &len ) const;

    // LIEN: group name: LIEN
    //       user is responsible for deleting LIEN memory allocation.
    void get_LIEN( const int &e, int * &LIEN, int &length ) const;

    // LE: group name: Local_Elem
    //     return elem_loc array and nlocalele in one function.
    //     user is responsible for deleting memory allocations.
    void get_LE( int * &elem_loc, int &nlocalele ) const;
    
    // LN: group name: Local_Node
    //     return the data in Local_Node group, including:
    //     int scalar: nlocalnode, nghostnode, nbadnode, nlocghonode, 
    //                 ntotalnode.
    //     int 1-D array: local_to_global, node_ghost, node_loc,
    //                    node_loc_original
    //     user is responsible for freeing the memory allocations.
    void get_LN( int &nlocalnode, int &nghostnode, int &nbadnode, 
        int &nlocghonode, int &ntotalnode, int * &local_to_global,
       int * &node_ghost, int * &node_loc, int * &node_loc_original ) const;
};
#endif
