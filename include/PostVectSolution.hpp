#ifndef POSTVECTSOLUTION_HPP
#define POSTVECTSOLUTION_HPP
// ============================================================================
// PostVectSolution.hpp
// ----------------------------------------------------------------------------
// This is a solution class that is used for postprocessing, including calculating 
// solution errors, generating vtk visualization files, etc. Notice, since the 
// solution vector generated from fem analysis code is usually in PETSc binary 
// format and reordered based on the node partition, this class has to read the 
// binary solutino and reorder the solution indices in order to make it compatible 
// with the partition of the postprocessing. Hence, there are three sets of numbering 
// involved: (a) the natural golbal numbering; (b) the analysis code's numbering; 
// (c) the postprocessing code's numbering.
//
// To do these, this code
// 1. read the PETSc vector into memory as a whole;
// 2. rearrange the numbering of the vector based on postproce code's node partition;
// 3. extract local vector for each processor in postprocessing.
//
// Author: Ju Liu
// Date: Dec 10 2013
// ============================================================================
#include "APart_Node.hpp"
#include "HDF5_Tools.hpp"
#include "petscvec.h"

class PostVectSolution
{
  public:
    // ------------------------------------------------------------------------
    // Constructor:
    //   \para solution_file_name: the name of PETSc binary file that
    //                             records the solution vector;
    //   \para analysis_node_mapping_file: the old_2_new and new_2_old 
    //                             mapping from analysis run;
    //   \para aNode_ptr: the partition of node information from
    //                    post_part files;
    //   \para input_dof: the degree of freedom of the input solution vector 
    // Output: Construct the vector that stores the solution vector
    //         which belongs to the local partition, i.e., the post
    //         processes' local_to_global indices. 
    // ------------------------------------------------------------------------
    PostVectSolution( const std::string &solution_file_name,
       const std::string &analysis_node_mapping_file,
       const std::string &post_node_mapping_file,
       const APart_Node * const &aNode_ptr,
       const int &in_nfunc, const int &input_dof );

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~PostVectSolution() = default;

    // ------------------------------------------------------------------------
    // Print the loc_solution vector
    // ------------------------------------------------------------------------
    void print_info() const;

    // ------------------------------------------------------------------------
    // get functions that give read permission to private data
    // ------------------------------------------------------------------------
    double get_locsol(const int &pos) const {return loc_solution[pos];}
    
    int get_dof() const {return dof_per_node;}
    
    int get_solsize() const {return loc_sol_size;}

    // ------------------------------------------------------------------------
    // get local element's solution
    //   \para field: the dof field that we need, field < dof_per_node
    //   \para nLocBas: the number of solution that needs to be passed out
    //   \para eien: the element ien array
    //   \para esol: the double array that passes the solution vector
    // Users are responsible for allocating and deleting eien & esol
    // ------------------------------------------------------------------------
    void get_esol( const int &field, const int &nLocBas, 
        const int * const &eien, double * const &esol) const;

  private:
    const int dof_per_node; // dof of the solution vector
    const int loc_sol_size; // nlocghonode x dof
    
    // loc_solution is the data structure that holds the solution vector of the
    // local postprocessing partition.
    std::vector<double> loc_solution;
    
    // ------------------------------------------------------------------------
    // ReadPETSc_vec: read a PETSc vector into memory as a double array
    //   \para vec_size: the total length of the PETSc vector. This is used
    //                   to check the correcness of the reading.
    //   \para veccopy: the double array that is used to allocate the solution
    //                  users are responsible for allocating and deleting this
    //                  dynamic array. This array should have length vec_size.
    // ------------------------------------------------------------------------
    void ReadPETSc_vec( const std::string &solution_file_name,
        const int &vec_size, double * const &veccopy );

    // ------------------------------------------------------------------------
    // ReadNodeMapping: reads the old_2_new or new_2_old array into the nodemap 
    //                  array. The users are responsible for allocating nFunc 
    //                  length for nodemap and deleting it after use.
    // \para node_mapping_file: the file that stores the mapping arrays
    // \para mapping_type: data_name in the file: new_2_old / old_2_new
    // \para node_size: the allocated length for nodemap
    // \para nodemap: the array that passes out the mapping.
    // ------------------------------------------------------------------------
    void ReadNodeMapping( const std::string &node_mapping_file,
        const char * const &mapping_type, const int &node_size,
        int * const &nodemap ) const;
};

#endif
