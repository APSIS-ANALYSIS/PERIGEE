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
#include "Vis_Tools.hpp"

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
       const std::vector<int> &analysis_node_mapping,
       const std::vector<int> &post_node_mapping,
       const APart_Node * const &aNode_ptr,
       int input_dof );

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
    double get_locsol(int pos) const {return loc_solution[pos];}
    
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
    void get_esol( int field, int nLocBas, 
        const int * const &eien, double * const &esol) const;

  private:
    const int dof_per_node; // dof of the solution vector
    const int loc_sol_size; // nlocghonode x dof
    
    // loc_solution is the data structure that holds the solution vector of the
    // local postprocessing partition.
    std::vector<double> loc_solution;
};

#endif
