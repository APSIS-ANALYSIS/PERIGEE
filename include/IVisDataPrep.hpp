#ifndef IVISDATAPREP_HPP
#define IVISDATAPREP_HPP
// ==================================================================
// IVisDataPrep.hpp
// Interface for visualization data preparation.
//
// This is an abstract class that includes two main functions that we
// need to handle visualization data: generating the types of data
// that we want to visualize, specify the array length and name; then
// read in the binary solution and obtain the correcponding vectors.
//
// Date: Dec. 17 2013
// ==================================================================
#include "PostVectSolution.hpp"

class IVisDataPrep
{
  public:
    // -------------------------------------------------------------------
    // ! empty constructor
    //   initialization for the private members should be done in the 
    //   constructor of inherited classes
    //   This Visualization data preparation class should only define
    //   the number of types of data that needs to be visualized, i.e.,
    //   arrayCompSize; their names, i.e., arrayNames; and their length,
    //   i.e., arraySizes.
    //   The one extra function, get_pointArray, is the reader that reads
    //   in petsc binary vector and transforms it into the pointArrays,
    //   which corresponds to the data named arrayNames.
    // -------------------------------------------------------------------
    IVisDataPrep() = default;

    virtual ~IVisDataPrep() = default;
  
    // -------------------------------------------------------------------
    // ! output: pointArrays, should have arrayCompSize for rows, 
    //   and row ii has length arraySizes[ii] * nlocghonode, 
    //   with name arrayNames[ii].
    //   
    //   Users are responsible for allocating and deleting the memory
    //   for pointArrays!
    //
    // ! input:
    // \para solution_file_name: the petsc binary solution vector file;
    // \para analysis_node_mapping_file: the old_2_new & new_2_old
    //                                   mapping from analysis run;
    // \para post_node_mapping_file: the old_2_new & new_2_old mapping
    //                            from the preprocess of postprocessor;
    // \para APart_Node: the node partition info;
    // \para input_nfunc: number of basis functions;
    // \para input_dof: degree of freedom of this problem.
    // -------------------------------------------------------------------
    virtual void get_pointArray(
        const std::string solution_file_name,
        const std::vector<int> &analysis_node_mapping,
        const std::vector<int> &post_node_mapping,
        const APart_Node * const &nNode_ptr,
        const int &input_nfunc,
        const int &input_dof,
        double ** &pointArrays ) const
    {
      SYS_T::print_fatal("Warning: get_pointArray is not implemented.\n");
    }

    // -------------------------------------------------------------------
    // ! input gives a list of solution file names to read.
    //   This is used for problems where we store the solutions separately.
    //   E.g. we may store the displacement, pressure, and velocity as
    //   three petsc binary vectors on disk. 
    // -------------------------------------------------------------------
    virtual void get_pointArray(
        const std::vector<std::string> solution_file_names,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const int &in_nfunc,
        double ** &pointArrays ) const
    {
      SYS_T::print_fatal("Warning: get_pointArray is not implemented.\n");
    }

    // ------------------------------------------------------------------------
    // ! input gives three solution names for the disp, pres, and velo
    // respectively. The kinematic variables (i.e. disp and velo) were
    // partitioned in a different manner from that of the pressure variable.
    // Thus, there are two sets of node mapping files.
    // ------------------------------------------------------------------------
    virtual void get_pointArray(
        const std::string &disp_solution_file_name,
        const std::string &pres_solution_file_name,
        const std::string &velo_solution_file_name,
        const std::string &an_v_node_mapping_file,
        const std::string &an_p_node_mapping_file,
        const std::string &pn_v_node_mapping_file,
        const std::string &pn_p_node_mapping_file,
        const APart_Node * const &pNode_v,
        const APart_Node * const &pNode_p,
        const int &input_nfunc_v,
        const int &input_nfunc_p,
        double ** &pointArrays ) const
    {
      SYS_T::print_fatal("Warning: get_pointArray is not implemented.\n");
    }

    // -------------------------------------------------------------------
    // ! give read permission for the private data: arrayCompSize
    //   i.e. the number of quantities to be visualized.
    // -------------------------------------------------------------------
    int get_arrayCompSize() const {return arrayCompSize;}

    // -------------------------------------------------------------------
    // ! give the point array size.
    //   In default, we assume that the visualization is for each input
    //   point array, then ptarray_size = arrayCompSize.
    //   On the other side, if we input diplacement vector, and wish to
    //   visualize displacement as well as det(F), then the ptarray_size
    //   is the input point array's size, aka, 1 (displacement), and the
    //   arrayCompSize = 2, aka, displacement vector with length 3 and 
    //   J, with length 1.
    // -------------------------------------------------------------------
    virtual int get_ptarray_size() const 
    {
      SYS_T::commPrint("Warning: ptarray_size is defaulted to arrayCompSize. \n");
      return arrayCompSize;
    }

    // -------------------------------------------------------------------
    // ! Return the ii's quantity's name
    // -------------------------------------------------------------------
    std::string get_arrayNames(const int &ii) const
    {
      if(ii>=arrayCompSize)
        SYS_T::print_fatal("Error: get_arrayNames argument is out of range. \n");
      
      return arrayNames[ii];
    }

    // -------------------------------------------------------------------
    // ! Return the ii's quantity's size: 
    //        scalar 1,    2D vector 2,    3D vector 3
    // -------------------------------------------------------------------
    int get_arraySizes(const int &ii) const
    {
      if(ii>=arrayCompSize)
        SYS_T::print_fatal("Error: get_arraySizes argument is out of range. \n");
      
      return arraySizes[ii];
    }

    // -------------------------------------------------------------------
    // ! Return the point array's ii-th component's length.
    //   0 <= ii < get_ptarray_size().
    //   E.g. if input is only a 3D displacement vector,
    // -------------------------------------------------------------------
    virtual int get_ptarray_comp_length( const int &ii ) const
    {
      if(ii>=get_ptarray_size())
        SYS_T::print_fatal("Error: get_ptarray_comp_length argument is out of range. \n");
      
      SYS_T::commPrint("Warning: get_ptarray_comp_length is defaulted to get_arraySizes. \n");

      return get_arraySizes(ii);
    }


    // -------------------------------------------------------------------
    // ! print_info:
    //   print the visualization data arrangement on screen
    // -------------------------------------------------------------------
    virtual void print_info() const
    {
      PetscPrintf(PETSC_COMM_WORLD, "======================================= \n");
      PetscPrintf(PETSC_COMM_WORLD, "Data to be visualized: \n");
      PetscPrintf(PETSC_COMM_WORLD, "-- %d set(s) of data:\n", get_arrayCompSize());
      for(int ii=0; ii<get_arrayCompSize(); ++ii)
      {
        PetscPrintf(PETSC_COMM_WORLD, "--- %s \t", get_arrayNames(ii).c_str() );
        PetscPrintf(PETSC_COMM_WORLD, "with size %d \n", get_arraySizes(ii)  );
      }
      PetscPrintf(PETSC_COMM_WORLD, "-- %d set(s) data to be read, with length: \n", 
          get_ptarray_size());
      for(int ii=0; ii<get_ptarray_size(); ++ii)
        PetscPrintf(PETSC_COMM_WORLD, "--- Entry %d : length %d \n", 
            ii, get_ptarray_comp_length(ii));
      PetscPrintf(PETSC_COMM_WORLD, "======================================= \n");
    }

  protected:
    // -------------------------------------------------------------------
    // ! the number of types of data we need to visualize
    // -------------------------------------------------------------------
    int arrayCompSize {};

    // -------------------------------------------------------------------
    // ! arrayNames: the names for each type data to be visualized
    //   length : 0 <= ii < arrayCompSize
    // -------------------------------------------------------------------
    std::vector<std::string> arrayNames {};

    // -------------------------------------------------------------------
    // ! arraySizes, should have length of arrayCompSize, it stores
    //   the size of each type data, e.g., 1 is scalar data, 3 is 
    //   3D vector data, etc.
    // -------------------------------------------------------------------
    std::vector<int> arraySizes {};
};

#endif
