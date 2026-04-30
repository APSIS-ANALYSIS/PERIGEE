#ifndef PDNSOLUTION_HPP
#define PDNSOLUTION_HPP
// ============================================================================
// PDNSolution.hpp
//
// This is a base class for my finite element solutions.
//
// Main functions: 1. save the solution on disk;
//                 2. perform linear operations on the solution vector;
//                 3. obtain the local element solution efficiently;
//                 4. print the solution on screen.
// 
// Data: public : Vec the solution
//       private: int nlocal, nghost, the number of local and ghost
//                nodes in the subdomain.
//                int dof_num, the number of degrees of freedom per
//                node, which could be different from dof in APart_Node.
//
// Author: Ju Liu
// Date: Nov. 23th 2013
// ============================================================================
#include "Sys_Tools.hpp"
#include "APart_Node.hpp"

class PDNSolution
{
  public:
    Vec solution; 
    
    // ------------------------------------------------------------------------
    // ! Construct a solution vec compatible with the analysis node partition,
    //   with the solution vector's dof value being equal to the APart_Node
    //   class's dof value.
    // ------------------------------------------------------------------------
    PDNSolution( const APart_Node * const &pNode );

    // ------------------------------------------------------------------------
    // ! Construct a solution vec compatible with the analysis node partition, 
    //   but with a different dof_num from pNode->dof. 
    //   The users specify a dof number for the solution class.
    // ------------------------------------------------------------------------
    PDNSolution( const APart_Node * const &pNode, int input_dof_num );
    
    // ------------------------------------------------------------------------
    // ! Copy constructors
    // ------------------------------------------------------------------------
    PDNSolution( const PDNSolution &INPUT );

    PDNSolution &operator=( const PDNSolution &INPUT ) noexcept;

    PDNSolution( PDNSolution &&INPUT ) noexcept;

    PDNSolution &operator=( PDNSolution &&INPUT ) noexcept;

    PDNSolution( const PDNSolution * INPUT_ptr );

    virtual ~PDNSolution() noexcept;

    // ------------------------------------------------------------------------
    // ! Construct and return a random solution vector.
    // ! If input_dof_num <= 0, use pNode->get_dof().
    // ------------------------------------------------------------------------
    static PDNSolution Gen_random( const APart_Node * const &pNode,
        int input_dof_num = -1 );

    // ------------------------------------------------------------------------
    // ! Construct and return a zero solution vector.
    // ! If input_dof_num <= 0, use pNode->get_dof().
    // ------------------------------------------------------------------------
    static PDNSolution Gen_zero( const APart_Node * const &pNode,
        int input_dof_num = -1 );

    // ------------------------------------------------------------------------
    // ! Copy the INPUT's Vec, nlocal, nghost to the current vector.
    //   (similar to the copy constructor)
    // ------------------------------------------------------------------------
    void Copy(const PDNSolution &INPUT);
    
    void Copy(const PDNSolution * INPUT_ptr);

    void CopyScale(const PDNSolution &src, double a);

    void CopyScale(const PDNSolution * src_ptr, double a);
    
    // ------------------------------------------------------------------------
    // ! Update the solution entries associated with the ghost nodes
    // ------------------------------------------------------------------------
    void GhostUpdate();

    // ------------------------------------------------------------------------
    // ! Compute 1-, 2-, and infinity- Norms of the solution vector 
    // ------------------------------------------------------------------------
    double Norm_1() const;

    double Norm_2() const;
		
    double Norm_inf() const;

    // ------------------------------------------------------------------------
    // ! Perform solution = solution + a * x 
    // ------------------------------------------------------------------------
    void PlusAX(const PDNSolution &x, double a);
    
    void PlusAX(const PDNSolution * x_ptr, double a);

    // ------------------------------------------------------------------------
    // ! Perform solution = solution + a * x
    //   here x is a plain PETSc Vec object. The user is responsible for making
    //   sure that the parallel layout is compatible between x and solution. 
    // ------------------------------------------------------------------------
    void PlusAX(const Vec &x, double a);

    // ------------------------------------------------------------------------
    // ! Perform uniform scaling operation : solution = a * solution
    // ------------------------------------------------------------------------
    void ScaleValue( double a );

    // ------------------------------------------------------------------------
    // ! Perform solution = a * x + b * y
    // ------------------------------------------------------------------------
    void LinearCombination(double a, const PDNSolution &x, double b, const PDNSolution &y);
    
    // ------------------------------------------------------------------------
    // ! Get the part of the solution vector that belongs to the local 
    //   nodes and ghost nodes.
    //   If one uses a dynamic array, one should allocate it with the
    //   size nlocal + nghost, and one is also responsible for freeing
    //   the memory allocation of the local_array pointer.
    // ------------------------------------------------------------------------
    void GetLocalArray( double * const &local_array ) const;
    
    std::vector<double> GetLocalArray() const;

    // ------------------------------------------------------------------------
    // ! Assembly the vector and update its ghost values. It is just a routine
    //   calling the following things.
    //          VecAssemblyBegin(solution);
    //          VecAssemblyEnd(solution);
    //          GhostUpdate();
    //   This is called after VecSetValues to finish the assembly of vector.
    // ------------------------------------------------------------------------
    void Assembly_GhostUpdate();

    // ------------------------------------------------------------------------
    // ! Print the vec solution on screen with or without the ghost part
    // ------------------------------------------------------------------------
    void PrintWithGhost() const;
    
    void PrintNoGhost() const;

    // ------------------------------------------------------------------------
    // ! Write and Read the solution vector in PETSc binary format
    // ------------------------------------------------------------------------
    void WriteBinary(const std::string &file_name) const;
    
    void ReadBinary(const std::string &file_name);

    // ------------------------------------------------------------------------
    // ! Get the number of local and ghost nodes for the parallel vector's 
    //   local portion.
    // ------------------------------------------------------------------------
    int get_nlocalnode() const noexcept {return nlocalnode;}

    int get_nghostnode() const noexcept {return nghostnode;}

    // ------------------------------------------------------------------------
    // ! Get the number of local and ghost entries in the parallel vector's 
    //   local portion.
    // ------------------------------------------------------------------------
    int get_nlocal() const noexcept {return nlocal;}
    
    int get_nghost() const noexcept {return nghost;}

    // ------------------------------------------------------------------------
    // ! Get the number of local plus ghost entries in the parallel vector's 
    //   local portion. This can be used for allocating the array size for 
    //   GetLocalArray.
    // ------------------------------------------------------------------------
    int get_nlgn() const noexcept {return nlocal + nghost;}

    // ------------------------------------------------------------------------
    // ! Get the number of degrees of freedom of this solution vector
    // ------------------------------------------------------------------------
    int get_dof_num() const noexcept {return dof_num;}

    // ------------------------------------------------------------------------
    // ! Compare the layout of the solution vector, that is,
    //                 nlocalnode & nghostnode. 
    //   Return true if the nlocalnode and nghostnode for the two inputs
    //   are the same; return false otherwise. 
    // ------------------------------------------------------------------------
    friend bool is_layout_equal( const PDNSolution &left, 
        const PDNSolution &right ) noexcept;

  protected:
    // ------------------------------------------------------------------------
    // dof_num default value is apart_node -> get_dof, but user may reset its value.
    //     nlocalnode := apart_node -> get_nlocalnode
    //     nghostnode := apart_node -> get_nghostnode
    //     nlocal     := apart_node -> get_nlocalnode * dof_num
    //     nghost     := apart_node -> get_nghostnode * dof_num
    // ------------------------------------------------------------------------
    const int dof_num, nlocalnode, nghostnode, nlocal, nghost; 
};

#endif
