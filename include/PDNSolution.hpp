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
#include "Math_Tools.hpp"
#include "APart_Node.hpp"
#include "FEANode.hpp"
#include "ALocal_InflowBC.hpp"

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
    PDNSolution( const APart_Node * const &pNode, const int &input_dof_num );
    
    // ------------------------------------------------------------------------
    // ! Copy constructors
    // ------------------------------------------------------------------------
    PDNSolution( const PDNSolution &INPUT );

    PDNSolution( const PDNSolution * const &INPUT_ptr );

    virtual ~PDNSolution();

    // ------------------------------------------------------------------------
    // ! Generate a random solution vector.
    // ------------------------------------------------------------------------
    virtual void Gen_random();

    // ------------------------------------------------------------------------
    // ! Copy the INPUT's Vec, nlocal, nghost to the current vector.
    //   (similar to the copy constructor)
    // ------------------------------------------------------------------------
    virtual void Copy(const PDNSolution &INPUT);
    
    virtual void Copy(const PDNSolution * const &INPUT_ptr);
    
    // ------------------------------------------------------------------------
    // ! Update the solution entries associated with the ghost nodes
    // ------------------------------------------------------------------------
    virtual void GhostUpdate();

    // ------------------------------------------------------------------------
    // ! Compute 1-, 2-, and infinity- Norms of the solution vector 
    // ------------------------------------------------------------------------
		virtual double Norm_1() const;
		
    virtual double Norm_2() const;
		
    virtual double Norm_inf() const;

    // ------------------------------------------------------------------------
    // ! Perform solution = solution + a * x 
    // ------------------------------------------------------------------------
    virtual void PlusAX(const PDNSolution &x, const double &a);
    
    virtual void PlusAX(const PDNSolution * const &x_ptr, const double &a);

    // ------------------------------------------------------------------------
    // ! Perform solution = solution + a * x
    //   here x is a plain PETSc Vec object. The user is responsible for making
    //   sure that the parallel layout is compatible between x and solution. 
    // ------------------------------------------------------------------------
    virtual void PlusAX(const Vec &x, const double &a);

    // ------------------------------------------------------------------------
    // ! Perform uniform scaling operation : solution = a * solution
    // ------------------------------------------------------------------------
    virtual void ScaleValue( const double &a );
    
    // ------------------------------------------------------------------------
    // ! Get the part of the solution vector that belongs to the local 
    //   nodes and ghost nodes.
    //   If one uses a dynamic array, one should allocate it with the
    //   size nlocal + nghost, and one is also responsible for freeing
    //   the memory allocation of the local_array pointer.
    // ------------------------------------------------------------------------
    virtual void GetLocalArray( double * const &local_array ) const;
    
    virtual std::vector<double> GetLocalArray() const;

    // ------------------------------------------------------------------------
    // ! Assembly the vector and update its ghost values. It is just a routine 
    //   calling the following things. 
    //          VecAssemblyBegin(solution);
    //          VecAssenblyEnd(solution);
    //          GhostUpdate();
    //   This is called after VecSetValues to finish the assembly of vector.
    // ------------------------------------------------------------------------
    virtual void Assembly_GhostUpdate();

    // ------------------------------------------------------------------------
    // ! Print the vec solution on screen with or without the ghost part
    // ------------------------------------------------------------------------
    virtual void PrintWithGhost() const;
    
    virtual void PrintNoGhost() const;

    // ------------------------------------------------------------------------
    // ! Write and Read the solution vector in PETSc binary format
    // ------------------------------------------------------------------------
    virtual void WriteBinary(const std::string &file_name) const;
    
    virtual void ReadBinary(const std::string &file_name) const;

    // ------------------------------------------------------------------------
    // ! Get the number of local and ghost nodes for the parallel vector's 
    //   local portion.
    // ------------------------------------------------------------------------
    virtual int get_nlocalnode() const {return nlocalnode;}

    virtual int get_nghostnode() const {return nghostnode;}

    // ------------------------------------------------------------------------
    // ! Get the number of local and ghost entries in the parallel vector's 
    //   local portion.
    // ------------------------------------------------------------------------
    virtual int get_nlocal() const {return nlocal;}
    
    virtual int get_nghost() const {return nghost;}

    // ------------------------------------------------------------------------
    // ! Get the number of local plus ghost entries in the parallel vector's 
    //   local portion. This can be used for allocating the array size for 
    //   GetLocalArray.
    // ------------------------------------------------------------------------
    virtual int get_nlgn() const {return nlocal + nghost;}

    // ------------------------------------------------------------------------
    // ! Get the number of degrees of freedom of this solution vector
    // ------------------------------------------------------------------------
    virtual int get_dof_num() const {return dof_num;}

    // ------------------------------------------------------------------------
    // ! Compare the layout of the solution vector, that is,
    //                 nlocalnode & nghostnode. 
    //   Return true if the nlocalnode and nghostnode for the two inputs
    //   are the same; return false otherwise. 
    // ------------------------------------------------------------------------
    friend bool is_layout_equal( const PDNSolution &left, const PDNSolution &right );

    virtual void randomly_perturbed_parabolic_inflow( const double &std_dev,
        const APart_Node * const &pNode,
        const FEANode * const &fNode,
        const ALocal_InflowBC * const &infbc )
    {
      SYS_T::print_fatal("Error: PDNSolution::randomly_perturbed is not implemented.\n");
      return;
    } 

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
