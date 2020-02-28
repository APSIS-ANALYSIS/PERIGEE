#ifndef PDNSOLUTION_HPP
#define PDNSOLUTION_HPP
// ==================================================================
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
//                node.
//
// Author: Ju Liu
// Date: Nov. 23th 2013
// ==================================================================
#include "APart_Node.hpp"

class PDNSolution
{
  public:
    Vec solution;  
    
    // --------------------------------------------------------------
    // Construct a solution vec compatible with the analysis node 
    // partition
    // --------------------------------------------------------------
    PDNSolution( const APart_Node * const &pNode );
    
    // --------------------------------------------------------------
    // Construct a solution vec compatible with the analysis node 
    // partition, but with a different dof_num from pNode->dof. 
    // The users specify a dof number for the solution class.
    // --------------------------------------------------------------
    PDNSolution( const APart_Node * const &pNode,
       const int &input_dof_num );
    
    // --------------------------------------------------------------
    // Copy constructor
    // --------------------------------------------------------------
    PDNSolution( const PDNSolution &INPUT );

    PDNSolution( const PDNSolution * const &INPUT_ptr );

    virtual ~PDNSolution();

    // --------------------------------------------------------------
    // ! Generate a random solution vector.
    //   This is used mainly for debugging.
    // --------------------------------------------------------------
    virtual void Gen_random();

    // --------------------------------------------------------------
    // ! Copy the INPUT's Vec, nlocal, nghost to the current vector.
    //   (similar to the copy constructor)
    // --------------------------------------------------------------
    virtual void Copy(const PDNSolution &INPUT);
    
    virtual void Copy(const PDNSolution * const &INPUT_ptr);
    
    // --------------------------------------------------------------
    // ! Update the ghost part
    //   This function is often used in initialization
    // --------------------------------------------------------------
    virtual void GhostUpdate();

    // --------------------------------------------------------------
    // ! Compute Norms of solution vector 
    // --------------------------------------------------------------
		virtual double Norm_1() const;
		virtual double Norm_2() const;
		virtual double Norm_inf() const;

    // --------------------------------------------------------------
    // ! Perform solution = solution + a * x 
    // --------------------------------------------------------------
    virtual void PlusAX(const PDNSolution &x, const double &a);
    
    virtual void PlusAX(const PDNSolution * const &x_ptr, const double &a);

    // --------------------------------------------------------------
    // ! Perform += a_i * x_i.
    //   x has the same length with the solution, and length of x is
    //   a multiple of (na + nb);
    //   a_i = a for (na+nb)*k    <= i < (na + nb)*k + na
    //   a_i = b for (na+nb)*k+na <= i < (na+nb) * (k+1)
    //   This function call is useful for the special Generalized-alpha
    //   time integration used for VMS-NS solver. 
    //   Refer to CMAME 197(2007), pp 182 for details
    // --------------------------------------------------------------
    virtual void PlusAiX(PDNSolution &x, const double &a,
        const double &b, const int &na, const int &nb );

    // -------------------------------------------------------------
    // ! Perform solution[ii*dofNum + jj] += a[jj] * x[ii*dofNum + jj].
    //   This is a generalized version for the PlusAiX( PDNSolution,
    //   const double, const double, const int, const int ).
    //   The users should make sure that aa.size() == dofNum of the
    //   solution vector. The PlusAiX(x,a,b, na,nb) is equivalent to
    //   setting aa = [a,..., a, b, ..., b].
    //                 na times   nb times
    // -------------------------------------------------------------
    virtual void PlusAiX( const PDNSolution &xx, 
        const std::vector<double> &aa );

    // --------------------------------------------------------------
    // ! Perform uniform scaling operation : solution = a * solution
    // --------------------------------------------------------------
    virtual void ScaleValue(const double &a);

    // --------------------------------------------------------------
    // ! Get the part of the solution vector that belongs to the local 
    //   nodes and ghost nodes. Here the local_array should have been 
    //   allocated with length equal to nlocal + nghost.
    //   The user is responsible for freeing the memory allocation 
    //   after the task is done.
    //   Note: to get the local array, one does NOT need the APart_Node
    //   class as an input. The GetLocalArray( double * const & ) is
    //   one that corrects this issue and use safer pointer usage. A
    //   safer way is to call the function with std::vector as output.
    // --------------------------------------------------------------
    virtual void GetLocalArray( double * &local_array,
        const APart_Node * const &pNode ) const;

    virtual void GetLocalArray( double * const &local_array ) const;
    
    virtual void GetLocalArray( std::vector<double> &local_array ) const;

    // --------------------------------------------------------------
    // ! Assembly the vector and update the ghost node values. It is
    //   just a routine calling the following things. This is called
    //   after VecSetValues calling to finish the assembly of vector.
    //          VecAssemblyBegin(solution);
    //          VecAssenblyEnd(solution);
    //          GhostUpdate();
    // --------------------------------------------------------------
    virtual void Assembly_GhostUpdate();

    // --------------------------------------------------------------
    // ! Print the vec solution on screen with or without ghost part
    // --------------------------------------------------------------
    virtual void PrintWithGhost() const;
    virtual void PrintNoGhost() const;

    // --------------------------------------------------------------
    // ! Write and Read the solution vector in PETSc binary format
    // --------------------------------------------------------------
    virtual void WriteBinary(const char * const &file_name) const;
    virtual void ReadBinary(const char * const &file_name) const;

    // --------------------------------------------------------------
    // ! Get the number of local and ghost entries in the parallel 
    //   vector's local portion.
    // --------------------------------------------------------------
    virtual int get_nlocal() const {return nlocal;}
    virtual int get_nghost() const {return nghost;}
    
    // --------------------------------------------------------------
    // ! Get the number of degrees of freedom 
    // --------------------------------------------------------------
    virtual int get_dof_num() const {return dof_num;}

  private:
   int nlocal, nghost, dof_num; 
};

#endif
