#ifndef IPLOCASSEM_HPP
#define IPLOCASSEM_HPP
// ==================================================================
// IPLocAssem.hpp
// Interface for parallel local assembly routine.
//
// This is the pure abstract class for different implementation of
// local assembly routines.
//
// The tangent matrix and the residual vector are stored in the array
// of PetscScalar type: Tangent, Residual.
//
// Author: Ju Liu
// Date: Dec. 3 2013
// ==================================================================
#include "FEAElement.hpp"
#include "AInt_Weight.hpp"

class IPLocAssem
{
  public:
    IPLocAssem()
    {
      Tangent  = nullptr;
      Residual = nullptr;
      
      sur_Tangent  = nullptr;
      sur_Residual = nullptr;
    }
    
    virtual ~IPLocAssem(){};

    PetscScalar * Tangent;
    
    PetscScalar * Residual;

    // -------------------------------------------------------------- 
    // Tangent and Residual on surface elements 
    // -------------------------------------------------------------- 
    PetscScalar * sur_Tangent;

    PetscScalar * sur_Residual;

    // -------------------------------------------------------------- 
    // ! Get degree of freedom of this problem. In segregated algorithms
    //   this dof returns the fully coupled multiphysics problem's dof.
    // -------------------------------------------------------------- 
    virtual int get_dof() const = 0;

    // -------------------------------------------------------------- 
    // ! Get degree of freedom of the matrix. In segregated algorithms, 
    //   this dof returns the actually implicit solver's dof per node.
    //   In fully coupled fashions, this defaults to the get_dof function.
    // -------------------------------------------------------------- 
    virtual int get_dof_mat() const {return get_dof();}

    // -------------------------------------------------------------- 
    // ! Get the number of ebc functions implemented inside this 
    //   local assembly routine
    // -------------------------------------------------------------- 
    virtual int get_num_ebc_fun() const
    {SYS_T::commPrint("Warning: IPLocAssem::get_num_ebc_fun is not implemented. \n");
      return 0;}

    // -------------------------------------------------------------- 
    // ! Assign all values in Tangent matrix 0.0
    //   Call this function before assembly to zero everything in container
    // -------------------------------------------------------------- 
    virtual void Zero_Tangent_Residual() = 0;

    virtual void Zero_sur_Tangent_Residual()
    {
      SYS_T::print_fatal("Error: Zero_sur_Tangent_Residual is not implemented.\n");
    }

    // -------------------------------------------------------------- 
    // ! Assign all values in Residual vector 0.0
    //   Call this function before assembly to zero everything in container
    // -------------------------------------------------------------- 
    virtual void Zero_Residual() = 0;

    virtual void Zero_sur_Residual()
    {
      SYS_T::print_fatal("Error: Zero_sur_Residual is not implemented. \n");
    }

    // -------------------------------------------------------------- 
    // ! Give nonzero pattern of the sparse matrix 
    // -------------------------------------------------------------- 
    virtual void Assem_Estimate() = 0;

    // -------------------------------------------------------------- 
    // ! Assembly element residual vector: Residual
    // \para vec_a: input vector a -- displacement / current solution
    // \para vec_b: input vector b -- velocity / next solution
    // \para element: the element quadrature info
    // \para eleCtrlPts: this element's control points
    // \para wight: the corresponding quadrature weights    
    // -------------------------------------------------------------- 
    virtual void Assem_Residual(
        double time, double dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}
    
    
    // -------------------------------------------------------------- 
    // ! Assembly element residual vector without precached quadrature info for
    //   3D element.
    //   Element quadrature info is computed inside the element assembly
    //   routine.
    // \para vec_a: input vector a -- displacement / current solution
    // \para vec_b: input vector b -- velocity / next solution
    // \para bs : Bernstein basis function precomputed in s direction
    // \para bt : Bernstein basis function precomputed in t direction
    // \para bu : Bernstein basis function precomputed in u direction
    // \para extractor : Bezier extraction operator
    // -------------------------------------------------------------- 
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}
  

    // \para element: the container for element basis functions
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}


    // \para element: the container for element basis functions
    // input contains the two groups of local solution vectors and their time
    // derivatives: vec_a (disp)
    //              vec_da (dot disp)
    //              vec_b (pres-velo)
    //              vec_db (dot_pres-velo)
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_da,
        const double * const &vec_b,
        const double * const &vec_db,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}


    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}

    
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}


    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    // \para vec_c : designed for prestress
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &vec_c,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual(...) is not implemented. \n");}


    // ! Assembly element residual vector and tangent matrix
    // \para vec_a: input vector a -- displacement / current solution
    // \para vec_b: input vector b -- velocity / next solution
    // \para element: the element quadrature info
    // \para eleCtrlPts: this element's control points
    // \para wight: the corresponding quadrature weights    
    virtual void Assem_Tangent_Residual(
        double time, double dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual(...) is not implemented. \n");}


    // ! Assembly element residual vector and tengent matrix without cached
    //   quadrature info for 3D element.
    // \para vec_a: input vector a -- displacement / current solution
    // \para vec_b: input vector b -- velocity / next solution
    // \para bs : Bernstein basis function precomputed in s direction
    // \para bt : Bernstein basis function precomputed in t direction
    // \para bu : Bernstein basis function precomputed in u direction
    // \para extractor : Bezier extraction operator
    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual(...) is not implemented. \n");}



    // \para element: the container for element basis functions
    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual(...) is not implemented. \n");}


    // input contains the two groups of local solution vectors and their time
    // derivatives: vec_a (disp)
    //              vec_da (dot disp)
    //              vec_b (pres-velo)
    //              vec_db (dot_pres-velo)
    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_da,
        const double * const &vec_b,
        const double * const &vec_db,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual(...) is not implemented. \n");}


    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual(...) is not implemented. \n");}


    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    // vec_c : designed for prestress
    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const double * const &vec_c,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual(...) is not implemented. \n");}


    // -------------------------------------------------------------------
    // ! Assembly the element residual vector resulting from the boundary
    //   integral on the top boundary of the 3D element.
    //   The boundary quadrature points will be calculated inside this local
    //   assembly routine. 
    //   \para time : the current time
    //   \para dt   : the time step
    //   \para vec_a : displacement / current solution 
    //   \para vec_b : velocity / next solution
    //   \para eleCtrlPts_ : element control points
    //   \para ext_ : extraction operator
    virtual void Assem_Residual_TopFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: Assem_Residual_TopFace() is not implemented. \n");}


    // ! Assembly the element residual vector resulting from the boundary
    //   integral on the top boundary of the 3D element.
    //   The boundary quadrature points will be calculated inside this local
    //   assembly routine. 
    //   \para time : the current time
    //   \para dt   : the time step
    //   \para vec_a : displacement / current solution 
    //   \para vec_b : velocity / next solution
    //   \para element : element holder
    //   \para eleCtrlPts_ : element control points
    //   \para ext_ : extraction operator
    virtual void Assem_Residual_TopFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: Assem_Residual_TopFace() is not implemented. \n");}


    // -------------------------------------------------------------------
    // ! Assembly the element residual vector resulting from the boundary
    //   integral on the bottom boundary of the 3D element.
    //   The boundary quadrature points will be calculated inside this local
    //   assembly routine. 
    //   \para time : the current time
    //   \para dt   : the time step
    //   \para vec_a : displacement / current solution 
    //   \para vec_b : velocity / next solution
    //   \para eleCtrlPts_ : element control points
    //   \para ext_ : extraction operator
    virtual void Assem_Residual_BotFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: Assem_Residual_BotFace() is not implemented. \n");}


    // ! Assembly the element residual vector resulting from the boundary
    //   integral on the bottom boundary of the 3D element.
    //   The boundary quadrature points will be calculated inside this local
    //   assembly routine. 
    //   \para time : the current time
    //   \para dt   : the time step
    //   \para vec_a : displacement / current solution 
    //   \para vec_b : velocity / next solution
    //   \para element : element holder
    //   \para eleCtrlPts_ : element control points
    //   \para ext_ : extraction operator
    virtual void Assem_Residual_BotFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: Assem_Residual_BotFace() is not implemented. \n");}


    // ------------------------------------------------------------------- 
    // ! Assembly the element residual vector resulting from the boundary
    //   integral on the left boundary of the 3D element.
    //   The boundary quadrature points will be calculated inside this local
    //   assembly routine. 
    //   \para time : the current time
    //   \para dt   : the time step
    //   \para vec_a : displacement / current solution 
    //   \para vec_b : velocity / next solution
    //   \para element : element holder
    //   \para eleCtrlPts_ : element control points
    //   \para ext_ : extraction operator
    virtual void Assem_Residual_LefFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: Assem_Residual_LefFace() is not implemented. \n");}



    // ------------------------------------------------------------------- 
    // ! Assembly the element residual vector resulting from the boundary
    //   integral on the right boundary of the 3D element.
    //   The boundary quadrature points will be calculated inside this local
    //   assembly routine. 
    //   \para time : the current time
    //   \para dt   : the time step
    //   \para vec_a : displacement / current solution 
    //   \para vec_b : velocity / next solution
    //   \para element : element holder
    //   \para eleCtrlPts_ : element control points
    //   \para ext_ : extraction operator
    virtual void Assem_Residual_RigFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: Assem_Residual_RigFace() is not implemented. \n");}


    // ------------------------------------------------------------------- 
    // ! Assembly the element residual vector resulting from the boundary
    //   integral on the front boundary of the 3D element.
    //   The boundary quadrature points will be calculated inside this local
    //   assembly routine. 
    //   \para time : the current time
    //   \para dt   : the time step
    //   \para vec_a : displacement / current solution 
    //   \para vec_b : velocity / next solution
    //   \para element : element holder
    //   \para eleCtrlPts_ : element control points
    //   \para ext_ : extraction operator
    virtual void Assem_Residual_FroFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: Assem_Residual_FroFace() is not implemented. \n");}


    // ------------------------------------------------------------------- 
    // ! Assembly the element residual vector resulting from the boundary
    //   integral on the back boundary of the 3D element.
    //   The boundary quadrature points will be calculated inside this local
    //   assembly routine. 
    //   \para time : the current time
    //   \para dt   : the time step
    //   \para vec_a : displacement / current solution 
    //   \para vec_b : velocity / next solution
    //   \para element : element holder
    //   \para eleCtrlPts_ : element control points
    //   \para ext_ : extraction operator
    virtual void Assem_Residual_BacFace(
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: Assem_Residual_BacFace() is not implemented. \n");}


    // -------------------------------------------------------------------
    // ! Assembly the mass matrix
    // \para element: the element quadrature info
    // \para wight: the corresponding quadrature weights    
    virtual void Assem_Mass(
        const FEAElement * const &element,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Mass(...) is not implemented. \n");}


    // ! Assembly the mass matrix without cached quadrature info for 3D element
    // \para bs : Bernstein basis function precomputed in s direction
    // \para bt : Bernstein basis function precomputed in t direction
    // \para bu : Bernstein basis function precomputed in u direction
    // \para extractor : Bezier extraction operator
    virtual void Assem_Mass(
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Mass(...) is not implemented. \n");}


    // ! Assembly the mass matrix with residual vector
    //   this function may be used for initialization consistent velocity 
    //   in generalized-alpha method
    // \para vec_a: input vector a -- displacement / current solution
    // \para element: the element quadrature info
    // \para eleCtrlPts: this element's control points
    // \para wight: the corresponding quadrature weights    
    virtual void Assem_Mass_Residual(
        const double * const &vec_a,
        const FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Mass_Residual(...) is not implemented. \n");}


    // ! Assembly the mass matrix and residual vector without cached quadrature
    //   info for 3D element.
    // \para mSize : mesh size info
    // \para bs : Bernstein basis function precomputed in s direction
    // \para bt : Bernstein basis function precomputed in t direction
    // \para bu : Bernstein basis function precomputed in u direction
    // \para extractor : Bezier extraction operator
    virtual void Assem_Mass_Residual(
        const double * const &vec_a,
        const int &eindex,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Mass_Residual(...) is not implemented. \n");}



    // \para element: the container for element basis functions
    virtual void Assem_Mass_Residual(
        const double * const &curr,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Mass_Residual(...) is not implemented. \n");}


    // input contains the two groups of local solution vectors and their time
    // derivatives: vec_a (disp)
    //              vec_da (dot disp)
    //              vec_b (pres-velo)
    //              vec_db (dot_pres-velo)
    virtual void Assem_Mass_Residual(
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &eleCtrlPts_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z,
        const AInt_Weight * const &weight )
    {SYS_T::commPrint("Warning: this Assem_Mass_Residual(...) is not implemented. \n");}

    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    virtual void Assem_Mass_Residual(
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Mass_Residual(...) is not implemented. \n");}


    // \para element: the container for classical element routine. It only
    //                requires the x-y-z coordinates for the nodes and the 
    //                volumetric quadrature routine to generate the basis
    //                functions
    // \vec_c : designed for prestress
    virtual void Assem_Mass_Residual(
        const double * const &vec_b,
        const double * const &vec_c,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Mass_Residual(...) is not implemented. \n");}


    // Perform Elemental BC surface integration for elemental BC id ebc_id.
    // Based on ebc_id, the traction forcing function will be called accordingly
    // inside the local assembly routine.
    // \para element: the container for the element, only requires the geometry
    //                information for the control points, and the quadrature
    //                info to generate basis function info.
    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC is not implemented.\n");}

    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC is not implemented.\n");}


    // Perform elemental BC surface integration for backflow stabilization
    // for the residual only.
    virtual void Assem_Residual_BackFlowStab(
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_BackFlowStab is not implemented.\n");}


    // Perform elemental BC surface integration for backflow stabilization
    // for the residual as well as the tangent matrix
    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual_BackFlowStab is not implemented.\n");}


    // Perform Elemental BC surface integration for elemental BC id ebc_id.
    // Based on ebc_id, the traction forcing function will be called accordingly
    // inside the local assembly routine.
    // \para in_x, in_y, in_z : interior point to the element
    // \para element: the container for the element, only requires the geometry
    //                information for the control points, and the quadrature
    //                info to generate basis function info.
    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        const double &in_x, const double &in_y, const double &in_z,
        const double * const &vec_a,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC is not implemented.\n");}


    // Perform Elemental BC surface integration for elemental BC id ebc_id and
    // for resistance type BC.
    // Based on ebc_id, the traction forcing function will be called accordingly
    // inside the local assembly routine.
    // \para element: the container for the element, only requires the geometry
    //                information for the control points, and the quadrature
    //                info to generate basis function info.
    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id,
        const double &flow_rate,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC_Resistance is not implemented.\n");}

    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id,
        const double &flow_rate,
        const double * const &vec_b,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC_Resistance is not implemented.\n");}


    // Perform elemental BC surface integration for the coupled momentum FSI method, in which
    // the fluid is coupled with a thin-walled membrane formulation for the vascular wall. 
    // \para dot_sol:       dot pressure, dot velocity
    // \para sol_wall_disp: wall displacement
    // \para element:       container for membrane element. Only requires nodal x,y,z-coordinates
    //                      and the surface quadrature rule to generate basis functions.
    // \para ele_thickness: wall thickness
    // \para ele_youngsmod: wall youngsmod
    // \para qua_prestress: prestress tensor at each quadrature point
    virtual void Assem_Residual_EBC_Wall(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &sol_wall_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &ele_thickness,
        const double * const &ele_youngsmod,
        const double * const &ele_springconst,
        const double * const &ele_dampingconst,
        const double * const &qua_prestress,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC_Wall is not implemented.\n");}

    virtual void Assem_Tangent_Residual_EBC_Wall(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        const double * const &sol_wall_disp,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const double * const &ele_thickness,
        const double * const &ele_youngsmod,
        const double * const &ele_springconst,
        const double * const &ele_dampingconst,
        const double * const &qua_prestress,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual_EBC_Wall is not implemented.\n");}

    
    // ! Get the model parameter 1
    //   This function is used to pass out the parameters appearing in the weak
    //   form, such as the Reynolds number, Capallarity number, etc.
    //   The definition of this function varies depending on the derived class
    virtual double get_model_para_1() const
    {
      SYS_T::commPrint("Warning: get_model_para_1() is not implemented. \n");
      return 0.0;
    }

    virtual double get_model_para_2() const
    {
      SYS_T::commPrint("Warning: get_model_para_2() is not implemented. \n");
      return 0.0;
    }

    // This is a function in local assembly that calculates the flow rate.
    virtual double get_flowrate( const double * const &vec, 
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {
      SYS_T::commPrint("Warning: get_flowrate() is not implemented. \n");
      return 0.0;
    }

    // This is a function in local assembly that calculates the pressure
    // integrated over surface: int_{Gamma} p dA
    // as well as the area of the surface: int_{Gamma} 1 dA
    virtual void get_pressure_area( const double * const &vec, 
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        double &pres, double &area )
    {
      SYS_T::commPrint("Warning: get_pressure_area() is not implemented. \n");
    }

    // Computes the Cauchy stress in the wall for the coupled momentum FSI method 
    virtual void get_Wall_CauchyStress(
        const double * const &sol_wall_disp,
        const FEAElement * const &element,
        const double * const &ele_youngsmod,
        std::vector<Matrix_3x3> &stress )
    {
      SYS_T::commPrint("Warning: get_Wall_CauchyStress() is not implemented. \n");
    }

};

#endif
