#ifndef FEAELEMENT_HPP
#define FEAELEMENT_HPP
// ==================================================================
// FEAElement.hpp
// Object: This is the interface for Finite Element classes.
// It records:
// 1. Element basic indices;
// 2. Element nodes through the IEN array;
// 3. reference element position / hx hy hz; 
// 4. extraction operators.
//
// It should be capable of representing H1, H-div, H-curl, L2 conforming
// elements.
//
// Author: Ju Liu
// Date created: Nov. 6 2013
// ==================================================================
#include <cassert>
#include "ALocal_IEN.hpp"
#include "IALocal_IEN.hpp"
#include "IALocal_meshSize.hpp"
#include "FEANode.hpp"
#include "BernsteinBasis_Array.hpp"
#include "IBernsteinBasis.hpp"
#include "Matrix_3x3.hpp"

class FEAElement
{
  public:
    // Constructor
    FEAElement(){};
    
    // Destructor
    virtual ~FEAElement(){};

    // Return this element's local ordering index
    virtual int get_elemIndex() const
    {SYS_T::commPrint("Warning: get_elemIndex is not implemented. \n"); return 0;}

    // Return this element's dimension
    virtual int get_elemDim() const = 0;

    // Return this element's Type, which defines the type of different 
    // elements defined on this single element domain.
    // H1 2D NURBS der0:      Type = 120
    // H1 2D NURBS der1:      Type = 121
    // H1 2D NURBS     :      Type = 122
    // H1 3D NURBS der0:      Type = 130
    // H1 3D NURBS der1:      Type = 131
    // H1 3D NURBS der2:      Type = 132
    // H1 3D NURBS der1_wJac: Type = 133
    // H1 3D NURBS der2_vms:  Type = 134
    // H-div B-spline:        Type = 220
    // H1 2D T-splines der 0: Type = 520
    // H1 2D T-splines der 1: Type = 521
    // H1 2D T-splines der 2: Type = 522
    // H1 2D T-splines der 1 lap: Type = 523
    // H1 2D B-splines der 0: Type = 620
    // H1 2D B-splines der 1: Type = 621
    // H1 2D B-splines der 2: Type = 622
    // H1 3D B-splines der 2: type = 632
    virtual int get_Type() const
    {SYS_T::commPrint("Warning: get_Type is not implemented. \n"); return -1;}

    // Return the number of different basis functions defined on this element.
    // H1 NURBS: 1
    // H-div 2D: 2
    // H-div 3D: 3
    virtual int get_numType() const = 0;

    // Return the number of nodes for each type elements
    virtual int get_nLocBas() const = 0;
    
    // Return the quadrature info
    virtual int get_numQuapts() const = 0;

    // Return true if element size is non-zero
    virtual bool is_sizeNonzero() const
    {
      SYS_T::commPrint("Warning: is_sizeNonzero is not implemented. \n"); 
      return false;
    }
  
    // Clear memory of basis function quadrature
    virtual void clearBasisCache()
    {SYS_T::commPrint("Warning: clearBasisCache is not implemented. \n");}
    
    // --------------------------------------------------------------
    // print function
    // -- note: in the old implementation, I use print(). In newer
    //          implementations, I use print_info().
    // --------------------------------------------------------------
    virtual void print() const 
    {SYS_T::commPrint("Warning: print is not implemented. \n");}
    
    virtual void print_info() const 
    {SYS_T::commPrint("Warning: print is not implemented. \n");}

    // Return the memory usage of this class in bytes
    virtual double get_memory_usage() const = 0; 

    // --------------------------------------------------------------    
    // Calculate the element size
    // --------------------------------------------------------------    
    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z ) const
    {SYS_T::commPrint("Warning: get_h is not implemented. \n"); return 0.0;}

    virtual double get_h( const double * const &ctrl_x,
        const double * const &ctrl_y ) const
    {SYS_T::commPrint("Warning: get_h is not implemented. \n"); return 0.0;}


    // --------------------------------------------------------------    
    // Build Basis function quadrature info
    // --------------------------------------------------------------    
    // Build basis function values at quadrature points for 3D elements
    // \para bs/t/u : Bernstein polynomial value at quadrature points
    // \para ctrl_x/y/z/w : control points and weights
    // \para ext_x/y/z : Bezier extraction operator in x/y/z direction
    virtual void buildBasis(
        const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z,
        const double * const &ctrl_w,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: this buildBasis is not implemented. \n");}


    // Build basis function values at quadrature points for 2D elements
    // \para hx, hy : parametric element size
    // \para bs, bt : Bezier elements
    // \para ctrl_x/y/w : control points' x-y-w value
    // \para ext_x/y : extraction operator
    virtual void buildBasis(
        const double &hx, const double &hy,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_w,
        const double * const &ext_x,
        const double * const &ext_y )
    {SYS_T::commPrint("Warning: this buildBasis is not implemented. \n");}

    // Build basis funtion values at quadrature points for 2D elements.
    // This is for irrational cases where weights can be ignored. 
    virtual void buildBasis(
        const double &hx, const double &hy,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ext_x,
        const double * const &ext_y )
    {SYS_T::commPrint("Warning: this buildBasis is not implemented. \n");}

    // Build 3D basis -- 3D B-spine case
    virtual void buildBasis( const double &hx, const double &hy, const double &hz,
        const BernsteinBasis_Array * const &bs,
        const BernsteinBasis_Array * const &bt,
        const BernsteinBasis_Array * const &bu,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z,
        const double * const &ext_x,
        const double * const &ext_y,
        const double * const &ext_z )
    {SYS_T::commPrint("Warning: buildBasis() is not implemented. \n");}


    // Build 3D basis -- FEM
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x,
        const double * const &ctrl_y,
        const double * const &ctrl_z )
    {SYS_T::commPrint("Warning: buildBasis() is not implemented. \n");}


    // Build 2D basis -- FEM
    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x, const double * const &ctrl_y )
    {SYS_T::commPrint("Warning: buildBasis() is not implemented. \n");}


    // ------------------------------------------------------------------------
    // Get functions : Obtain the value of basis functions and their derivatives.
    // ------------------------------------------------------------------------
    // Return basis function value using a dynamic double array.
    // Users are responsible for allocating proper memory for basis,
    // and delete the pointer after use.
    virtual void get_R( const int &quaindex, double * const &basis ) const
    {SYS_T::commPrint("Warning: get_R is not implemented. \n");} 

    virtual std::vector<double> get_R( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_R is not implemented. \n");
      return {};
    } 

    // ------------------------------------------------------------------------    
    // Return basis function and all 1st order derivatives value 
    // using a dynamic double array.
    // Users are responsible for allocating proper memory for basis,
    // and delete the pointer after use.
    // ------------------------------------------------------------------------    
    // 3D case
    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y, double * const &basis_z )
      const {SYS_T::commPrint("Warning: get_gradR is not implemented. \n");} 

    // 2D case
    virtual void get_gradR( const int &quaindex, double * const &basis_x,
        double * const &basis_y )
      const {SYS_T::commPrint("Warning: get_gradR is not implemented. \n");} 

    // 2D case:
    virtual void get_R_gradR( const int &quaindex, double * const &basis, 
        double * const &basis_x, double * const &basis_y )
      const {SYS_T::commPrint("Warning: get_R_gradR is not implemented. \n");} 

    // 3D case:
    virtual void get_R_gradR( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y, double * const &basis_z ) 
      const {SYS_T::commPrint("Warning: get_R_gradR is not implemented. \n");}

    virtual std::vector<double> get_dR_dx( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_dR_dx is not implemented. \n");
      return {};
    }

    virtual std::vector<double> get_dR_dy( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_dR_dy is not implemented. \n");
      return {};
    }

    virtual std::vector<double> get_dR_dz( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_dR_dz is not implemented. \n");
      return {};
    }

    // ------------------------------------------------------------------------    
    // R, gradR, and Laplacian R
    // ------------------------------------------------------------------------    
    virtual void get_3D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, double * const &basis_y,
        double * const &basis_z, double * const &basis_xx,
        double * const &basis_yy, double * const &basis_zz )
      const {SYS_T::commPrint("Warning: get_3DLaplacianR is not implemented. \n");}

    virtual void get_2D_R_gradR_LaplacianR( const int &quaindex,
        double * const &basis, double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy )
      const {SYS_T::commPrint("Warning: get_2DLaplacianR is not implemented. \n");}

    // ------------------------------------------------------------------------    
    // R, gradR, and grad gradR
    // ------------------------------------------------------------------------    
    virtual void get_2D_R_dR_d2R( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y,
        double * const &basis_xx, double * const &basis_yy, double * const &basis_xy
        ) const
    {SYS_T::commPrint("Warning: get_2D_R_dR_d2R is not implemented. \n");}

    virtual void get_3D_R_dR_d2R( const int &quaindex, double * const &basis,
        double * const &basis_x, double * const &basis_y, double * const &basis_z,
        double * const &basis_xx, double * const &basis_yy, double * const &basis_zz,
        double * const &basis_xy, double * const &basis_xz, double * const &basis_yz ) 
      const {SYS_T::commPrint("Warning: get_3D_R_dR_d2R is not implemented. \n");}

    virtual std::vector<double> get_d2R_dxx( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_d2R_dxx is not implemented. \n");
      return {};
    }

    virtual std::vector<double> get_d2R_dyy( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_d2R_dyy is not implemented. \n");
      return {};
    }

    virtual std::vector<double> get_d2R_dzz( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_d2R_dzz is not implemented. \n");
      return {};
    }

    virtual std::vector<double> get_d2R_dxy( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_d2R_dxy is not implemented. \n");
      return {};
    }

    virtual std::vector<double> get_d2R_dxz( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_d2R_dxz is not implemented. \n");
      return {};
    }

    virtual std::vector<double> get_d2R_dyz( const int &quaindex ) const
    {
      SYS_T::commPrint("Warning: get_d2R_dyz is not implemented. \n");
      return {};
    }

    // ------------------------------------------------------------------------    
    // Return the Jacobian determinant
    // ------------------------------------------------------------------------    
    virtual double get_detJac(const int &quaindex) const = 0;

    // ------------------------------------------------------------------------    
    // Return the Jacobian matrix in rows, i.e. dx_dxi, dx_deta, dx_dzeta,
    // dy_dxi, ... dz_deta, dz_dzeta. The size of teh return vector is 9 for 3d,
    // and 4 for 2d, for exampel. 
    // ------------------------------------------------------------------------    
    virtual void get_Jacobian(const int &quaindex, double * const &jac_value) 
      const {SYS_T::commPrint("Warning: get_Jacobian is not implemented. \n");} 

    // ------------------------------------------------------------------------    
    // Return the inversion of the Jacobian matrix, the 9(4) components are
    // given. The output array dxi_dx is a 9 or 4 component array for 3D / 2D
    // case.
    // ------------------------------------------------------------------------    
    virtual void get_invJacobian(const int &quaindex, double * const &dxi_dx)
      const {SYS_T::commPrint("Warning: get_invJacobian is not implemented. \n");} 

    // ------------------------------------------------------------------------    
    // Return the global-to-lamina rotation matrix for transforming between the
    // global and local lamina coordinates. This routine is used in membrane
    // and shell elements.
    // Reference: TJRH Linear finite element book page 386.
    // ------------------------------------------------------------------------    
    virtual Matrix_3x3 get_rotationMatrix( const int &quaindex ) const 
    {
      SYS_T::commPrint("Warning: get_rotationMatrix is not implemented. \n"); 
      return Matrix_3x3();
    }

    // ------------------------------------------------------------------------    
    // Unit outward normal vector:
    // The following virtual function returns the unit normal vector on the
    // face of a 2D/3D element. 
    // In a 3D hexahedron-shaped element, there are six faces. In its
    //    reference coordinate (s,t,u), the six faces are
    //      bottom   (s,t,0)     top      (s,t,1)
    //      left     (s,0,u)     right    (s,1,u)
    //      front    (1,t,u)     back     (0,t,u)
    // The unit outward normal vector on each face is calculated by taking
    // curl of the Jacobian matrix dx_ds evaluated at the quadrature point.
    // Unit outward normal vector and its norm
    // This group of get-functions will return the unit normal vector as well as
    // the norm of r_s x r_t (r_s x r_u, r_t x r_u ...). The returned norm give 
    // the measure of the surface element, which is necessary for boundary
    // integral.
    //    Input: the quadrature point index.
    //    Output: Cartesian components of the 3D normal vector.
    // ------------------------------------------------------------------------
    virtual void get_3d_normal_bottom( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const
    {SYS_T::commPrint("Warning: get_3d_normal_bottom is not implemented. \n");}

    virtual void get_3d_normal_top( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const
    {SYS_T::commPrint("Warning: get_3d_normal_top is not implemented. \n");}

    virtual void get_3d_normal_left( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const
    {SYS_T::commPrint("Warning: get_3d_normal_left is not implemented. \n");}

    virtual void get_3d_normal_right( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const
    {SYS_T::commPrint("Warning: get_3d_normal_right is not implemented. \n");}

    virtual void get_3d_normal_front( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const
    {SYS_T::commPrint("Warning: get_3d_normal_front is not implemented. \n");}

    virtual void get_3d_normal_back( const int &quaindex,
        double &nx, double &ny, double &nz, double &surfaceArea ) const
    {SYS_T::commPrint("Warning: get_3d_normal_back is not implemented. \n");}

    // ------------------------------------------------------------------------
    // Unit outward normal vector and boundary line measure for 2D
    // In a 2D quad element, there are four faces:
    //              left (s, 0)     right (s, 1)
    //              back (0, t)     front (1, t)
    // in a (s,t) in [0,1]^2 square.
    // The ourward normal are calculated and the boundary line length is
    // calculated.
    // ------------------------------------------------------------------------
    virtual void get_2d_normal_front( const int &quaindex,
        double &nx, double &ny, double &line ) const
    {SYS_T::commPrint("Warning: get_2d_normal_front is not implemented. \n");}

    virtual void get_2d_normal_back( const int &quaindex,
        double &nx, double &ny, double &line ) const
    {SYS_T::commPrint("Warning: get_2d_normal_back is not implemented. \n");}

    virtual void get_2d_normal_left( const int &quaindex,
        double &nx, double &ny, double &line ) const
    {SYS_T::commPrint("Warning: get_2d_normal_left is not implemented. \n");}

    virtual void get_2d_normal_right( const int &quaindex,
        double &nx, double &ny, double &line ) const
    {SYS_T::commPrint("Warning: get_2d_normal_right is not implemented. \n");}


    // ------------------------------------------------------------------------
    // Unit outward normal vector and boundary area measure.
    // In a oriented surface, the outward normal direction and the area
    // measure can be calculated using the cross product of vector x_,xi and
    // vector x_,eta. Assuing that the nodal points are properly ordered, this
    // product will point in the outward direction with respect to the volume
    // element. To ensure the ordering, some special procedure may need to be
    // performed in the preprocessor to check and order the nodal indices.
    // See, ElemBC_3D_tet4::resetTriIEN_outwardnormal() function for example
    // for the ordering of a triangle element that belongs to a tetrahedron
    // surface.
    // This function is called in FEAElement_Triangle3_3D_der0.
    // ------------------------------------------------------------------------
    virtual Vector_3 get_2d_normal_out( const int &quaindex, double &area ) const
    {
      SYS_T::commPrint("Warning: get_2d_normal_out is not implemented. \n");
      return Vector_3();
    }

    // ------------------------------------------------------------------------
    // Unit outward normal vector and boundary area measure.
    // In a line element, the tangential vector is well defined. Given
    // an interior point's coordinate, the outward normal vector is defined.
    // intpt_x/y/z gives the interior point's coordinate.
    // quaindex gives the index of the quadrature point, from which one can
    // calculate the "root" of the tangential vector.
    // This function is, for example, called in FEAElement_Line2_3D_der0.
    // ------------------------------------------------------------------------
    virtual void get_normal_out( const int &quaindex,
        const double &sur_pt_x, const double &sur_pt_y,
        const double &sur_pt_z, const double &intpt_x, 
        const double &intpt_y, const double &intpt_z,
        double &nx, double &ny, double &nz, double &length ) const
    {SYS_T::commPrint("Warning: get_normal_out is not implemented. \n");}


    // ------------------------------------------------------------------------
    // Reset function
    // Reset the parameters of the element object. 
    // These functions are utilized when one construct the element class once
    // and reset the element properties over the element loop. In this way, the
    // container are never destructed during the element loop (potentially this
    // is beneficial for the memory pool.) The values are recomputed repeatedly
    // during the buildBasis calls. 
    // ------------------------------------------------------------------------
    virtual void reset_degree( const int &new_sdeg, const int &new_tdeg )
    {SYS_T::commPrint("Warning: reset_degree is not implemented. \n");}

    virtual void reset_degree( const int &new_sdeg, const int &new_tdeg, const int &new_udeg )
    {SYS_T::commPrint("Warning: reset_degree is not implemented. \n");}

    virtual void reset_nLocBas( const int &new_nlocbas )
    {SYS_T::commPrint("Warning: reset_nLocBas is not implemented. \n");}

    virtual void reset_numQua( const int &new_squa, const int &new_tqua )
    {SYS_T::commPrint("Warning: reset_numQua is not implemented. \n");}

    virtual void reset_numQua( const int &new_squa, const int &new_tqua, const int &new_uqua )
    {SYS_T::commPrint("Warning: reset_numQua is not implemented. \n");}

    virtual void reset_numQuapts( const int &new_numQuapts )
    {SYS_T::commPrint("Warning: reset_numQuapts is not implemented. \n");}

    virtual void resize_container()
    {SYS_T::commPrint("Warning: resize_container() is not implemented. \n");}

};

#endif
