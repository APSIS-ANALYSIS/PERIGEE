#ifndef POST_PROBE_HPP
#define POST_PROBE_HPP
// ============================================================================
// Post_probe.hpp
// ----------------------------------------------------------------------------
// This is the header file for the postprocess -- numerical probe object.
// This object is utilized to calculate quantities' values at given locations
// and time.
//
// Date: Jan. 14 2016
// ============================================================================
#include "HDF5_Reader.hpp"
#include "QuadPts_Gauss.hpp"
#include "QuadPts_debug.hpp"
#include "AInt_Weight.hpp"
#include "FEAElement.hpp"
#include "ALocal_meshSize_2D_NURBS.hpp"
#include "AExtractor_2D_NURBS_xy.hpp"

class Post_probe
{
  public:
    Post_probe( const std::string &part_bname_in,
        const std::string &sol_bname_in,
        const int &in_sdeg, const int &in_tdeg, const int &in_udeg,
        const int &nlocbas, const int &in_nfunc, const int &in_dof );

    ~Post_probe();

    void print_info() const;


    // ------------------------------------------------------------------------
    // read the PETSc binary solution vector;
    // check the solution vector length, if compatible with the mesh file, then
    // return a double vector containing the full solution vection
    // ------------------------------------------------------------------------
    void readPETSc_full( const int &sol_index, std::vector<double> &out ) const;


    // ------------------------------------------------------------------------
    // Perform line integral
    // ------------------------------------------------------------------------
    void get_lineAve_forNu( const int &sol_index,
        const std::vector<double> &eta_list,
        const std::vector<int> &ey_list,
        const int &in_nqps, const int &nElem_x, 
        const int &nElem_y, FEAElement * const &element,
        std::vector<double> &uzt, std::vector<double> &the, std::vector<double> &yyy ) const;


    // ------------------------------------------------------------------------
    // Perform line integral over a time-averaged solution
    // ------------------------------------------------------------------------
    void get_lineAve_forNu( const int &sol_index_start,
        const int &sol_index_step, const int &sol_index_end,
        const std::vector<double> &eta_list,
        const std::vector<int> &ey_list,
        const int &in_nqps, const int &nElem_x, 
        const int &nElem_y, FEAElement * const &element,
        std::vector<double> &uzt, std::vector<double> &the, std::vector<double> &yyy ) const;


    // ------------------------------------------------------------------------
    // Perform line integral for a TNSK solution to obtain density and
    // temperature, and the y-coordinate.
    // ------------------------------------------------------------------------
    void get_lineAve_NSK_rho_theta( const int &sol_index,
        const std::vector<double> &eta_list,
        const std::vector<int> &ey_list,
        const int &in_nqps, const int &nElem_x, 
        const int &nElem_y, FEAElement * const &element,
        std::vector<double> &rho, std::vector<double> &the, 
        std::vector<double> &pre, std::vector<double> &yyy ) const;


    // ------------------------------------------------------------------------
    // Given the element index and the \xi coordinate in the reference domain
    // [0,1]^2 (xi_x, xi_y), interpolate the quantities' value val and the 
    // physical coordinate (xx,yy).
    // ------------------------------------------------------------------------
    void get_val_cood( const int sol_index,
        const int &ee,
        const double &xi_s, const double &xi_t, 
        FEAElement * const &element,
        std::vector<double> &val, double &xx, double &yy ) const;

  private:
    const std::string partfile_bname;
    const std::string sol_bname;

    const int s_degree, t_degree, u_degree;
    const int nLocBas;
    const int nFunc;
    const int dof;

    int * elepart;
    int nelement;
   
    // ---------- Private Functions ----------
    // ------------------------------------------------------------------------
    // Given the element global index, determine which cpu does it belong to
    // based on the epart info
    // ------------------------------------------------------------------------
    int get_elem_cpu(const int &ee ) const;


    // ------------------------------------------------------------------------
    // read the PETSc binary solution vector and return the solution
    // corresponding to gloIEN 
    // the veccopy pointer should allocate enough space. 
    // ------------------------------------------------------------------------
    void readPETSc_vec( const int &sol_index,
        const std::vector<int> &gloIEN,
        double * const &veccopy ) const;
    
    // ------------------------------------------------------------------------
    // read a list of PETSc binary solution vectors and return their arithemtic
    // averages.
    // out is double vector that returns the 
    //    (SOL_sol_index_start + ... + SOL_sol_index_end) / N
    // ------------------------------------------------------------------------
    void readPETSc_average( const int &sol_index_start, const int &sol_index_step,
      const int &sol_index_end, std::vector<double> &out  ) const;

};

#endif
