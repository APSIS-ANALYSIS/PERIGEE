#ifndef GENBC_CORONARY_HPP
#define GENBC_CORONARY_HPP
// ============================================================================
// GenBC_Coronary.hpp
// This class is used to provide outflow boundary conditions for an
// open loop coronary lumped parameter model (LPM). The RCR model is also
// handled by this class in case coronary and RCR BCs are used together.
// The coronary artery is modeled by an RCRCR circuit with intramyocardial
// pressures:
//
// Ra---Ramicro---Rv---Pd
//    |         |
//   Ca        Cim
//              |
//             Pim
//
// Ra, Ramicro and Rv are resistors for the coronary arteries, coronary 
// microvasculature and coronary veins, respectively.
// Ca and Cim are capacitors for proximal and distal vasculature, 
// respectively.
// Pd and Pim are the distal and intramyocardial pressures, respectively.
// Intramyocardial pressure Pim is applied to capacitor Cim to model 
// restricted coronary flow during systole.
//
// The GenBC_Coronary input file uses the following format,
// -----------------------
// Coronary <num_outlets>
// <face_id0> <Ra> <Ca> <Ramicro> <Cim> <Rv> <Pd> <num_Pim_data> <Pim_scaling_coef>
// <t1> <P1>
// <t2> <P2>
//   ....
// <tn> <Pn>
// <face_id1> <Ra> <Ca> <Ramicro> <Cim> <Rv> <Pd> <num_Pim_data> <Pim_scaling_coef>
// <t1> <P1>
// <t2> <P2>
//   ....
// <tn> <Pn>
//-----------------------
// where n=num_Pim_data. Pim_scaling_coef can be used to scale original pressure Pim.
// If an RCR is applied to an outlet, set num_Pim_data=0 and omit Pim data as follows,
//
// <outlet_id> <Rp> <C> <Rd> <C> <Rd> <Pd> 0 1.0
//
// ============================================================================
#include <iomanip>
#include "IGenBC.hpp"
#include "GenBC_Tools.hpp"

class GenBC_Coronary : public IGenBC
{
  public:
    GenBC_Coronary( const std::string &lpn_filename, const int &in_N, const double &dt3d,
        const int &in_index, const std::string &in_lpn_sol_file="lpn_coronary_sol.txt" );

    virtual ~GenBC_Coronary() = default;

    virtual void print_info() const;

    virtual int get_num_ebc() const { return num_ebc; }

    virtual double get_m( const int &ii, const double &in_dot_Q, const double &in_Q ) const;

    virtual double get_n( const int &ii, const double &in_dot_Q, const double &in_Q ) const 
    { return 0.0; }

    // Obtain P for the ii-th outlet surface (coronary or RCR).
    virtual double get_P( const int &ii, const double &in_dot_Q, const double &in_Q,
       const double &time = 0.0 ) const;

    // Get initial P for the ii-th outlet face at the beginninging of LPM ODE integration.
    // 0 <= ii < num_ebc
    virtual double get_P0( const int &ii ) const;

    // Set initial values for the LPM ODE integration.
    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
       const double &in_P_0, const double &curr_time, const bool &is_restart );

    // Write 0D solutions into a file for restart
    virtual void write_0D_sol( const int &curr_index, const double &curr_time ) const;

  private:
    const int num_odes; // Number of ODEs in the model

    const int N; // ODE integrator's number of time steps

    const double h; // delta t = Nh

    // file to store 0D solutions at each 3D time step
    const std::string lpn_sol_file;

    // Parameters used to define difference quotient for get_m.
    const double absTol, relTol;

    // Total number of outlet surfaces
    int num_ebc;

    // Vectors storing the Ra, Ca, Ra_micro, Cim, Rv, Pd, and alpha_Pim.
    // alpha_Pim stores the scaling values for all coronary outlet faces. 
    // The length of the vectors is num_ebc
    std::vector<double> Ra, Ca, Ra_micro, Cim, Rv, Pd, alpha_Pim;

    // Number of intramyocardial pressure Pim data points for each outlet face
    // The vector length is num_ebc.
    // Note: num_Pim_data=0 indicates an RCR outlet.
    std::vector<int> num_Pim_data;

    // Time_data and Pim_data for user-provided intramyocardial pressure 
    // waveform (time-pressure) for each coronary outlet face
    // der_Pim_data stands for the corresponding dPim/dt for each coronary 
    // outlet face.
    // Their sizes are num_ebc x num_Pim_data[ii] with 0 <= ii < num_ebc 
    std::vector< std::vector<double> > Time_data, Pim_data, der_Pim_data;

    // precomputed dPim/dt needed by RK4, 
    // dPimdt_k1 has size num_ebc x N+1
    // dPimdt_k2/3 has size num_ebc x N
    std::vector< std::vector<double> > dPimdt_k1, dPimdt_k2, dPimdt_k3;

    // prev_0D_sol records solutions when each ODE integration is completed.
    mutable std::vector< std::vector<double> > prev_0D_sol;

    // restart_0D_sol stores 0D solutions read from lpn_sol_file for a
    // restart job from a previous time step.
    std::vector< std::vector<double> > restart_0D_sol;

    // Vectors for outlet initial flow and capacitor pressures (2 capacitors)
    // The vector length is num_ebc
    std::vector<double> Q0;

    // The size of Pi0 is 2 x num_ebc
    std::vector< std::vector<double> > Pi0;

    // PRIVATE FUNCTIONS:
    // Evaluate the coronary LPM (2 first order ODEs) for the ii-th outlet 
    // face (which is a coronary outlet) and output the ODE derivatives to K.
    void F_coronary( const int &ii, const std::vector<double> &pi, const double &q, 
        const double &dPimdt, std::vector<double> &K ) const;

    // Evaluate the RCR ODE for the ii-th outlet face and return the result.
    double F_RCR( const int &ii, const double &pi, const double &q ) const;

    // Pre-compute dPim/dt at the begining of the ODE integration for the 
    // ii-th outlet face, 0 <= ii < num_ebc
    void get_dPim_dt( const int &ii, const double &time_start, const double &time_end );
};

#endif
