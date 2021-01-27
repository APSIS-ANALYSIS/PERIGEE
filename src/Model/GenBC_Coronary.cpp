#include "GenBC_Coronary.hpp"

GenBC_Coronary::GenBC_Coronary( const char * const &lpn_filename, 
    const int &in_N, const double &dt3d )
: N( in_N ), h( dt3d/static_cast<double>(N) ),
  absTol( 1.0e-8 ), relTol( 1.0e-5 )
{
  // Now read the lpn input file for num_ebc and coronary model 
  // parameters (Ra, Ca, Ra_micro, Cim, Rv, Pd and Pim)
  std::string temp_name( lpn_filename );
  SYS_T::file_check( temp_name ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename, std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  std::string bc_type;
  tstart=0.0;
  tend=N*h;
  
  // The first non-commented line should be
  // Coronary num_ebc
  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str(sline);
      sstrm >> bc_type;
      sstrm >> num_ebc;
      sstrm.clear();
      break;
    }
  }

  // Check the file's bc_type matches Coronary
  if( bc_type.compare("Coronary") ==0
      || bc_type.compare("CORONARY") == 0
      || bc_type.compare("coronary") == 0 )
  {
    Ra.resize( num_ebc ); 
    Ca.resize( num_ebc );
    Ra_micro.resize( num_ebc ); 
    Cim.resize( num_ebc );
    Rv.resize( num_ebc ); 
    Pd.resize( num_ebc );
    alpha_Pim.resize( num_ebc );
    Q0.resize( num_ebc );
    Pi0.resize( num_ebc ); 
    num_Pimdata.resize( num_ebc );
    tdata.resize( num_ebc );
    Pimdata.resize( num_ebc );
    Pimderdata.resize( num_ebc );
    prev_0D_sol.resize( num_ebc );
    dPimdt_k1.resize( num_ebc );
    dPimdt_k2.resize( num_ebc );
    dPimdt_k3.resize( num_ebc );

    for(int ii =0;ii<num_ebc;++ii)
    {
      prev_0D_sol[ii].resize(2);
      Pi0[ii].resize(2);
    }

    for(int ii=0;ii<num_ebc;++ii)
    {
      dPimdt_k1[ii].resize(N+1);
      dPimdt_k2[ii].resize(N);
      dPimdt_k3[ii].resize(N);
    }
  }
  else SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Coronary.\n", lpn_filename);

  // Read files for each ebc to set the values of Ra, Ca, Ra_micro, 
  // Cim,Rv, Pd Pim and alpha_Pim
  int counter = 0;
  int data_size=1;
  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str( sline );
      int face_id;
      sstrm >> face_id;

      // Make sure the face_id, the first column in the file are listed
      // from 0 to ebc_id - 1
      SYS_T::print_fatal_if( face_id != counter, "Error: GenBC_Coronary the input file %s has wrong format in the face id column (the first column). \n", lpn_filename);

      sstrm >> Ra[ counter ];
      sstrm >> Ca[ counter ];
      sstrm >> Ra_micro[ counter ];
      sstrm >> Cim[ counter ];
      sstrm >> Rv[ counter ];
      sstrm >> Pd[ counter ];
      sstrm >> num_Pimdata[ counter ];
      sstrm >> alpha_Pim[ counter ];

      SYS_T::print_fatal_if(num_Pimdata[counter]<=2 && num_Pimdata[counter]!=0, 
          "Error: num of Pim data needs to be 0 for RCR or >2 for coronary  \n");

      if(num_Pimdata[counter]>0)
      {
       data_size=num_Pimdata[counter];
      }
      else
      {
       data_size=1;
      }

      tdata[counter].resize(data_size);
      Pimdata[counter].resize(data_size);
      Pimderdata[counter].resize(data_size);

      sstrm.clear();

      for (int ii =0;ii<num_Pimdata[counter];++ii)
      {
        getline(reader, sline);
        sstrm.str( sline );
        sstrm>>tdata[counter][ii];

        sstrm>>Pimdata[counter][ii];
        Pimdata[counter][ii]=Pimdata[counter][ii]*alpha_Pim[counter];
        sstrm.clear();
      }

      if(num_Pimdata[counter]>0)
      {
        set_phcip(counter);
        get_dPimdt(counter);
      }

      SYS_T::print_fatal_if(tdata[counter][0]>0.0, "Error: Pim data do not start from 0.\n");
      counter += 1;
    }
  }

  if(counter != num_ebc ) SYS_T::print_fatal("Error: GenBC_Coronary the input file %s does not contain complete data for outlet faces. \n", lpn_filename);

  reader.close();

  SYS_T::commPrint( "===> GenBC_Coronary data are read in from %s.\n", lpn_filename );

  // Set a zero initial value. They should be reset based on the initial
  // 3D solutions.
  for(int ii=0; ii<num_ebc; ++ii)
  {
    Q0[ii] = 0.0;
    Pi0[ii][0] = 0.0;
    Pi0[ii][1]=0.0;
    prev_0D_sol[ii][0]=0.0;
    prev_0D_sol[ii][1]=0.0;

    // Make sure C and R are nonzero
    SYS_T::print_fatal_if(Ca[ii]==0.0, "Error: GenBC_Coronary Ca cannot be zero.\n");
    SYS_T::print_fatal_if(Cim[ii]==0.0, "Error: GenBC_Coronary Cim cannot be zero.\n");
    SYS_T::print_fatal_if(Ra_micro[ii]==0.0, "Error: GenBC_Coronary Ra_micro cannot be zero.\n");
    SYS_T::print_fatal_if(Rv[ii]==0.0, "Error: GenBC_Coronary Rv cannot be zero.\n");
  }
}

GenBC_Coronary::~GenBC_Coronary()
{}

void GenBC_Coronary::print_info() const
{
  SYS_T::commPrint( "     Coronary model: N = %d, h = %e, num_ebc = %d \n", N, h, num_ebc );

  for(int ii=0; ii<num_ebc; ++ii)
    SYS_T::commPrint( "     ebcid = %d, Ra = %e, Ca = %e, Ra_micro = %e,Rim = %e, Rv = %e, Pd = %e \n", ii, Ra[ii], Ca[ii],Ra_micro[ii],Cim[ii], Rv[ii], Pd[ii] );
}


// EOF
