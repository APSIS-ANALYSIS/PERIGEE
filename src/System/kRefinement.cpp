#include "kRefinement.hpp"

void hRefine_newKnot_Generator( const std::vector<double> &existKnot,
    std::vector<double> &insertKnot, const int num_inserted,
    double &h_max, double &h_min )
{
  // This flag is for debug use
  bool isPrint = false;

  // Initialize insertKnot vector
  insertKnot.clear();
  insertKnot.reserve(num_inserted);

  // Check the given vector format
  if(existKnot.size()<2)
  {
    std::cerr<<"Error: The knot vector size is less than 2."<<std::endl;
    exit(1);
  }
  if( existKnot[0] != 0.0 )
  {
    std::cerr<<"ERROR: The knot vector does not start with 0.0, but with "<<existKnot[0]<<std::endl;
    exit(1);
  }
  if( existKnot.back() != 1.0 )
  {
    std::cerr<<"ERROR: The knot vector does not end with 1.0, but with "<<existKnot.back()<<std::endl;
    exit(1);
  }

  // count the num of knot spans and corresponding lengths in existKnot
  int num_knot_span = 0;
  std::vector<double> knot_span_length;
  knot_span_length.clear();
  for(unsigned int ii = 0; ii<existKnot.size()-1; ++ii)
  {
    if( existKnot[ii] > existKnot[ii+1] )
    {
      std::cerr<<"ERROR: The knot vector is not in ascending order at "
        <<ii<<std::endl;
      exit(1);
    }

    if( existKnot[ii] < existKnot[ii+1] )
    {
      num_knot_span += 1;
      knot_span_length.push_back(existKnot[ii+1] - existKnot[ii]);
    }
  }
  VEC_T::shrink2fit( knot_span_length );

  if( isPrint )
  {
    std::cout<<"# of know spans: "<<num_knot_span<<std::endl;
    for(int ii=0; ii<(int)knot_span_length.size(); ++ii)
      std::cout<<knot_span_length[ii]<<'\t';
    std::cout<<std::endl;
  }


  // quasi-even distrubtion of insert knot in each span
  if(num_inserted < num_knot_span)
  {
    std::cerr<<"ERROR: There are "<<num_inserted<<" new knots to be inserted into "
      <<num_knot_span<<" knot spans."<<std::endl;
    exit(1); 
  }

  int inst_res = num_inserted % num_knot_span;
  int inst_step = (num_inserted - inst_res) / num_knot_span;
  std::vector<int> inserted_knot_in_span; 
  for(int ii=0; ii<num_knot_span; ++ii)
  {
    if( ii < inst_res )
      inserted_knot_in_span.push_back(inst_step + 1);
    else
      inserted_knot_in_span.push_back( inst_step );
  }
  VEC_T::shrink2fit(inserted_knot_in_span);


  if(isPrint)
  {
    for(int ii=0; ii<num_knot_span; ++ii)
      std::cout<<inserted_knot_in_span[ii]<<'\t';
    std::cout<<std::endl;
  }

  double span_start_knot = 0.0;
  h_max = 0.0; h_min = 1.0;
  for(int ii=0; ii<num_knot_span; ++ii)
  {
    double span_h_insert = knot_span_length[ii] / (1+inserted_knot_in_span[ii]);
    
    if(span_h_insert > h_max)
      h_max = span_h_insert;
    if(span_h_insert < h_min)
      h_min = span_h_insert;
    
    for(int jj=1; jj<1+inserted_knot_in_span[ii]; ++jj)
    {
      insertKnot.push_back( span_start_knot + jj * span_h_insert );
    }
    span_start_knot += knot_span_length[ii];
  }
  
  VEC_T::shrink2fit( insertKnot );
}


void hRefine_newKnot_Generator( const std::vector<double> &existKnot,
    std::vector<double> &insertKnot, const std::vector<int> &num_inserted,
    double &h_max, double &h_min )
{
  // This flag is for debug use
  bool isPrint = false;

  // Initialize insertKnot vector
  insertKnot.clear();

  // Check the given vector format
  if(existKnot.size()<2)
  {
    std::cerr<<"Error: The knot vector size is less than 2."<<std::endl;
    exit(1);
  }
  if( existKnot[0] != 0.0 )
  {
    std::cerr<<"ERROR: The knot vector does not start with 0.0, but with "<<existKnot[0]<<std::endl;
    exit(1);
  }
  if( existKnot.back() != 1.0 )
  {
    std::cerr<<"ERROR: The knot vector does not end with 1.0, but with "<<existKnot.back()<<std::endl;
    exit(1);
  }
  for(unsigned int ii=0; ii<num_inserted.size(); ++ii)
  {
    if(num_inserted[ii]<0)
    {
      std::cerr<<"ERROR: The num_inserted is in wrong format at "<<
        ii<<" with value "<<num_inserted[ii]<<std::endl;
      exit(1);
    }
  }

  // count the num of knot spans and corresponding lengths in existKnot
  int num_knot_span = 0;
  std::vector<double> knot_span_length;
  knot_span_length.clear();
  for(unsigned int ii = 0; ii<existKnot.size()-1; ++ii)
  {
    if( existKnot[ii] > existKnot[ii+1] )
    {
      std::cerr<<"ERROR: The knot vector is not in ascending order at "
        <<ii<<std::endl;
      exit(1);
    }

    if( existKnot[ii] < existKnot[ii+1] )
    {
      num_knot_span += 1;
      knot_span_length.push_back(existKnot[ii+1] - existKnot[ii]);
    }
  }
  VEC_T::shrink2fit( knot_span_length );

  if( isPrint )
  {
    std::cout<<"# of know spans: "<<num_knot_span<<std::endl;
    for(int ii=0; ii<(int)knot_span_length.size(); ++ii)
      std::cout<<knot_span_length[ii]<<'\t';
    std::cout<<std::endl;
  }

  if( (int)num_inserted.size() != num_knot_span )
  {
    std::cerr<<"ERROR: The size of given inserted knot number vector"
      <<" does not match the number of knot spans in exsiting knot vector"
      <<std::endl;
    exit(1);
  }

  double span_start_knot = 0.0;
  h_max = 0.0; h_min = 1.0;
  for(int ii=0; ii<num_knot_span; ++ii)
  {
    double span_h_insert = knot_span_length[ii] / (1+num_inserted[ii]);
    
    if(span_h_insert > h_max)
      h_max = span_h_insert;
    if(span_h_insert < h_min)
      h_min = span_h_insert;
    
    for(int jj=1; jj<1+num_inserted[ii]; ++jj)
    {
      insertKnot.push_back( span_start_knot + jj * span_h_insert );
    }
    span_start_knot += knot_span_length[ii];
  }

  VEC_T::shrink2fit( insertKnot );
}

int knotVec_check(const std::vector<double> &knotVec,
   const int degree )
{
  int span_num = 0;
  std::vector<double> copy_knotVec(knotVec);
 
  if(degree < 1)
  {
    std::cerr<<"ERROR: The given degree is less than 1."<<std::endl;
    return 0;
  }

  if((int)copy_knotVec.size() < 2*(degree+1) )
  {
    std::cerr<<"ERROR: The knot vector size is less than 2*p+2."<<std::endl;
    return 0;
  }

  if(copy_knotVec[0] != 0.0)
  {
    std::cerr<<"ERROR: The knot vector does not start with 0.0."<<std::endl;
    return 0;
  }
  if(copy_knotVec[copy_knotVec.size()-1] != 1.0)
  {
    std::cerr<<"ERROR: The knot vector does not end with 1.0"<<std::endl;
    return 0;
  }

  if(copy_knotVec[degree+1] == 0.0)
  {
    std::cerr<<"ERROR: The given degree is wrong."<<std::endl;
    return 0;
  }
  if(copy_knotVec[copy_knotVec.size()-2-degree] == 1.0)
  {
    std::cerr<<"ERROR: The given degree is wrong."<<std::endl;
    return 0;
  }

  for(int ii=1; ii<degree+1; ++ii)
  {
    if(copy_knotVec[ii] != 0.0)
    {
      std::cerr<<"ERROR: The given degree is wrong."<<std::endl;
      return 0;
    }
    if(copy_knotVec[copy_knotVec.size()-1-ii] != 1.0)
    {
      std::cerr<<"ERROR: The given degree is wrong."<<std::endl;
      return 0;
    }
  }

  std::vector<double>::iterator it;
  it = unique(copy_knotVec.begin(), copy_knotVec.end());
  copy_knotVec.resize( distance( copy_knotVec.begin(), it ) );

  for(unsigned int ii=0; ii<copy_knotVec.size()-1; ++ii)
  {
    if(copy_knotVec[ii] >= copy_knotVec[ii+1])
    {
      std::cerr<<"ERROR: The knot vector is not in ascending order."<<std::endl;
      return 0;
    }
  }


  span_num = copy_knotVec.size() - 1; 

  return span_num;
}

void kRefinement( const int &addSDegree, 
    const int &addTDegree, const int &addUDegree,
    const std::vector<double> &insertSKnot,
    const std::vector<double> &insertTKnot,
    const std::vector<double> &insertUKnot,
    std::vector<double> &sknots, 
    std::vector<double> &tknots, 
    std::vector<double> &uknots,
    std::vector<double> &ctrlPts, 
    const int &dim, int &sdegree, int &tdegree, int &udegree )
{
  // s-dir order elevation
  if(addSDegree >0)
  {
    if(NURBS_T::degreeElevateVolume( sknots, sdegree, tknots, tdegree,
          uknots, udegree, dim, ctrlPts, addSDegree, 's') != 0)
    {
      std::cerr<<"ERROR: Failed to elevate the degree in s direction. \n";
      exit(1);
    }
    sdegree = sdegree + addSDegree;
  }
  
  // t-dir order elevation
  if(addTDegree >0)
  {
    if(NURBS_T::degreeElevateVolume( sknots, sdegree, tknots, tdegree,
          uknots, udegree, dim, ctrlPts, addTDegree, 't') != 0)
    {
      std::cerr<<"ERROR: Failed to elevate the degree in t direction. \n";
      exit(1);
    }
    tdegree = tdegree + addTDegree;
  }

  // u-dir order elevation
  if(addUDegree >0)
  {
    if(NURBS_T::degreeElevateVolume( sknots, sdegree, tknots, tdegree,
          uknots, udegree, dim, ctrlPts, addUDegree, 'u') != 0)
    {
      std::cerr<<"ERROR: Failed to elevate the degree in u direction. \n";
      exit(1);
    }
    udegree = udegree + addUDegree;
  }

  // s-dir knot insertion
  if( insertSKnot.size() > 0)
  {
    if( NURBS_T::knotRefinementVolume( sknots, sdegree, tknots, tdegree, uknots,
          udegree, ctrlPts, insertSKnot, 's', dim ) != 0 )
    {
      std::cerr<<"ERROR: Failed to insert knots in s direction. \n";
      exit(1);
    }
  }

  // t-dir knot insertion
  if( insertTKnot.size() > 0)
  {
    if( NURBS_T::knotRefinementVolume( sknots, sdegree, tknots, tdegree, uknots,
          udegree, ctrlPts, insertTKnot, 't', dim ) != 0 )
    {
      std::cerr<<"ERROR: Failed to insert knots in t direction. \n";
      exit(1);
    }
  }

  // u-dir knot insertion
  if( insertUKnot.size() > 0)
  {
    if( NURBS_T::knotRefinementVolume( sknots, sdegree, tknots, tdegree, uknots,
          udegree, ctrlPts, insertUKnot, 'u', dim ) != 0 )
    {
      std::cerr<<"ERROR: Failed to insert knots in u direction. \n";
      exit(1);
    }
  }

}


void kRefinement( const int &addSDegree, const int &addTDegree,
    const std::vector<double> &insertSKnot,
    const std::vector<double> &insertTKnot,
    std::vector<double> &sknots, std::vector<double> &tknots,
    std::vector<double> &ctrlPts, 
    const int &dim, int &sdegree, int &tdegree )
{
  // s-dir order elevation
  if(addSDegree >0)
  {
    if(NURBS_T::degreeElevateSurface( sknots, sdegree, tknots, tdegree,
          dim, ctrlPts, addSDegree, 's') != 0)
    {
      std::cerr<<"ERROR: Failed to elevate the degree in s direction. \n";
      exit(1);
    }
    sdegree = sdegree + addSDegree;
  }
  
  // t-dir order elevation
  if(addTDegree >0)
  {
    if(NURBS_T::degreeElevateSurface( sknots, sdegree, tknots, tdegree,
          dim, ctrlPts, addTDegree, 't') != 0)
    {
      std::cerr<<"ERROR: Failed to elevate the degree in t direction. \n";
      exit(1);
    }
    tdegree = tdegree + addTDegree;
  }


  // s-dir knot insertion
  if( insertSKnot.size() > 0)
  {
    if( NURBS_T::knotRefinementSurface( sknots, sdegree, tknots, tdegree,
          dim, ctrlPts, insertSKnot, 's' ) != 0 )
    {
      std::cerr<<"ERROR: Failed to insert knots in s direction. \n";
      exit(1);
    }
  }

  // t-dir knot insertion
  if( insertTKnot.size() > 0)
  {
    if( NURBS_T::knotRefinementSurface( sknots, sdegree, tknots, tdegree,
          dim, ctrlPts, insertTKnot, 't' ) != 0 )
    {
      std::cerr<<"ERROR: Failed to insert knots in t direction. \n";
      exit(1);
    }
  }

}


void pRefinement(  const int &dim, const int &addSDegree,
    const int &addTDegree, const int &addUDegree,
    std::vector<double> &sKnots,
    std::vector<double> &tKnots,
    std::vector<double> &uKnots,
    std::vector<double> &ctrlPts,
    int &sdegree, int &tdegree, int &udegree )
{
  SYS_T::print_fatal_if(addSDegree != addTDegree || addSDegree != addUDegree || addTDegree != addUDegree, "Error: pRefinement, add degree in s-t-u directions are not equal.\n");

  // s-dir order elevation
  if( addSDegree > 0 )
  {
    if(NURBS_T::degreeElevateVolume( sKnots, sdegree, tKnots, tdegree,
          uKnots, udegree, dim, ctrlPts, addSDegree, 's') != 0)
    {
      std::cerr<<"ERROR: Failed to elevate the degree in s direction. \n";
      exit(1);
    }
    sdegree = sdegree + addSDegree;
  }

  // t-dir order elevation
  if(addTDegree >0)
  {
    if(NURBS_T::degreeElevateVolume( sKnots, sdegree, tKnots, tdegree,
          uKnots, udegree, dim, ctrlPts, addTDegree, 't') != 0)
    {
      std::cerr<<"ERROR: Failed to elevate the degree in t direction. \n";
      exit(1);
    }
    tdegree = tdegree + addTDegree;
  }

  // u-dir order elevation
  if(addUDegree >0)
  {
    if(NURBS_T::degreeElevateVolume( sKnots, sdegree, tKnots, tdegree,
          uKnots, udegree, dim, ctrlPts, addUDegree, 'u') != 0)
    {
      std::cerr<<"ERROR: Failed to elevate the degree in u direction. \n";
      exit(1);
    }
    udegree = udegree + addUDegree;
  }
}

// EOF
