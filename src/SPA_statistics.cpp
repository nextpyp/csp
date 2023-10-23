#include <iostream>
#include <fstream>
#include <limits>
#include <vector>

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array/stencil-et.h>

using namespace blitz;

int main( int argc, char ** argv ) {
  // if (argc != 25) {
    // cout << "USAGE: " << argv[0]
        // << " tilt_angle axis_rotation axis_correction box_correction[2] normal[0] normal[1] normal[2] m[0-15] \n"
        // << endl;
    // return 0;
  // }
  
	ifstream ifs( argv[1] );	
    if (ifs.bad())
    {
        cerr << "Unable to open file: " << argv[1] << endl;
        exit(1);
    }

	Array< float, 2 > A( atoi(argv[2]), atoi(argv[3]) );
	Array< float, 2 > :: iterator iter = A.begin();
	
	if (ifs.is_open())
  {
    while ( ifs.good() )
    {
	  float number;
	  ifs >> number;
		*iter = number;
		iter++;
    }
    ifs.close();
  }
  
	firstIndex i;
	secondIndex j;
	cout << A << endl;
	cout << mean(A(Range::all(),0)) << endl;
	cout << mean(A(Range::all(),1)) << endl;
	Array< float, 1 > means( mean(A(j,i),j) );
	
	// sqrt( sum( pow2( T_bin4 - meanInside ) ) / T_bin4.size() );
	
	// Array< float, 1 > variances( mean(A(j,i),j) );
 	// cout << B << endl;
	
  
  // printf( "%.2f %.2f %.2f %.2f %.2f\n", frealign[0], frealign[1], frealign[2], shift[0], shift[1] );

  return 0;
}
