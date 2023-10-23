/* randn()
* 
* Normally (Gaussian) distributed random numbers, using the Box-Muller 
* transformation.  This transformation takes two uniformly distributed deviates
* within the unit circle, and transforms them into two independently 
* distributed normal deviates.  Utilizes the internal rand() function; this can
* easily be changed to use a better and faster RNG.
* 
* The parameters passed to the function are the mean and standard deviation of 
* the desired distribution.  The default values used, when no arguments are 
* passed, are 0 and 1 - the standard normal distribution.
* 
* 
* Two functions are provided:
* 
* The first uses the so-called polar version of the B-M transformation, using
* multiple calls to a uniform RNG to ensure the initial deviates are within the
* unit circle.  This avoids making any costly trigonometric function calls.
* 
* The second makes only a single set of calls to the RNG, and calculates a 
* position within the unit circle with two trigonometric function calls.
* 
* The polar version is generally superior in terms of speed; however, on some
* systems, the optimization of the math libraries may result in better 
* performance of the second.  Try it out on the target system to see which
* works best for you.  On my test machine (Athlon 3800+), the non-trig version
* runs at about 3x10^6 calls/s; while the trig version runs at about
* 1.8x10^6 calls/s (-O2 optimization).
* 
* 
* Example calls:
* randn_notrig();        //returns normal deviate with mean=0.0, std. deviation=1.0
* randn_notrig(5.2,3.0);        //returns deviate with mean=5.2, std. deviation=3.0
* 
* 
* Dependencies - requires <cmath> for the sqrt(), sin(), and cos() calls, and a
* #defined value for PI.
*/
 
//#include <stdio.h>
#include <stdlib.h>
// #include <stdarg.h>
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

double randn_notrig( double mu, double sigma, double seed ) {
  
  static bool deviateAvailable=false; // flag
  static float storedDeviate; // deviate from previous calculation
  double polar, rsquared, var1, var2;
  
  srand ( seed );
  // choose pairs of uniformly distributed deviates, discarding
  // those that don't fall within the unit circle
  do {
    var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
    var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
    rsquared=var1*var1+var2*var2;
  } while ( rsquared>=1.0 || rsquared == 0.0);
  
  // calculate polar tranformation for each deviate
  polar=sqrt(-2.0*log(rsquared)/rsquared);
    
  // return second deviate
  return var2 * polar * sigma + mu;
}

int main( int argc, char * argv[] ) {

  double mu;
  double sigma;
  double output;
  double seed;

  if (argc == 1) {
    mu = 0;
    sigma = 1;
  } else if (argc == 4) {
    mu = atof( argv[1] );
    sigma = atof( argv[2] );
    seed = atof( argv[3] );
  } else {
    // std::cerr << "Usage: " << argv[0] << " mean sigma" << std::endl;
    return 1;
  }

  output = randn_notrig( mu, sigma, seed );
  cout << output << endl;
  
  return 0;
}
