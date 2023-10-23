#include <iostream>
#include <limits>

typedef double EulerAngleType;
typedef double CoordinateType;

void SPAEulerAngles( EulerAngleType psi, // particle orientation
                     EulerAngleType theta,
                     EulerAngleType phi,
                     EulerAngleType tilt_angle, // micrograph angles
                     EulerAngleType tilt_axis_angle,
                     EulerAngleType normal[3], // spike normal
                     EulerAngleType matrix[16], // spike refinement
                     EulerAngleType frealign[3],
                     CoordinateType shift[2] );

using namespace std;

int main( int argc, char ** argv ) {
  if (argc != 25) {
    cout << "USAGE: " << argv[0]
        << " tilt_angle axis_rotation axis_correction box_correction[2] normal[0] normal[1] normal[2] m[0-15] \n"
        << endl;
    return 0;
  }
  
  // tilt angle
  EulerAngleType tilt_angle = atof( argv[1] );

  // tilt axis rotation
  EulerAngleType tilt_axis_angle = atof( argv[2] );
  EulerAngleType axis_correction = atof( argv[3] );
  EulerAngleType tilt_X_correction = atof( argv[4] );
  EulerAngleType tilt_Y_correction = atof( argv[5] );

  // spike normal
  EulerAngleType normal[3] = { atof( argv[6] ), atof( argv[7] ), atof( argv[8] ) };

  EulerAngleType matrix[16];
  matrix[0] = atof(argv[9]);
  matrix[1] = atof(argv[10]);
  matrix[2] = atof(argv[11]);
  matrix[3] = atof(argv[12]);
  matrix[4] = atof(argv[13]);
  matrix[5] = atof(argv[14]);
  matrix[6] = atof(argv[15]);
  matrix[7] = atof(argv[16]);
  matrix[8] = atof(argv[17]);
  matrix[9] = atof(argv[18]);
  matrix[10] = atof(argv[19]);
  matrix[11] = atof(argv[20]);
  matrix[12] = atof(argv[21]);
  matrix[13] = atof(argv[22]);
  matrix[14] = atof(argv[23]);
  matrix[15] = atof(argv[24]);
  
  EulerAngleType frealign[3];
  CoordinateType shift[2];
  
  SPAEulerAngles( 0, 0, 0, tilt_angle, tilt_axis_angle, normal, matrix, frealign, shift );

  printf( "%.2f %.2f %.2f %.2f %.2f\n", frealign[0], frealign[1], frealign[2], shift[0], shift[1] );

  return 0;
}
