/*  This example program will convert a native format tipsy file
 *  into a standard format file.
 */
#include <stdlib.h>
#include "ftipsy.hpp"

int main( int argc, char *argv[] ) {
    ifTipsy in;          // The input file
    ofTipsy out;         // The output file

    TipsyHeader       h; // The header structure
    TipsyGasParticle  g; // A gas particle
    TipsyDarkParticle d; // A dark particle
    TipsyStarParticle s; // A star particle
    uint32_t i;

    // Make sure we have two parameters, a native and a standard file.
    if ( argc != 3 ) {
        fprintf( stderr, "Usage: %s <native> <standard>\n", argv[0] );
        exit(1);
    }

    // Open the native file and abort if there is an error.
    in.open(argv[1],"native");
    if ( ! in.is_open() ) {
        fprintf( stderr, "Unable to open Tipsy binary %s\n", argv[1] );
        exit(2);
    }

    // Open the output file.
    out.open(argv[2],"standard");
    if ( ! out.is_open() ) {
        fprintf( stderr, "Unable to create Native binary %s\n", argv[2] );
        exit(2);
    }

    // Read the header from the input and write it to the output.
    in >> h; out << h;

    // Read every particle and write it to the output file.
    for( i=0; i<h.h_nSph;  i++ ) { in >> g; out << g; }
    for( i=0; i<h.h_nDark; i++ ) { in >> d; out << d; }
    for( i=0; i<h.h_nStar; i++ ) { in >> s; out << s; }

    // Close the files.
    out.close();
    in.close();
}
