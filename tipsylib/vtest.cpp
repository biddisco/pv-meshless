/*  This example program tests the vTipsy class.
 */

#include "vtipsy.hpp"

int main( int argc, char *argv[] ) {
    ifTipsy in;				// The input file
    ofTipsy out;			// The output file

    vTipsy v1;				// Our first in-memory array
    vTipsy v2;				// The second in-memory array

    vTipsy::gas_type::iterator  gi;	// Gas particle iterator
    vTipsy::dark_type::iterator di;	// Dark particle iterator
    vTipsy::star_type::iterator si;	// Star particle iterator

    // Make sure we have two parameters, two standard tipsy files.
    if ( argc != 3 ) {
        fprintf( stderr, "Usage: %s <infile> <outfile>\n", argv[0] );
        exit(1);
    }

    // Open the input file and abort if there is an error.
    in.open(argv[1]);
    if ( ! in.is_open() ) {
        fprintf( stderr, "Unable to open Tipsy binary %s\n", argv[1] );
        exit(2);
    }

    // Open the output file.
    out.open(argv[2]);
    if ( ! out.is_open() ) {
        fprintf( stderr, "Unable to create Tipsy binary %s\n", argv[2] );
        exit(2);
    }

    // Read the input file directly into the first array
    v1 << in;

    // We can close the input file now, because we are done with it.
    in.close();

    // Set the expansion factor (simulation time).  Copy it from v1 to v2.
    v2.expFactor = v1.expFactor;

    // Now copy each particle type from the first to the second array.
    for( gi=v1.gas.begin(); gi!=v1.gas.end(); gi++ ) {
#ifdef REVERSE_DARK_VELOCITY
	// An example of manipulating a particle.
	gi->vel[0] = -gi->vel[0];
	gi->vel[1] = -gi->vel[1];
	gi->vel[2] = -gi->vel[2];
#endif
	v2 << *gi;
    }

    for( di=v1.dark.begin(); di!=v1.dark.end();  di++ ) {
	v2 << *di;
    }

    for( si=v1.star.begin(); si!=v1.star.end();  si++ ) {
	v2 << *si;
    }

    // Write the second array to the output file
    v2 >> out;

    // Close the output file.
    out.close();
}
