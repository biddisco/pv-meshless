/*  This example program will convert a native format tipsy file
 *  into a standard format file.
 */

#include "ftipsy.hpp"

static inline float adjust(float v,float c) {
    v -= c;
    if ( v < -0.5 ) v += 1.0;
    if ( v >= 0.5 ) v -= 1.0;
    return v;
    }

static inline void recenter(TipsyBaseParticle &p, float *c) {
    p.pos[0] = adjust(p.pos[0],c[0]);
    p.pos[1] = adjust(p.pos[1],c[1]);
    p.pos[2] = adjust(p.pos[2],c[2]);
    }

int main( int argc, char *argv[] ) {
    ifTipsy in;          // The input file
    ofTipsy out;         // The output file

    TipsyHeader       h; // The header structure
    TipsyGasParticle  g; // A gas particle
    TipsyDarkParticle d; // A dark particle
    TipsyStarParticle s; // A star particle
    uint32_t i;
    float c[3];


    // Make sure we have two parameters, a native and a standard file.
    if ( argc != 6 ) {
        fprintf( stderr, "Usage: %s <in >out cx cy cz\n", argv[0] );
        exit(1);
    }

    in.open(argv[1],"standard");
    if ( ! in.is_open() ) {
        fprintf( stderr, "Unable to open Tipsy binary %s\n", argv[1] );
        exit(2);
    }

    out.open(argv[2],"standard");
    if ( ! out.is_open() ) {
        fprintf( stderr, "Unable to create Native binary %s\n", argv[2] );
        exit(2);
    }

    c[0] = atof(argv[3]);
    c[1] = atof(argv[4]);
    c[2] = atof(argv[5]);

    // Read the header from the input and write it to the output.
    in >> h; out << h;

    // Read every particle and write it to the output file.
    for( i=0; i<h.h_nSph;  i++ ) { in >> g; recenter(g,c); out << g; }
    for( i=0; i<h.h_nDark; i++ ) { in >> d; recenter(d,c); out << d; }
    for( i=0; i<h.h_nStar; i++ ) { in >> s; recenter(s,c); out << s; }

    // Close the files.
    out.close();
    in.close();
}
