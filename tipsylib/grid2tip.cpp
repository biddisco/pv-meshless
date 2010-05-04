#include <iostream>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>
#include "grid.hpp"
#include "ftipsy.hpp"


int main( int argc, char *argv[] ) {
    bool bError = false;
    Grid gridA;
    int i, j, k;
    float x,y,z;
    int S, E;
    TipsyDarkParticle D;
    TipsyHeader H;

    //! Parse command line
    for(;;) {
        int c, option_index=0;

        static struct option long_options[] = {
            { 0,             0, 0, 0 }
        };

        c = getopt_long( argc, argv, "",
                         long_options, &option_index );
        if ( c == -1 ) break;
        switch(c) {
	default:
	    bError = true;
	    break;
	}
    }

    if ( !gridA.loadGrid( std::cin ) ) {
	std::cerr << "Unable to open GRID input" << std::endl;
	exit(11);
    }

    ofTipsy out("test.std");

    S = gridA.getGridSize() / 2;
    H.h_nDark = H.h_nBodies = S*S*S;
    E = S + S/2 + 488;
    S = S/2 + 488;

    H.h_time = 1.0;
    H.h_nDims = 3;
    H.h_nSph = H.h_nStar = 0;
    out << H;

    D.vel[0] = D.vel[1] = D.vel[2] = D.phi = D.density = D.eps = 0.0;

//"-2.189652e-01,-4.229806e-01,-3.679137e-01"
//    -107,-206,-180

    for( i=S-107; i<E-107; i++ ) {
	//std::cerr << i << std::endl;
	D.pos[0] = static_cast<float>(i) / gridA.getGridSize() - 0.5 + 2.189652e-01;
	if ( D.pos[0] >= 0.5 ) D.pos[0] -= 1.0;
	if ( D.pos[0] < -0.5 ) D.pos[0] += 1.0;
	for( j=S-206; j<E-206; j++ ) {
	    D.pos[1] = static_cast<float>(j) / gridA.getGridSize() - 0.5 + 4.229806e-01;
	    if ( D.pos[1] >= 0.5 ) D.pos[1] -= 1.0;
	    if ( D.pos[1] < -0.5 ) D.pos[1] += 1.0;

	    for( k=S-180; k<E-180; k++ ) {
		D.pos[2] = static_cast<float>(k) / gridA.getGridSize() - 0.5 + 3.679137e-01;
		if ( D.pos[2] >= 0.5 ) D.pos[2] -= 1.0;
		if ( D.pos[2] < -0.5 ) D.pos[2] += 1.0;

		D.mass = gridA.getGrid(i%488,j%488,k%488) ? 1.0 : 0.0;
		out << D;
	    }
	}
    }
    



}
