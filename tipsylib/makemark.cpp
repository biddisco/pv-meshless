#include <iostream>
#include <stdint.h>
#include <getopt.h>
#include <math.h>
#include "ftipsy.hpp"
#include "grid.hpp"

#define OPT_GRID 'g'


void outIf( Grid &grid, TipsyBaseParticle &b, uint32_t iIndex )
{
    if ( grid.inGrid( b.pos[0], b.pos[1], b.pos[2] ) ) {
	std::cout << iIndex << std::endl;
    }
}

int main( int argc, char *argv[] ) {
    ifTipsy in;          // The input file

    const char *gridName = 0;
    const char *tipsyName = 0;

    TipsyHeader       h; // The header structure
    TipsyGasParticle  g; // A gas particle
    TipsyDarkParticle d; // A dark particle
    TipsyStarParticle s; // A star particle
    uint32_t i, iIndex;

    Grid grid;

    //! Parse command line
    for(;;) {
        int c, option_index=0;

        static struct option long_options[] = {
            { "grid",        1, 0, OPT_GRID },
            { 0,             0, 0, 0 }
        };

        c = getopt_long( argc, argv, "g:",
                         long_options, &option_index );
        if ( c == -1 ) break;
        switch(c) {
        case OPT_GRID:
	    gridName = optarg;
            break;
	}
    }

    if ( optind < argc ) {
        tipsyName = argv[optind++];
    }
    else {
        std::cerr << "Missing tipsy input file" << std::endl;
        exit(2);
    }

    in.open(tipsyName,"standard");
    if ( ! in.is_open() ) {
	std::cerr <<"Unable to open Tipsy binary " << tipsyName << std::endl;
        exit(2);
    }

    if ( gridName ) {
	grid.loadGrid( gridName );
    }
    else {
	grid.loadGrid( std::cin );
    }

    // Read the header from the input and write it to the output.
    in >> h;

    std::cout << h.h_nDark << " " << h.h_nSph << " " << h.h_nStar << std::endl;


    // Read every particle and write it to the output file.
    iIndex = 1;
    for( i=0; i<h.h_nSph;  i++,iIndex++ ) { in >> g; outIf(grid,g,iIndex); }
    for( i=0; i<h.h_nDark; i++,iIndex++ ) { in >> d; outIf(grid,d,iIndex); }
    for( i=0; i<h.h_nStar; i++,iIndex++ ) { in >> s; outIf(grid,s,iIndex); }

    // Close the file.
    in.close();
}
