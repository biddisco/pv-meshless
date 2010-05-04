#include <iostream>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>
#include "grid.hpp"

#define OPT_GRID_BOUND  'b'
#define OPT_VERBOSE     'v'
#define OPT_GRID_FILL   'f'

int main( int argc, char *argv[] ) {
    int32_t GridBound = 0;
    int32_t GridFill  = 0;
    int32_t Verbose = 0;
    Grid gridA, gridB;//, gridC;

    //! Parse command line
    for(;;) {
        int c, option_index=0;

        static struct option long_options[] = {
            { "boundary",    1, 0, OPT_GRID_BOUND },
            { "fill",        1, 0, OPT_GRID_FILL },
            { 0,             0, 0, 0 }
        };

        c = getopt_long( argc, argv, "b:vf:",
                         long_options, &option_index );
        if ( c == -1 ) break;
        switch(c) {
        case OPT_GRID_BOUND:
	    assert( optarg != 0 );
	    GridBound = atoi(optarg);
            break;
	case OPT_VERBOSE:
	    Verbose++;
	    break;
        case OPT_GRID_FILL:
	    assert( optarg != 0 );
	    GridFill = atoi(optarg);
            break;
	default:
	    assert(0);
	    break;
	}
    }
    if ( GridFill==0 && GridBound==0 ) {
	std::cerr << "Specify a fill or a boundary" << std::endl;
	exit(11);
    }


    if ( !gridA.loadGrid( std::cin ) ) {
	std::cerr << "Unable to open GRID input" << std::endl;
	exit(11);
    }

    if ( Verbose )
	std::clog << "Before: " << gridA.getGridCount() << " cells marked"
		<< std::endl;
    if ( GridBound )
	gridB.addBorder(gridA,GridBound);
    else {
	// Why twice?  For the famous cross problem
	//gridC.fill(gridA,GridFill);
	//gridB.fill(gridC,2);
	gridB.fill(gridA,GridFill);

	//gridC.fill(gridB,GridFill);
	//gridB.fill(gridC,GridFill);
    }

    if ( Verbose )
	std::clog << "After:  " << gridB.getGridCount() << " cells marked"
		<< std::endl;

    if ( !gridB.saveGrid( std::cout ) ) {
	std::cerr << "Unable to create GRID output" << std::endl;
	exit(11);
    }
}
