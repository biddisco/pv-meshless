/*  This  program takes a Tipsy file and bins particles
 */

#include <iostream>
#include <math.h>
#include <getopt.h>
#include "ftipsy.hpp"
#include "binner.hpp"

#define OPT_HELP   'h'
#define OPT_RADIUS 'r'
#define OPT_CENTER 'c'

static void Usage( const char *name ) {
    std::clog << "Usage: " << name
              << " [-h] -r <radius> [-c <center>] <input> <output>"
              << std::endl
              << "  -h,--help\t\tShow help (this text)" << std::endl
              << "  -r,--radius=0.1\tSet radius at which to bin" << std::endl
	      << "  -c,--center=x,y,z\tSet the new center" << std::endl
	      << std::endl;
}

static void BadCenter(void) {
    std::clog << "Specify center as --center=x,y,z" << std::endl;
    exit(1);
}

static void BadRadius(void) {
    std::clog << "Specify radius as --radius=R or --radius=Rx,Ry,Rz"
	      << std::endl;
    exit(1);
}

static inline float convertFloat( const char *p )
{
    if ( p == NULL ) BadCenter();
    return atof(p);
}


int main( int argc, char *argv[] ) {
    bool bHelp  = false;    //!< Help was requested on the command line
    bool bError = false;    //!< An error occurred on the command line
    float Rx = 0.0;
    float Ry = 0.0;
    float Rz = 0.0;
    float cx=0.0;
    float cy=0.0;
    float cz=0.0;
    int nLevels = 2;
    int iB, iN, iC;
    double dMassIn, dMassOut;
    std::string inName;
    std::string outName;

    ifTipsy in;          // The input file
    ofTipsy out;         // The output file

    //Tipsy::SphericalBinner binner;
    Tipsy::RadialBinner binner;

    TipsyHeader       h; // The header structure
    TipsyGasParticle  g; // A gas particle
    TipsyDarkParticle d; // A dark particle
    TipsyStarParticle s; // A star particle
    uint32_t i;

    std::clog << "WARNING: This tool is experimental" << std::endl;

    //! Parse command line
    for(;;) {
        int c, option_index=0;
	char *p;

        static struct option long_options[] = {
            { "help",        0, 0, OPT_HELP },
            { "radius",      1, 0, OPT_RADIUS },
            { "center",      1, 0, OPT_CENTER },
            { "centre",      1, 0, OPT_CENTER },
        };

        c = getopt_long( argc, argv, "hc:r:",
                         long_options, &option_index );
        if ( c == -1 ) break;
        switch(c) {
        case OPT_HELP:
            bHelp = true;
            break;
	case OPT_RADIUS:
	    assert( optarg != NULL );
	    if ( Rx > 0 ) {
		std::clog << "Specify radius only once" << std::endl;
		bError = true;
	    }
	    p = strtok(optarg,",");
	    assert( p != NULL);
	    Rx = atof(p);
	    p = strtok(NULL,",");
	    if ( p == NULL ) {
		Ry = Rz = Rx;
	    }
	    else {
		Ry = atof(p);
		p = strtok(NULL,",");
		if ( p == NULL ) BadRadius();
		Rz = atof(p);
		p = strtok(NULL,",");
		if ( p != NULL ) BadRadius();
	    }

	    if ( Rx <= 0 || Ry <= 0 || Rz <= 0 ) {
		std::clog << "Radius must be positive (not zero)" << std::endl;
		bError = true;
	    }
	    break;
	case OPT_CENTER:
	    assert( optarg != NULL );
	    p = strtok(optarg,",");
	    cx = convertFloat(p);
	    p = strtok(NULL,",");
	    cy = convertFloat(p);
	    p = strtok(NULL,",");
	    cz = convertFloat(p);
	    p = strtok(NULL,",");
	    if ( p != NULL ) BadCenter();
	    break;
        default:
            bError = true;
            break;
        }
    }

    if ( bHelp || bError ) {
        Usage(argv[0]);
        exit(1);
    }

    if ( Rx <= 0 ) {
	std::clog << "Specify the binning radius with --radius" << std::endl;
	exit(1);
    }

    if ( optind < argc ) {
	inName = argv[optind++];
    }
    else {
	std::cerr << "Missing input file" << std::endl;
	exit(2);
    }

    if ( optind < argc ) {
	outName = argv[optind++];
    }
    else {
	std::cerr << "Missing output file" << std::endl;
	exit(2);
    }

    // Open the input file for binning
    in.open(inName.c_str());
    if ( ! in.is_open() ) {
        fprintf( stderr, "Unable to open Tipsy binary %s\n", argv[1] );
        exit(2);
    }

    // Open the output file.
    out.open(outName.c_str());
    if ( ! out.is_open() ) {
        fprintf( stderr, "Unable to create Tipsy binary %s\n", argv[2] );
        exit(2);
    }

    iB = iN = iC = 0;
    dMassIn = dMassOut = 0.0;

    // Read the header from the input and write it to the output.
    // This won't be correct, but we will fix it up at the end.
    in >> h; out << h;


    // Calculate the binning
    float cr = powf(h.h_nBodies,1.0/3.0);
    nLevels = 0;
    do { nLevels++; } while( (cr*=0.5) >= 100 );
    std::clog << "Grid: " << (int)cr << ", levels: " << nLevels << std::endl;


    //binner.setTheta(0.01);
    //binner.setBox(cx,cy,cz);
    //binner.setRadius(Radius);
    binner.setBox(cx,cy,cz);
    //binner.setBinning( (int)cr, nLevels, Rx, Ry, Rz );
    binner.setBinning( nLevels, (int)cr, (int)cr, (int)cr, Rx, Ry, Rz );

    // Read every particle and write it to the output file.
    for( i=0; i<h.h_nSph;  i++ ) {
	in >> g; out << g;
    }
    for( i=0; i<h.h_nDark; i++ ) {
	in >> d;
	dMassIn += d.mass;
	if ( !binner.Bin(d) ) {
	    iN++;
	    out << d;
	    dMassOut += d.mass;
	}
	else {
	    iB++;
	}
    }
    for( i=0; i<h.h_nStar; i++ ) {
	in >> s;
	out << s;
    }


    // Write the binned, dark particles
    while( binner.getParticle(d) ) {
	iC++;
	out << d;
	dMassOut += d.mass;
    }

    // Update the header
    h.h_nDark = iN + iC;
    h.h_nBodies = h.h_nDark;

    out.seekp( tipsypos( tipsypos::header, 0 ) );
    out << h;

    // Close the files.
    out.close();
    in.close();

    std::clog << "Binned " << iB << " particles into " << iC << " and wrote " << iN << std::endl;

    std::clog << "Mass check: "
	      << " in: " << dMassIn
	      << " out: " << dMassOut
	      << std::endl;
}
