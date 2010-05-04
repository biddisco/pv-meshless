/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "ftipsy.hpp"

bool Same(float a, float b, float eps = 1e-5) {
    if ( fabs(a-b) > eps ) return false;
    return true;
}


void Display( int n, TipsyDarkParticle &p1, TipsyDarkParticle &p2 ) {
    bool same = true;

    same &= Same(p1.pos[0],p2.pos[0], 1e-7);
    same &= Same(p1.pos[1],p2.pos[1], 1e-7);
    same &= Same(p1.pos[2],p2.pos[2], 1e-7);
    same &= Same(p1.vel[0],p2.vel[0]);
    same &= Same(p1.vel[1],p2.vel[1]);
    same &= Same(p1.vel[2],p2.vel[2]);
    same &= Same(p1.mass,  p2.mass,  0.0f);

    if ( ! same ) {
        std::cout
            << n << " "
            << fabs(p1.pos[0]-p2.pos[0]) << " " 
            << fabs(p1.pos[1]-p2.pos[1]) << " " 
            << fabs(p1.pos[2]-p2.pos[2]) << " " 
            << fabs(p1.vel[0]-p2.vel[0]) << " " 
            << fabs(p1.vel[1]-p2.vel[1]) << " " 
            << fabs(p1.vel[2]-p2.vel[2]) << " "
            << fabs(p1.mass-p2.mass) << std::endl;
    }
}
void Display( int n, TipsyGasParticle &p1, TipsyGasParticle &p2 ) {
    bool same = true;

    same &= Same(p1.pos[0],p2.pos[0], 1e-7);
    same &= Same(p1.pos[1],p2.pos[1], 1e-7);
    same &= Same(p1.pos[2],p2.pos[2], 1e-7);
    same &= Same(p1.vel[0],p2.vel[0]);
    same &= Same(p1.vel[1],p2.vel[1]);
    same &= Same(p1.vel[2],p2.vel[2]);
    same &= Same(p1.mass,  p2.mass,  0.0f);

    if ( ! same ) {
        std::cout 
            << n << " "
            << fabs(p1.pos[0]-p2.pos[0]) << " " 
            << fabs(p1.pos[1]-p2.pos[1]) << " " 
            << fabs(p1.pos[2]-p2.pos[2]) << " " 
            << fabs(p1.vel[0]-p2.vel[0]) << " " 
            << fabs(p1.vel[1]-p2.vel[1]) << " " 
            << fabs(p1.vel[2]-p2.vel[2]) << " "
            << fabs(p1.mass-p2.mass) << std::endl;
    }
}
void Display( int n, TipsyStarParticle &p1, TipsyStarParticle &p2 ) {
    bool same = true;

    same &= Same(p1.pos[0],p2.pos[0], 1e-7);
    same &= Same(p1.pos[1],p2.pos[1], 1e-7);
    same &= Same(p1.pos[2],p2.pos[2], 1e-7);
    same &= Same(p1.vel[0],p2.vel[0]);
    same &= Same(p1.vel[1],p2.vel[1]);
    same &= Same(p1.vel[2],p2.vel[2]);
    same &= Same(p1.mass,  p2.mass,  0.0f);

    if ( ! same ) {
        std::cout 
            << n << " "
            << (p1.pos[0]-p2.pos[0]) << " " 
            << (p1.pos[1]-p2.pos[1]) << " " 
            << (p1.pos[2]-p2.pos[2]) << " " 
            << fabs(p1.vel[0]-p2.vel[0]) << " " 
            << fabs(p1.vel[1]-p2.vel[1]) << " " 
            << fabs(p1.vel[2]-p2.vel[2]) << " "
            << fabs(p1.mass-p2.mass) << std::endl;
    }
}

int main( int argc, char *argv[] ) {
    ifTipsy in1, in2;

    TipsyHeader       h;
    TipsyGasParticle  g1, g2;
    TipsyDarkParticle d1, d2;
    TipsyStarParticle s1, s2;

    uint32_t i;

    if ( argc < 3 ) {
        fprintf( stderr, "Usage: %s <instd1> <instd2>\n", argv[0] );
        exit(1);
    }

    in1.open(argv[1]);
    if ( ! in1.is_open() ) {
        fprintf( stderr, "Unable to open Tipsy binary %s\n", argv[1] );
        exit(2);
    }
    in2.open(argv[2]);
    if ( ! in2.is_open() ) {
        fprintf( stderr, "Unable to open Tipsy binary %s\n", argv[2] );
        exit(2);
    }

    in1 >> h;
    in2 >> h;

    std::cout	<< "    dx         dy         dz     "
                << "    dVx        dVy        dVz      dMass" << std::endl;

    std::cout.flags(std::ios::right | std::ios::scientific | std::ios::showpos );
    std::cout.precision(3);

    for( i=0; i<h.h_nSph;  i++ ) { in1 >> g1; in2 >> g2; Display(i,g1,g2); }
    for( i=0; i<h.h_nDark; i++ ) { in1 >> d1; in2 >> d2; Display(i,d1,d2); }
    for( i=0; i<h.h_nStar; i++ ) { in1 >> s1; in2 >> s2; Display(i,s1,s2); }

    in1.close();
    in2.close();
}
