/**
 *  @file
 *  @brief Limit Tracker.
 *  @author Doug Potter
 */

/*
 *  The limitTracker class will track box limits (max/min x,y,z) and takes
 *  into account the periodic boundaries.  If the size of the structure
 *  being tracked exceeds half the box size, then the resulting limits will
 *  be correct.  If not, the center of the limits is always inside the real
 *  box, but one of the limits may be outside the physical box.
 */

#ifndef TIPLIM_H
#define TIPLIM_H

namespace Tipsy {

    class limitTracker;
    class oneLimit;

    class oneLimit {
    protected:
	float minBox, maxBox;   //!< Physical box size
	float centerBox;

	float minVal, maxVal;   //!< Limits of the structure
	float minAlt, maxAlt;   //!< Alternate limits

    public:
	void setLimits( float left, float right ) {
	    maxVal = maxAlt = left;
	    minVal = minAlt = right;
	    minBox = left;
	    maxBox = right;
	    centerBox = (minBox+maxBox) * 0.5;
	}

	oneLimit( float left = -0.5, float right = 0.5 ) {
	    setLimits(left,right);
	}

	float getMin() { return minVal; }
	float getMax() { return maxVal; }
	float getSize() { return maxVal-minVal; }
	float getCenter() { return (minVal+maxVal) * 0.5; }

	// Use with caution
	void  setVals(float v1, float v2 ) { minVal = v1; maxVal = v2; }

	// Change the minimum and maximum values to reflect reality
	void Adjust( float val ) {
	    // Allow periodic inputs (fix somehow perhaps)
	    //assert( val >= minBox && val < maxBox );

	    // Just stuff in the value if this is the first one
	    if ( minVal > maxVal ) {
		minVal = maxVal = val;
	    }
	    else {
		// Now adjust the limits
		if ( val < minVal ) minVal = val;
		if ( val > maxVal ) maxVal = val;
	    }
	    if ( val >= centerBox ) {
		if ( val < minAlt ) minAlt = val;
	    }
	    else {
		if ( val > maxAlt ) maxAlt = val;
	    }
	}

	// Fix up the values internally
	void Normalize( bool bAlternate=true ) {

	    // Okay, there was nothing to track -- make it sensible
	    if ( minVal > maxVal )
		minVal = maxVal = 0.0;

	    // Fix up the coordinates properly
	    else {
		float sizeBox = (maxBox-minBox);
		float halfBox = sizeBox * 0.5;

		// If the alternate limits are best, then use them
		if ( maxAlt + sizeBox - minAlt < halfBox
		     && maxVal - minVal > halfBox && bAlternate ) {
		    float center;

		    minVal = minAlt;
		    maxVal = maxAlt + sizeBox;

		    // Make sure the center is inside the real box
		    center = (minVal + maxVal) * 0.5;
		    if ( center > maxBox ) {
			minVal -= sizeBox;
			maxVal -= sizeBox;
		    }
		}
	    }
	}

	// Returns the value adjusted for the periodic boundary
	float boxValue( float val ) {
	    float sizeBox = (maxBox-minBox);
	    if ( val > maxVal ) val -= sizeBox;
	    else if ( val < minVal ) val += sizeBox;
	    return val;
	}

	// Returns true if the value is inside the range.  Adjusts for
	// periodic boundary.
	bool inRange( float val ) {
	    val = boxValue( val );
	    return val >= minVal && val <= maxVal;
	}



    };

    class limitTracker {
    public:
	oneLimit X, Y, Z;

	limitTracker( float left = -0.5, float right = 0.5 )
	    : X(left,right), Y(left,right), Z(left,right)
	    { }

	void setLimits( float left, float right ) {
	    X.setLimits(left,right);
	    Y.setLimits(left,right);
	    Z.setLimits(left,right);
	}

	void Adjust( float x, float y, float z ) {
	    X.Adjust(x);
	    Y.Adjust(y);
	    Z.Adjust(z);
	}

	void Adjust( const float *pos ) {
	    Adjust( pos[0], pos[1], pos[2] );
	}

	void Normalize(bool bAlternate=true ) {
	    X.Normalize(bAlternate);
	    Y.Normalize(bAlternate);
	    Z.Normalize(bAlternate);
	}

	// Returns true if the specified coordinates are inside the box.
	bool inBox( float x, float y, float z ) {
	    return X.inRange(x) && Y.inRange(y) && Z.inRange(z);
	}

	bool inBox( const float *pos ) {
	    return inBox( pos[0], pos[1], pos[2] );
	}



    };

}
#endif
