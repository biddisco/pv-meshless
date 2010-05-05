/*=========================================================================

  Project: meshless plugin

  Copyright (c) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing this
  copyright notice appears on all copies of source code and an
  acknowledgment appears with any substantial usage of the code.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

  This code is derived from an earlier work and is distributed
  with permission from, and thanks to

=========================================================================*/
/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkTipsyReader.cxx,v $
=========================================================================*/

/*=========================================================================
                                                                                
Copyright (c) 2007, Los Alamos National Security, LLC

All rights reserved.

Copyright 2007. Los Alamos National Security, LLC. 
This software was produced under U.S. Government contract DE-AC52-06NA25396 
for Los Alamos National Laboratory (LANL), which is operated by 
Los Alamos National Security, LLC for the U.S. Department of Energy. 
The U.S. Government has rights to use, reproduce, and distribute this software. 
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
If software is modified to produce derivative works, such modified software 
should be clearly marked, so as not to confuse it with the version available 
from LANL.
 
Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
-   Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer. 
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution. 
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software 
    without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
                                                                                
=========================================================================*/

// .NAME TipsyParticleDistribute - distribute particles to processors
//
// .SECTION Description
// TipsyParticleDistribute takes a series of data files containing RECORD style
// .cosmo data or Gadget style BLOCK data
// along with parameters defining the box size for the data and for
// determining halos within the particle data.  It distributes the data
// across processors including a healthy dead zone of particles belonging
// to neighbor processors.  By definition all halos can be determined
// completely for any processor because of this dead zone.  The serial
// halo finder is called on each processor.
//

#ifndef TipsyParticleDistribute_h
#define TipsyParticleDistribute_h

#include "ParticleDistribute.h"
#include "tipsylib/ftipsy.hpp" // functions take tipsy particle objects

class TipsyParticleDistribute : public ParticleDistribute {
public:
   TipsyParticleDistribute();
  ~TipsyParticleDistribute();

  // Read particle files per processor and share round robin with others
  // extracting only the alive particles
  virtual void readParticlesRoundRobin(int reserveQ=0);
  virtual void partitionInputFiles();

  // Read one particle file per processor with alive particles 
  // and correct topology
  virtual void readParticlesOneToOne(int reserveQ=0);

  // Get particle counts for allocating buffers
  virtual void findFileParticleCount();

  void readFromTipsyFile(
        ifTipsy* inStream,     // Stream to read from
        int firstParticle,      // First particle index to read in this chunk
        int numberOfParticles,  // Number of particles to read in this chunk
        POSVEL_T* fblock,       // Buffer for read in data
        ID_T* iblock,           // Buffer for read in data
        Message* message);      // Message buffer for distribution

  // One to one version of read is simpler with no MPI buffering
  // Not implemented yet for Tipsy files
  void readFromTipsyFile();

  // When set, the bounding box is recentred from the origin to {0.5,0.5,0.5}
  int RecentreBoundingBox;

private:

};

#endif
