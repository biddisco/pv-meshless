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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <sys/types.h>

#ifdef _WIN32
#include "winDirent.h"
#else
#include <dirent.h>
#endif

#include "Partition.h"
#include "TipsyParticleDistribute.h"

#ifdef USE_VTK_COSMO
#include "vtkStdString.h"
#include "vtkSetGet.h"
#endif

#include "tipsylib/ftipsy.hpp" // functions take tipsy particle objects

using namespace std;

const int  MAX_TIPSY_READ = 1000000;
const bool ENFORCE_TIPSY_MAX_READ = true;

/////////////////////////////////////////////////////////////////////////
//
// Particle data space is partitioned for the number of processors
// which currently is a factor of two but is easily extended.  Particles
// are read in from files where each processor reads one file into a buffer,
// extracts the particles which really belong on the processor (ALIVE) and
// those in a buffer region around the edge (DEAD).  The buffer is then
// passed round robin to every other processor so that all particles are
// examined by all processors.  All dead particles are tagged with the
// neighbor zone (26 neighbors in 3D) so that later halos can be associated
// with zones.
//
/////////////////////////////////////////////////////////////////////////

TipsyParticleDistribute::TipsyParticleDistribute()
{
  this->RecentreBoundingBox = 1;
}

TipsyParticleDistribute::~TipsyParticleDistribute()
{
}


/////////////////////////////////////////////////////////////////////////
//
// Set box sizes for determining if a particle is in the alive or dead
// region of this processor.  Data space is a DIMENSION torus.
//
/////////////////////////////////////////////////////////////////////////


/*
void TipsyParticleDistribute::setParticles(vector<POSVEL_T>* xLoc,
                                      vector<POSVEL_T>* yLoc,
                                      vector<POSVEL_T>* zLoc,
                                      vector<POSVEL_T>* xVel,
                                      vector<POSVEL_T>* yVel,
                                      vector<POSVEL_T>* zVel,
                                      vector<POSVEL_T>* mass,
                                      vector<ID_T>* id)
{
  this->xx = xLoc;
  this->yy = yLoc;
  this->zz = zLoc;
  this->vx = xVel;
  this->vy = yVel;
  this->vz = zVel;
  this->ms = mass;
  this->tag = id;
}
*/

/////////////////////////////////////////////////////////////////////////
//
// Each processor reads 0 or more files, a buffer at a time, and shares
// the particles by passing the buffer round robin to every other processor
//
/////////////////////////////////////////////////////////////////////////

void TipsyParticleDistribute::readParticlesRoundRobin(int reserveQ)
{
  // Find how many input files there are and deal them between the processors
  // Calculates the max number of files per processor and max number of
  // particles per file so that buffering can be done
  // For round robin sharing determine where to send and receive buffers from
  this->partitionInputFiles();

  // Compute the total number of particles in the problem
  // Compute the maximum number of particles in any one file to set buffer size
  // The tipsy file is opened in here
  this->findFileParticleCount();

  if (!tipsyInfile.is_open()) {
    cout << " Error in file open" << std::endl;
    cout.flush();
  }

  // MPI buffer size might limit the number of particles read from a file
  // and passed round robin
  // Largest file will have a number of buffer chunks to send if it is too large
  // Every processor must send that number of chunks even if its own file
  // does not have that much information

  if (ENFORCE_TIPSY_MAX_READ==true && this->LocalCount>MAX_TIPSY_READ) {
    this->maxRead = MAX_TIPSY_READ;
    this->maxReadsPerFile = (this->LocalCount / this->maxRead) + 1;
    cout << "Rank " << this->myProc << " maxReadsPerFile " << this->maxReadsPerFile << endl;
    cout.flush();
  } else {
    this->maxRead = this->maxParticles;
    this->maxReadsPerFile = 1;
  }

  // Allocate space to hold buffer information for reading of files
  // Mass is constant use that float to store the tag
  // Number of particles is the first integer in the buffer
  int bufferSize = (1 * sizeof(int)) +          // number of particles
                   (this->maxRead *
                   ((TIPSY_SIZE * sizeof(POSVEL_T)) + (COSMO_INT * sizeof(ID_T))));

  Message* message1 = new Message(bufferSize);
  Message* message2 = new Message(bufferSize);

  // Allocate space for the data read from the file
  POSVEL_T *fBlock = 0;
  POSVEL_T *lBlock = 0;
  POSVEL_T *vBlock = 0;
  ID_T* iBlock = 0;

  // CUSTOM/TIPSY format reads one particle at a time
  fBlock = new POSVEL_T[TIPSY_SIZE];
  iBlock = new ID_T[TIPSY_INT];

  // Reserve particle storage to minimize reallocation
  int reserveSize = (unsigned long) (this->maxFiles * this->maxParticles * DEAD_FACTOR);

  // If multiple processors are reading the same file we can reduce size
  reserveSize /= this->processorsPerFile;

  if(1 || reserveQ) {
#ifndef DEBUG_MESSAGES
    cout << "readParticlesRoundRobin reserving vectors " << reserveSize << endl;
    cout.flush();
#endif
    this->xx->reserve(reserveSize);
    this->yy->reserve(reserveSize);
    this->zz->reserve(reserveSize);
    this->vx->reserve(reserveSize);
    this->vy->reserve(reserveSize);
    this->vz->reserve(reserveSize);
    this->ms->reserve(reserveSize);
    this->tag->reserve(reserveSize);
  }

  // Running total and index into particle data on this processor
  this->particleCount = 0;

  // Using the input files assigned to this processor, read the input
  // and push round robin to every other processor
  // this->maxFiles is the maximum number to read on any processor
  // Some processors may have no files to read but must still participate 
  // in the round robin distribution

  int firstParticle = this->LocalStartIndex;
  int numberOfParticles = 0;
  int remainingParticles = 0;

#ifndef DEBUG_MESSAGES
  cout << "Rank " << this->myProc << " open file " << inFiles[0]
       << " with " << this->fileParticles[0] << " particles" << endl;
    cout.flush();
#endif

  // Number of particles read at one time depends on MPI buffer size
  numberOfParticles = this->fileParticles[0];
  if (numberOfParticles > this->maxRead)
    numberOfParticles = this->maxRead;

  // If a file is too large to be passed as an MPI message divide it up
  remainingParticles = this->fileParticles[0];

  for (int piece = 0; piece < this->maxReadsPerFile; piece++) {

    // Reset each MPI message for each file read
    message1->reset();
    message2->reset();

    // Processor has a file to read and share via round robin with others
    readFromTipsyFile(&tipsyInfile, firstParticle, numberOfParticles, fBlock, iBlock, message1);
    firstParticle      += numberOfParticles;
    remainingParticles -= numberOfParticles;
    if (remainingParticles <= 0) numberOfParticles = 0;
    else if (remainingParticles < numberOfParticles) numberOfParticles = remainingParticles;

    // Particles belonging to this processor are put in vectors
    distributeParticles(message1, message2);
  }

  // Can delete the read buffers as soon as last file is read because
  // information has been transferred into the double buffer1
  delete [] fBlock;
  delete [] iBlock;

  // After all particles have been distributed to vectors the double
  // buffers can be deleted
  delete message1;
  delete message2;

  tipsyInfile.close();

  // Count the particles across processors
  long totalAliveParticles = 0;
#ifdef USE_SERIAL_COSMO
  totalAliveParticles = this->numberOfAliveParticles;
#else
  MPI_Allreduce((void*) &this->numberOfAliveParticles, 
                (void*) &totalAliveParticles, 
                1, MPI_LONG, MPI_SUM, Partition::getComm());
#endif

#ifndef DEBUG_MESSAGES
  cout << "Rank " << setw(3) << this->myProc 
       << " #alive = " << this->numberOfAliveParticles << endl;
    cout.flush();
 
  if (this->myProc == MASTER) {
    cout << "TotalAliveParticles " << totalAliveParticles << endl;
    cout.flush();
  }
#endif
}

/////////////////////////////////////////////////////////////////////////
//
// There is only one file with tipsy, so 
// how many input files there are.  Parcel those files between all the
// processors which will be responsible for actually reading 0 or more.
//
/////////////////////////////////////////////////////////////////////////

void TipsyParticleDistribute::partitionInputFiles()
{
  // All processors are reading the same file
  this->inFiles.push_back(this->baseFile);
  this->numberOfFiles = 1;

  // Assign the round robin circle (last circle is bigger than others)
  this->processorsPerFile       = this->numProc / this->numberOfFiles;

  // Number of round robin sends to share all the files
  this->numberOfFileSends = this->numProc - 1;
  this->maxFileSends      = this->numberOfFileSends;

  // Where is the file sent, and where is it received
  if (this->myProc == this->numProc - 1)
    this->nextProc = 0;
  else
    this->nextProc = this->myProc + 1;

  if (this->myProc == 0)
    this->prevProc = this->numProc - 1;
  else
    this->prevProc = this->myProc - 1;

}

/////////////////////////////////////////////////////////////////////////
//
// Open the tipsy file and find how many particles there are in total
// and how many this process is reading.
// @TODO : this could be done on just one process since we all read the same file
//
/////////////////////////////////////////////////////////////////////////

void TipsyParticleDistribute::findFileParticleCount()
{
  // Compute the total number of particles in the problem
  // Compute the maximum number of particles read on any one process to set buffer size

  // Open the tipsy standard file and abort if there is an error.
  this->tipsyInfile.open(this->inFiles[0].c_str(),"standard");
  if (!this->tipsyInfile.is_open()) 
    {
    cout << "Error opening file " << this->baseFile.c_str() << endl;
    cout.flush();
    return;
    }

  // Reading in the header
  tipsyInfile >> this->tipsyHeader;

  unsigned long pieceSize  = floor(this->tipsyHeader.h_nBodies*1.0/this->numProc);
  this->LocalStartIndex    = this->myProc*pieceSize;
  this->LocalEndIndex      = (myProc == this->numProc - 1) ? (this->tipsyHeader.h_nBodies-1) : (this->myProc+1)*pieceSize-1;
  this->LocalCount         = 1 + this->LocalEndIndex - this->LocalStartIndex;
  this->fileParticles.push_back(this->LocalCount);
  this->maxFiles = 1;

#ifdef USE_SERIAL_COSMO
  this->maxParticles = maxNumberOfParticles;
#else
  MPI_Allreduce((void*) &this->LocalCount,
                (void*) &this->maxParticles,
                1, MPI_LONG, MPI_MAX, Partition::getComm());
#endif

  cout << "Rank " << this->myProc << " opened " << this->inFiles[0].c_str() 
    << " total : " << this->tipsyHeader.h_nBodies << " particles" << endl;
  cout << "Rank " << this->myProc << " will read " << this->LocalStartIndex << " to " << this->LocalEndIndex << endl;
  cout << "Rank " << this->myProc << " maxParticles " << this->maxParticles << endl;
  cout.flush();
}

/////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------
tipsypos::section_type SeekToIndex(unsigned long index, TipsyHeader& tipsyHeader, ifTipsy& tipsyInfile)
{
	if(index < tipsyHeader.h_nSph)
		{
		// we are seeking a gas particle
		tipsyInfile.seekg(tipsypos(tipsypos::gas,index));	
		return tipsypos::gas;
		}
	else if(index < tipsyHeader.h_nDark+tipsyHeader.h_nSph)
		{
		// we are seeking a dark particle, which begin at index 0 after
		// gas particles
		index-=tipsyHeader.h_nSph;
		tipsyInfile.seekg(tipsypos(tipsypos::dark,index));
		return tipsypos::dark;
		}
	else if(index < tipsyHeader.h_nBodies)
		{
		// we are seeking a star particle, which begin at index zero after gas
		// and dark particles
		index-=(tipsyHeader.h_nSph+tipsyHeader.h_nDark);
		tipsyInfile.seekg(tipsypos(tipsypos::star,index));	
		return tipsypos::star;
		}
	else
		{
//		vtkErrorMacro("An index is greater than the number of particles in the 	file, unable to read");
		return tipsypos::invalid;
		}
}
/////////////////////////////////////////////////////////////////////////////
void TipsyParticleDistribute::readFromTipsyFile(
      ifTipsy* tipsyInfile,   // file to read from
      int firstParticle,      // First particle index to read in this chunk
      int numberOfParticles,  // Number of particles to read in this chunk
      POSVEL_T* fBlock,       // Buffer for read in data
      ID_T* iBlock,           // Buffer for read in data
      Message* message)       // Message buffer for distribution
{
  cout << "Rank " << this->myProc << " reading particles from " << firstParticle << "  to " <<
       firstParticle + numberOfParticles-1 << " " << endl;
    cout.flush();

  // Store number of particles used in first position
  message->putValue(&numberOfParticles);
  if (numberOfParticles == 0) return;

  // Store each particle location, velocity, mass and tag (as float) in buffer
  int changeCount = 0;
  for (int p = 0; p < numberOfParticles; p++) {

    iBlock[0] = p;

    TipsyGasParticle  g;
    TipsyDarkParticle d;
    TipsyStarParticle s;
	  tipsypos::section_type particleSection=SeekToIndex(p, this->tipsyHeader, *tipsyInfile);
    static bool error_reported = false;
	  switch(particleSection)
		  {
     case tipsypos::gas:
       *tipsyInfile >> g;
      fBlock[0] = -g.pos[0];
      fBlock[2] = -g.pos[1];
      fBlock[4] = -g.pos[2];
      fBlock[1] = g.vel[0];
      fBlock[3] = g.vel[1];
      fBlock[5] = g.vel[2];
      fBlock[6] = g.mass;
       break;
     case tipsypos::dark:
       *tipsyInfile >> d;
      fBlock[0] = -d.pos[0];
      fBlock[2] = -d.pos[1];
      fBlock[4] = -d.pos[2];
      fBlock[1] = d.vel[0];
      fBlock[3] = d.vel[1];
      fBlock[5] = d.vel[2];
      fBlock[6] = d.mass;
       break;
     case tipsypos::star:
       *tipsyInfile >> s;
      fBlock[0] = -s.pos[0];
      fBlock[2] = -s.pos[1];
      fBlock[4] = -s.pos[2];
      fBlock[1] = s.vel[0];
      fBlock[3] = s.vel[1];
      fBlock[5] = s.vel[2];
      fBlock[6] = s.mass;
       break;
     default:
       // only do this once
       if (!error_reported) {
        error_reported = true;
       }
       break;
      }


    // If the location is not within the bounding box wrap around
    for (int i = 0; i <= 4; i = i + 2) {
      if (this->RecentreBoundingBox) {
        fBlock[i] += this->boxSize/2.0;
      }
//      fBlock[i] = (rand()%20)/20.0;

      if (fBlock[i] >= this->boxSize) {
#ifndef DEBUG_MESSAGES
        cout << "Location at " << i << " changed from " << fBlock[i] << endl;
    cout.flush();
#endif
        fBlock[i] -= this->boxSize;
        changeCount++;
      }
    }

    // Store location and velocity and mass in message buffer
    // Reorder so that location vector is followed by velocity vector
    message->putValue(&fBlock[0]);
    message->putValue(&fBlock[2]);
    message->putValue(&fBlock[4]);
    message->putValue(&fBlock[1]);
    message->putValue(&fBlock[3]);
    message->putValue(&fBlock[5]);
    message->putValue(&fBlock[6]);

    // Store the integer tag
    message->putValue(&iBlock[0]);
  }

#ifndef DEBUG_MESSAGES
  cout << "Rank " << this->myProc << " wrapped around " << changeCount
       << " particles" << endl;
    cout.flush();
#endif
}

/////////////////////////////////////////////////////////////////////////
//
// Each processor reads 1 file or gets a pointer to data eventually
// As the particle is read it will be stored as an alive particle on this
// processor and will be checked about neighbor ranges to see if it must
// be exchanged
//
/////////////////////////////////////////////////////////////////////////
  
void TipsyParticleDistribute::readParticlesOneToOne(int reserveQ)
{
  // File name is the base file name with processor id appended
  // Because MPI Cartesian topology is used the arrangement of files in
  // physical space must follow the rule of last dimension varies fastest
  ostringstream fileName;
  fileName << this->baseFile << "." << this->myProc;
  this->inFiles.push_back(fileName.str());

  // Compute the total number of particles in the problem
  // Compute the maximum number of particles in any one file to set buffer size
  findFileParticleCount();
  
  // Reserve particle storage to minimize reallocation
  int reserveSize = (int) (this->maxParticles * DEAD_FACTOR);
  
  if(reserveQ) {
#ifndef DEBUG_MESSAGES
    cout << "readParticlesOneToOne reserving vectors" << endl;
    cout.flush();
#endif
    this->xx->reserve(reserveSize);
    this->yy->reserve(reserveSize);
    this->zz->reserve(reserveSize);
    this->vx->reserve(reserveSize);
    this->vy->reserve(reserveSize);
    this->vz->reserve(reserveSize);
    this->ms->reserve(reserveSize);
    this->tag->reserve(reserveSize);
  }

  // Running total and index into particle data on this processor
  this->particleCount = 0;
  
  // Read the input file storing particles immediately because all are alive
  if (this->inputType == RECORD) {
    readFromRecordFile();
  } else if (this->inputType == BLOCK) {
    readFromBlockFile();
  } else if (this->inputType == CUSTOM) {
    readFromTipsyFile();
  }
}


void TipsyParticleDistribute::readFromTipsyFile()
{
}
