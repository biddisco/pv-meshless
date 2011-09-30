/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkTipsyReader.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTipsyReader - Read points from a Tipsy standard binary file
// .SECTION Description
// Read points from a Tipsy standard binary file. Fully parallel. Has ability
// to read in additional attributes from an ascii file, and to only load in
// marked particles but both these functions are serial only.
#ifndef __vtkTipsyReader_h
#define __vtkTipsyReader_h

#include "vtkPolyDataAlgorithm.h" // superclass

#include "vtkSmartPointer.h"
#include "tipsylib/ftipsy.hpp" // functions take Tipsy particle objects
#include <vtkstd/vector>

class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkFloatArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;
class FileSeriesFinder;

class VTK_EXPORT vtkTipsyReader : public vtkPolyDataAlgorithm
{
public:
  static vtkTipsyReader* New();
  vtkTypeRevisionMacro(vtkTipsyReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to read the marked points.
  vtkSetStringMacro(MarkFileName);
  vtkGetStringMacro(MarkFileName);
  // Description:
  // Set/Get the name of the file from which to read points.
  vtkSetStringMacro(FileName);
   vtkGetStringMacro(FileName);

  // Description:
  // Get/Set whether to distribute data
  vtkSetMacro(DistributeDataOn,int);
  vtkGetMacro(DistributeDataOn,int);

  // Description:
  // Get the number of timesteps in the file
  vtkGetMacro(NumberOfTimeSteps,int);

  //BTX
  void GetTimeStepValues(std::vector<double> &steps);
  //ETX


  // Description:
  // An H5Part file may contain multiple arrays
  // a GUI (eg Paraview) can provide a mechanism for selecting which data arrays
  // are to be read from the file. The PointArray variables and members can
  // be used to query the names and number of arrays available
  // and set the status (on/off) for each array, thereby controlling which
  // should be read from the file. Paraview queries these point arrays after
  // the (update) information part of the pipeline has been updated, and before the
  // (update) data part is updated.
  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAll();
  void        EnableAll();
  void        Disable(const char* name);
  void        Enable(const char* name);
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }

// The BTX, ETX comments bracket the portion of the code which should not be
// attempted to wrap for use by python, specifically the code which uses
// C++ templates as this code is unable to be wrapped. DO NOT REMOVE. 
//BTX
protected:
  vtkTipsyReader();
  ~vtkTipsyReader();
  //
  char   *MarkFileName;
  char   *FileName;
  int     DistributeDataOn;
  int     NumberOfTimeSteps;
  
  //
  int RequestInformation(vtkInformation*,  vtkInformationVector**,
    vtkInformationVector*);

  int RequestData(vtkInformation*,vtkInformationVector**,
    vtkInformationVector*);

  vtkIdType                       ParticleIndex;
  vtkSmartPointer<vtkIdTypeArray> GlobalIds;
  vtkSmartPointer<vtkPoints>      Positions;
  vtkSmartPointer<vtkCellArray>   Vertices;

  vtkSmartPointer<vtkFloatArray>   Potential;
  vtkSmartPointer<vtkFloatArray>   Mass;
  vtkSmartPointer<vtkFloatArray>   EPS;
  vtkSmartPointer<vtkFloatArray>   RHO;
  vtkSmartPointer<vtkFloatArray>   Hsmooth;
  vtkSmartPointer<vtkFloatArray>   Temperature;
  vtkSmartPointer<vtkFloatArray>   Metals;
  vtkSmartPointer<vtkFloatArray>   Tform;
  vtkSmartPointer<vtkFloatArray>   Type;
  vtkSmartPointer<vtkFloatArray>   Velocity;

  //
  int           UpdatePiece;
  int           UpdateNumPieces;

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

private:
  vtkTipsyReader(const vtkTipsyReader&);  // Not implemented.
  void operator=(const vtkTipsyReader&);  // Not implemented.
  /* Help functions for reading */
  // Description:
  // Reads the Tipsy header. 
  TipsyHeader ReadTipsyHeader(ifTipsy& tipsyInfile);
  // Description:
  // Reads all particles of this piece from the Tipsy file
  void ReadAllParticles(TipsyHeader& tipsyHeader,
    ifTipsy& tipsyInfile,int piece,int numPieces,vtkPolyData* output);
  // Description:
  // reads in a particle (either gas, dark or star as appropriate) 
  //from the tipsy in file of this class
  vtkIdType ReadParticle(unsigned long index, TipsyHeader& tipsyHeader,
    ifTipsy& tipsyInfile, vtkPolyData* output);
  // Description:
  // reads variables common to all particles
  vtkIdType ReadBaseParticle(vtkPolyData* output, TipsyBaseParticle& b);
  // Description:
  // reads variables common to all gas particles
  vtkIdType ReadGasParticle(vtkPolyData* output, TipsyGasParticle& g);
  // Description:
  // reads variables common to all star particles
  vtkIdType ReadStarParticle(vtkPolyData* output, TipsyStarParticle& s);
  // Description:
  // reads variables common to all dark particles
  vtkIdType ReadDarkParticle(vtkPolyData* output, TipsyDarkParticle& d);
  // Description:
  // Reads only Marked particles from the tipsy file.
  // Must be called after function ReadMarkedParticleIndices.
  void ReadMarkedParticles(
    vtkstd::vector<unsigned long>& markedParticleIndices,
    TipsyHeader& tipsyHeader,ifTipsy& tipsyInfile,vtkPolyData* output);
  // Description:
  // Helper function to read seek to a given index before reading
  tipsypos::section_type SeekToIndex(unsigned long index,
    TipsyHeader& tipsyHeader, ifTipsy& tipsyInfile);
  // Description:
  // reads in an array of the indices of marked particles from a file, 
  // returns a queue of marked particles
  // which is empty if reading was unsucessful.
  vtkstd::vector<unsigned long> ReadMarkedParticleIndices(
    TipsyHeader& tipsyHeader,  ifTipsy& tipsyInfile);
  /* Helper functions for storing data in output vector*/
  // Description:
  // allocates all vtk arrays for Tipsy variables and places them 
  // in the output vector
  void AllocateAllTipsyVariableArrays(vtkIdType numBodies,
    vtkPolyData* output);

  FileSeriesFinder         *Finder;
  vtkstd::vector<double>    TimeStepValues;
//ETX

};
#endif
