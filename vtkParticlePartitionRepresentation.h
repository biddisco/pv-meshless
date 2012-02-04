/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkParticlePartitionRepresentation.h
  Revision of last commit : $Rev: 155 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2006-07-13 10:23:31 +0200 #$

  Copyright (c) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing this
  copyright notice appears on all copies of source code and an
  acknowledgment appears with any substantial usage of the code.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
// .NAME vtkParticlePartitionRepresentation - Generate boxes from partitioned regions of data
//
// .SECTION Description
//
// .SECTION See Also
// 

#ifndef _vtkParticlePartitionRepresentation_h
#define _vtkParticlePartitionRepresentation_h

#include "vtkPolyDataAlgorithm.h"

class VTK_EXPORT vtkParticlePartitionRepresentation : public vtkPolyDataAlgorithm {
  public:
    // Description:
    // Standard Type-Macro
    vtkTypeRevisionMacro(vtkParticlePartitionRepresentation,vtkPolyDataAlgorithm);

    // Description:
    // Create an instance of vtkParticlePartitionRepresentation
    static vtkParticlePartitionRepresentation *New();

    // Description:
    // By default, each process generates only the box for its own rank
    // when this option is enabled, each node produces all N boxes which may be
    // useful for some debugging purposes
    vtkSetMacro(AllBoxesOnAllProcesses,int);
    vtkGetMacro(AllBoxesOnAllProcesses,int);
    vtkBooleanMacro(AllBoxesOnAllProcesses,int);

    // Description:
    // Specify a factor (default=1) to increase (>1) or decrease (<1)
    // the generated boxes by.
    vtkSetMacro(InflateFactor,double);
    vtkGetMacro(InflateFactor,double);
    
  protected:
     vtkParticlePartitionRepresentation(void);
    ~vtkParticlePartitionRepresentation();
    //
    virtual int RequestInformation (vtkInformation *,
                                    vtkInformationVector **,
                                    vtkInformationVector *);
    //
    virtual int RequestData(vtkInformation *request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);
    //
    virtual int FillInputPortInformation(int port, vtkInformation* info);

    int    AllBoxesOnAllProcesses;
    double InflateFactor;
private:
  // Not implemented.
  vtkParticlePartitionRepresentation(const vtkParticlePartitionRepresentation&);
  void operator=(const vtkParticlePartitionRepresentation&);
};

#endif

