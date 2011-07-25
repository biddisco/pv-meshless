/*=========================================================================

  Project:   RPD
  Module:    $RCSfile: vtkRegularGridSource.h,v $
  Date:      $Date: 2004/02/06 07:31:49 $
  Version:   $Revision: 1.2 $

  Copyright (C) 2000-2002 Skipping Mouse Software Ltd.
  All Rights Reserved.

  Source code from Skipping Mouse Software is supplied under the terms of a
  license agreement and may not be copied or disclosed except in accordance
  with the terms of that agreement. This file is subject to the license
  found in the file Copyright.txt supplied with the software.

=========================================================================*/
// .NAME vtkRegularGridSource - Generate points in a box-like volume
//
// .SECTION Description
// vtkRegularGridSource creates points along a regular grid
// in3 dimensions. The grid is defined by the X/Y/Z Spacing fields
// and a supplied bounding volume. In general the bounding volume
// is taken from a supplied input dataset.
//
// .SECTION See Also
// vtkGeneratePointsFilter

#ifndef _vtkRegularGridSource_h
#define _vtkRegularGridSource_h

#include "vtkDataSetAlgorithm.h"

class VTK_EXPORT vtkRegularGridSource : public vtkDataSetAlgorithm {
  public:
    // Description:
    // Standard Type-Macro
    static vtkRegularGridSource *New();
    vtkTypeRevisionMacro(vtkRegularGridSource, vtkDataSetAlgorithm);

    // Description:
    // The source can be controlled by manually, by setting the axes
    // points. If ou wish the source to ignore the points and simply
    // use the bounds of an input to generate a resampling grid, then
    // set AutoPlacement to true
    vtkSetMacro(UseAutoPlacement,int);
    vtkGetMacro(UseAutoPlacement,int);
    vtkBooleanMacro(UseAutoPlacement,int);

    // Description:
    // By default, the filter will create an array of points which fill
    // the grid - the output will be a PolyData object with vertices at 
    // each point. When GenerateConnectedCells is set, the output will
    // be a PolyData Plane in plane mode, or a StructuredGrid of 3D cells
    // in Box mode.
    vtkSetMacro(GenerateConnectedCells,int);
    vtkGetMacro(GenerateConnectedCells,int);
    vtkBooleanMacro(GenerateConnectedCells,int);

    // Description
    // This point defines the origin of the box of points
    // if no input is provided - if an input is supplied, 
    // the bounds of the input are taken
    vtkSetVector3Macro(Origin, double);
    vtkGetVector3Macro(Origin, double);

    // Description
    // This point defines x axis of the box of points
    // if no input is provided - if an input is supplied, 
    // the bounds of the input are taken
    vtkSetVector3Macro(Point1, double);
    vtkGetVector3Macro(Point1, double);

    // Description
    // This point defines the y axis of the box of points 
    // if no input is provided - if an input is supplied, 
    // the bounds of the input are taken
    vtkSetVector3Macro(Point2, double);
    vtkGetVector3Macro(Point2, double);

    // Description
    // This point defines the z axis of the box of points 
    // if no input is provided - if an input is supplied, 
    // the bounds of the input are taken
    vtkSetVector3Macro(Point3, double);
    vtkGetVector3Macro(Point3, double);

    // Description
    // Get the Normal
    virtual double *GetNormal();

    // Description:
    // Specify the point spacing on the X/Y/Z axis
    vtkSetVector3Macro(Spacing, double);
    vtkGetVector3Macro(Spacing, double);

    // Description
    // If Resolution[X/Y/Z]are all Non-zero, then
    // the spacing is ignored and the box defined by the points is 
    // sampled using the specified resolutions
    vtkSetVector3Macro(Resolution, int);
    vtkGetVector3Macro(Resolution, int);

    // Description
    // Add a delta to each axis to expand the box
    vtkSetMacro(Delta, double);
    vtkGetMacro(Delta, double);

    void ComputeOBB();

  protected:
    vtkRegularGridSource();

    // Description:
    // This is a special function which attaches a variant array to the
    // field data of the output. The array contains the origin/normal 
    // and x/y size of the data when a plane is generated. 
    // This can be used by other custom filters to generate an implicit
    // function for slicing instead of using the plane geometry.
    void TagDataSet(vtkDataSet *output, char *name);

    // Pipeline mechanism
    virtual int FillInputPortInformation(int port, vtkInformation *info);
    virtual int FillOutputPortInformation(int port, vtkInformation* info);

    // Description:
    // A utility function to help us decide what type of data to export
    virtual int RequiredDataType();

    // Description:
    // This is called within ProcessRequest when a request asks for 
    // Information. Typically an algorithm provides whatever lightweight 
    // information about its output that it can here without doing any 
    // lengthy computations. This happens in the first pass of the pipeline
    // execution.
    virtual int RequestInformation(vtkInformation*, 
                                   vtkInformationVector**, 
                                   vtkInformationVector*);
    
    //// Description:
    virtual int RequestUpdateExtent(vtkInformation*,
                                    vtkInformationVector**,
                                    vtkInformationVector*);

    // Description:
    // This is called within ProcessRequest to when a request asks the
    // algorithm to create empty output data objects. This typically happens
    // early on in the execution of the pipeline. The default behavior is to 
    // create an output DataSet of the same type as the input for each 
    // output port. This method can be overridden to change the output 
    // data type of an algorithm. This happens in the third pass of the 
    // pipeline execution.
    virtual int RequestDataObject(vtkInformation* request, 
                                  vtkInformationVector** inputVector, 
                                  vtkInformationVector* outputVector);
    
    // Description:
    // This is called within ProcessRequest when a request asks the algorithm
    // to do its work. This is the method you should override to do whatever the
    // algorithm is designed to do. This happens during the fourth pass in the
    // pipeline execution process.
    virtual int RequestData(vtkInformation*, 
                            vtkInformationVector**, 
                            vtkInformationVector*);

    virtual int ComputeInformation(vtkInformation *,
                                   vtkInformationVector **,
                                   vtkInformationVector *);

    virtual void ComputeAxesFromBounds(vtkDataSet *inputData, double lengths[3], bool inflate);
    virtual void ComputeAxesFromPoints(double lengths[3], bool inflate);
    virtual void BoundsToExtent(double *bounds, int *extent, int updatePiece);

    // properties which may be set by user
    double   Spacing[3];
    int      Resolution[3];
    double   Origin[3];
    double   Point1[3];
    double   Point2[3];
    double   Point3[3];
    double   Delta;
    int      GenerateConnectedCells;
    int      UseAutoPlacement;

    // internal values which may differ from above ones
    int      WholeDimension[3];
    double   origin[3];
    double   axesvectors[3][3];
    double   scaling[3];
    double   spacing[3];

    // only for export, do not use
    double   normal[3];
private:
  vtkRegularGridSource(const vtkRegularGridSource&);  // Not implemented.
  void operator=(const vtkRegularGridSource&);  // Not implemented.
};

#endif
