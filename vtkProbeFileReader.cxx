
#include "vtkProbeFileReader.h"

#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkObjectFactory.h"
#include "vtkFloatArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkFieldData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkMath.h"
#include "vtkDiskSource.h"
#include "vtkPlaneSource.h"
#include "vtkRegularGridSource.h"
#include "vtkStructuredGrid.h"
#include "vtkLineSource.h"
#include "vtkPointSource.h"
#include "vtkDelimitedTextReader.h"
#include "vtkTable.h"
#include "vtkVariantArray.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkPlane.h"
#include "vtkBox.h"
#include "vtkImplicitBoolean.h"

#include <vtksys/SystemTools.hxx>

#define GRID_MODE

#ifdef GRID_MODE
 #define PlaneSourceType vtkRegularGridSource
#else
 #define PlaneSourceType vtkPlaneSource
#endif

#define EPSILON 0.000001


//----------------------------------------------------------------------------
vtkSmartPointer<vtkTransform> RotationMatrix(double from[3], double to[3], double axis[3])
{
  double rotVector[3], theta;
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();

  // Compute rotation vector using a transformation matrix.
  // Note that if normals are parallel then the rotation is either
  // 0 or 180 degrees.
  double dp = vtkMath::Dot(from, to);
  if ( dp >= 1.0 )
    {
    return transform; //zero rotation
    }
  else if ( dp <= -1.0 )
    {
    theta = 180.0;
    rotVector[0] = axis[0];
    rotVector[1] = axis[1];
    rotVector[2] = axis[2];
    }
  else
    {
    vtkMath::Cross(from, to, rotVector);
    theta = vtkMath::DegreesFromRadians(acos((double)dp));
    }

  // create rotation matrix
  transform->PostMultiply();
  transform->RotateWXYZ(theta, rotVector);

  return transform;
}
//----------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkProbeFileReader, "$Revision: 557 $");
vtkStandardNewMacro(vtkProbeFileReader);
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkSmartPointer<vtkDataSet> vtkProbeFileReader::ESPHI_Copy(vtkDataSet *d) {
  vtkSmartPointer<vtkDataSet> result;
  result.TakeReference(d->NewInstance());
  result->ShallowCopy(d);
  result->CopyInformation(d);
  result->SetSource(NULL);
  return result;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkDataSet> vtkProbeFileReader::ESPHI_CopyOutput(vtkAlgorithm *a, int port) {
  a->Update();
  vtkDataSet *dobj = vtkDataSet::SafeDownCast(a->GetOutputDataObject(port));
  return ESPHI_Copy(dobj);
}
//----------------------------------------------------------------------------
void TagDataSet(vtkSmartPointer<vtkDataSet> output, vtkSmartPointer<vtkVariantArray> row, char *name)
{
  // tag the dataset with the data so we can use an implicit function 
  // if needed when probing mesh based data
  //
  vtkSmartPointer<vtkVariantArray> rowCopy = vtkSmartPointer<vtkVariantArray>::New();
  rowCopy->DeepCopy(row);
  rowCopy->SetName(name);
  output->GetFieldData()->AddArray(rowCopy);
}
//---------------------------------------------------------------------------
vtkSmartPointer<vtkTransform> TransformFromVectorToVector(double from[3], double to[3])
{
  double xaxis[3] = {1,0,0};
  vtkSmartPointer<vtkTransform> trans = RotationMatrix(from, to, xaxis);
  return trans;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkProbeFileReader::ESPHI_TransformCopy(vtkAlgorithm *a, int port, vtkTransform *trans) {
  vtkSmartPointer<vtkTransformPolyDataFilter> tranfilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  tranfilter->SetInputConnection(a->GetOutputPort(port));
  tranfilter->SetTransform(trans);
  vtkSmartPointer<vtkPolyData> result = vtkPolyData::SafeDownCast(ESPHI_CopyOutput(tranfilter, 0));
  return result;
}
//---------------------------------------------------------------------------
//----------------------------------------------------------------------------
#define BOXWIDGET_XMY(x,y,z) \
  z[0] = x[0] - y[0]; \
  z[1] = x[1] - y[1]; \
  z[2] = x[2] - y[2]; 
#define BOXWIDGET_AXPY(a,x,y,z) \
  z[0] = a*x[0] + y[0]; \
  z[1] = a*x[1] + y[1]; \
  z[2] = a*x[2] + y[2]; 
#define BOXWIDGET_XMAY(a,x,y,z) \
  z[0] = x[0] - a*y[0]; \
  z[1] = x[1] - a*y[1]; \
  z[2] = x[2] - a*y[2]; 
#define BOXWIDGET_PAX(a,x,z) \
  z[0] += a*x[0]; \
  z[1] += a*x[1]; \
  z[2] += a*x[2]; 
#define BOXWIDGET_PAXPBY(a,x,b,y,z) \
  z[0] += a*x[0] + b*y[0]; \
  z[1] += a*x[1] + b*y[1]; \
  z[2] += a*x[2] + b*y[2]; 
#define BOXWIDGET_AVERAGE(a,b,c) \
  c[0] = (a[0] + b[0])/2.0; \
  c[1] = (a[1] + b[1])/2.0; \
  c[2] = (a[2] + b[2])/2.0;
//---------------------------------------------------------------------------
//----------------------------------------------------------------------------
void CreateAxes(double centre[3], double v0[3], double v1[3], double norm[3])
{
  double vec1[3], vec2[3];
  vtkMath::Perpendiculars(norm, vec1, vec2, 0.0);
  //
  double sx,sy;
  sx = vtkMath::Norm(v0);
  sy = vtkMath::Norm(v1);
  BOXWIDGET_AXPY(sx, vec1, centre, v0);
  BOXWIDGET_AXPY(sy, vec2, centre, v1);
  // The centre is passed in, not the corner, 
  // so offset all 3 by half the vector lengths
  BOXWIDGET_XMAY(sx/2.0,centre,vec1,centre);
  BOXWIDGET_XMAY(sy/2.0,centre,vec2,centre);
  BOXWIDGET_XMAY(sx/2.0,v0,vec1,v0);
  BOXWIDGET_XMAY(sy/2.0,v0,vec2,v0);
  BOXWIDGET_XMAY(sx/2.0,v1,vec1,v1);
  BOXWIDGET_XMAY(sy/2.0,v1,vec2,v1);
}
//---------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkProbeFileReader::vtkProbeFileReader()
{
  this->ProbeFileName = NULL;
  this->Resolution1 = 20;
  this->Resolution2 = 20;
  this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
vtkProbeFileReader::~vtkProbeFileReader()
{
  delete [] this->ProbeFileName;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkImplicitFunction> vtkProbeFileReader::CreateImplicitPlane(
  vtkSmartPointer<vtkVariantArray> row)
{
  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
  for (int i=0; i<9; i++) if (row->GetValue(i).ToString().size()==0) {
    return plane;
  }
  vtkVariant label = row->GetValue(0);
  double O[3] = {
    row->GetValue(1).ToDouble(),
    row->GetValue(2).ToDouble(),
    row->GetValue(3).ToDouble() };
  double N[3] = { 
    row->GetValue(4).ToDouble(), 
    row->GetValue(5).ToDouble(), 
    row->GetValue(6).ToDouble() };
  double sx = row->GetValue(7).ToDouble();
  double sy = row->GetValue(8).ToDouble();
  //
  plane->SetNormal(N);
  plane->SetOrigin(O);
  //
  return plane;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkImplicitFunction> vtkProbeFileReader::CreateImplicitBox(
  vtkSmartPointer<vtkVariantArray> row, double bounds[6], double scalefactor)
{
  vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
  for (int i=0; i<9; i++) if (row->GetValue(i).ToString().size()==0) {
    return box;
  }
  double centre[3] = {
    row->GetValue(1).ToDouble(),
    row->GetValue(2).ToDouble(),
    row->GetValue(3).ToDouble() };
  double N[3] = { 
    row->GetValue(4).ToDouble(), 
    row->GetValue(5).ToDouble(), 
    row->GetValue(6).ToDouble() };
  double sx = row->GetValue(7).ToDouble();
  double sy = row->GetValue(8).ToDouble();
  //
  // The vtkBox implicit function is axis aligned
  // so we setup the box using unit axes, then apply a transform
  // to bring it along our user defined axes
  //
  double xaxis[3] = {1,0,0}, yaxis[3] = {0,1,0}, zaxis[3] = {0,0,1};
  double thickness = 0.5*sqrt(sx*sx + sy*sy)/scalefactor; // just to make sure it's non zero
  double mincorner[3] = {centre[0], centre[1], centre[2]-thickness};
  double maxcorner[3] = {centre[0], centre[1], centre[2]+thickness};
  BOXWIDGET_PAXPBY(-sx/2.0,xaxis,-sy/2.0,yaxis, mincorner);
  BOXWIDGET_PAXPBY( sx/2.0,xaxis, sy/2.0,yaxis, maxcorner);
  //
  box->SetXMin(mincorner);
  box->SetXMax(maxcorner);
  //
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->Translate( centre[0],  centre[1],  centre[2]);
  vtkSmartPointer<vtkLinearTransform> rotation = RotationMatrix(zaxis, N, xaxis);
  transform->Concatenate(rotation);
  transform->Translate(-centre[0], -centre[1], -centre[2]);
  vtkSmartPointer<vtkLinearTransform> inverse = transform->GetLinearInverse();
  box->SetTransform(inverse);
  //
  // Compute the bounding box in transformed coordinates so that the filters
  // using this box can do a fast BBox test for quick rejection.
  // Unit length sided box
  double p0[3] = {-0.5, -0.5, -0.5};
  double p1[3] = { 0.5,  0.5,  0.5};
  vtkSmartPointer<vtkTransform> bboxTransform = vtkSmartPointer<vtkTransform>::New();
  // we can make our box much tighter than the clipping box as we are doing true 
  // geometry tests - the clipper needs more tolerance to fully include cells
  bboxTransform->Scale( sx,  sy,  thickness/5.0); 
  bboxTransform->Concatenate(rotation);
  bboxTransform->Translate(centre[0], centre[1], centre[2]);
  bboxTransform->TransformPoint(p0,p0);
  bboxTransform->TransformPoint(p1,p1);
  bounds[0] = p0[0]<p1[0] ? p0[0] : p1[0];
  bounds[1] = p0[0]<p1[0] ? p1[0] : p0[0];
  bounds[2] = p0[1]<p1[1] ? p0[1] : p1[1];
  bounds[3] = p0[1]<p1[1] ? p1[1] : p0[1];
  bounds[4] = p0[2]<p1[2] ? p0[2] : p1[2];
  bounds[5] = p0[2]<p1[2] ? p1[2] : p0[2];
  return box;
}
//----------------------------------------------------------------------------
/*
vtkSmartPointer<vtkImplicitFunction> vtkProbeFileReader::CreateImplicitBox(
  vtkSmartPointer<vtkVariantArray> row)
{
  vtkSmartPointer<vtkImplicitFunction> result;
  for (int i=0; i<9; i++) if (row->GetValue(i).ToString().size()==0) {
    return NULL;
  }
  //
  vtkVariant label = row->GetValue(0);
  double O[3] = {
    row->GetValue(1).ToDouble(),
    row->GetValue(2).ToDouble(),
    row->GetValue(3).ToDouble() };
  double centre[3] = {O[0], O[1], O[2]};
  double N[3] = { 
    row->GetValue(4).ToDouble(), 
    row->GetValue(5).ToDouble(), 
    row->GetValue(6).ToDouble() };
  double sx = row->GetValue(7).ToDouble();
  double sy = row->GetValue(8).ToDouble();
  //
  double xvec[3], yvec[3];
  vtkMath::Perpendiculars(N, xvec, yvec, 0.0);

  // The vtkBox implicit function is axis aligned
  // so we setup the box using unit axes, then apply a transform
  // to bring it along our user defined axes
  vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
  double xaxis[3] = {1,0,0}, yaxis[3] = {0,1,0}, zaxis[3] = {0,0,1};
  double thickness = sqrt(sx*sx + sy*sy)/10.0; // just to make sure it's non zero
  double mincorner[3] = {centre[0], centre[1], centre[2]-thickness};
  double maxcorner[3] = {centre[0], centre[1], centre[2]+thickness};
  BOXWIDGET_PAXPBY(-sx/2.0,xaxis,-sy/2.0,yaxis, mincorner);
  BOXWIDGET_PAXPBY( sx/2.0,xaxis, sy/2.0,yaxis, maxcorner);
  //
  box->SetXMin(mincorner);
  box->SetXMax(maxcorner);
  //
  double norm[3];
  vtkMath::Cross(xvec, yvec, norm);
  vtkMath::Normalize(norm);
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->Translate(centre[0], centre[1], centre[2]);
  vtkSmartPointer<vtkLinearTransform> rotation = RotationMatrix(zaxis, norm, xaxis);
  transform->Concatenate(rotation);
  transform->Translate(-centre[0], -centre[1], -centre[2]);
  vtkSmartPointer<vtkLinearTransform> inverse = transform->GetLinearInverse();
  box->SetTransform(transform);
  return box;
}
*/
//----------------------------------------------------------------------------
vtkSmartPointer<vtkTable> vtkProbeFileReader::ReadProbeTableData()
{
  vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
  reader->SetFileName(this->ProbeFileName);
  reader->SetFieldDelimiterCharacters(" ");
  reader->SetMergeConsecutiveDelimiters(true);
  reader->Update();
  return reader->GetOutput();
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkDataSet> vtkProbeFileReader::MakeOneRectProbe(vtkSmartPointer<vtkVariantArray> row)
{
  for (int i=0; i<9; i++) if (row->GetValue(i).ToString().size()==0) {
    return NULL;
  }
  //
  vtkVariant label = row->GetValue(0);
  double O[3] = {
    row->GetValue(1).ToDouble(),
    row->GetValue(2).ToDouble(),
    row->GetValue(3).ToDouble() };
  double N[3] = { 
    row->GetValue(4).ToDouble(), 
    row->GetValue(5).ToDouble(), 
    row->GetValue(6).ToDouble() };
  double sx = row->GetValue(7).ToDouble();
  double sy = row->GetValue(8).ToDouble();
  //
  double p0[3] = {sx, 0.0, 0.0 };
  double p1[3] = {0.0, sy, 0.0};
  CreateAxes(O, p0, p1, N);
  //
  vtkSmartPointer<PlaneSourceType> plane = vtkSmartPointer<PlaneSourceType>::New();
  plane->SetOrigin(O);
  plane->SetPoint1(p0);
  plane->SetPoint2(p1);
#ifdef GRID_MODE
  int resolution[3] = { this->Resolution1, this->Resolution2, 1};
  plane->SetResolution(resolution);
  plane->SetGenerateConnectedCells(1);
#else
  plane->SetXResolution(this->Resolution1);
  plane->SetYResolution(this->Resolution2);
  TagDataSet(result, row, "ImplicitPlaneData");
#endif
  vtkSmartPointer<vtkDataSet> result = vtkDataSet::SafeDownCast(ESPHI_CopyOutput(plane, 0));
  return result;
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkProbeFileReader::MakeOneDiskProbe(vtkSmartPointer<vtkVariantArray> row)
{
  for (int i=0; i<8; i++) if (row->GetValue(i).ToString().size()==0) {
    return NULL;
  }
  //
  vtkVariant label = row->GetValue(0);
  double O[3] = {
    row->GetValue(1).ToDouble(),
    row->GetValue(2).ToDouble(),
    row->GetValue(3).ToDouble() };
  double N[3] = { 
    row->GetValue(4).ToDouble(), 
    row->GetValue(5).ToDouble(), 
    row->GetValue(6).ToDouble() };
  double r = row->GetValue(7).ToDouble();
  double zvec[3] = {0,0,1};
  vtkMath::Normalize(N);
  vtkSmartPointer<vtkTransform> trans = TransformFromVectorToVector(zvec, N);
  trans->Translate(O[0], O[1], O[2]);
  //
  vtkSmartPointer<vtkDiskSource> disk = vtkSmartPointer<vtkDiskSource>::New();
  disk->SetRadialResolution(this->Resolution1);
  disk->SetCircumferentialResolution(this->Resolution2);
  disk->SetInnerRadius(0.0);
  disk->SetOuterRadius(r);
  //
  vtkSmartPointer<vtkPolyData> result = ESPHI_TransformCopy(disk, 0, trans);
  TagDataSet(result, row, "ImplicitDiskData");
  return result;
}
//---------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkProbeFileReader::MakeOneLineProbe(vtkSmartPointer<vtkVariantArray> row)
{
  for (int i=0; i<7; i++) if (row->GetValue(i).ToString().size()==0) {
    return NULL;
  }
  //
  vtkVariant label = row->GetValue(0);
  double P1[3] = {
    row->GetValue(1).ToDouble(),
    row->GetValue(2).ToDouble(),
    row->GetValue(3).ToDouble() };
  double P2[3] = { 
    row->GetValue(4).ToDouble(), 
    row->GetValue(5).ToDouble(), 
    row->GetValue(6).ToDouble() };
  vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();
  line->SetPoint1(P1);
  line->SetPoint2(P2);
  line->SetResolution(this->Resolution1);
  //
  vtkSmartPointer<vtkPolyData> result = vtkPolyData::SafeDownCast(ESPHI_CopyOutput(line, 0));
  TagDataSet(result, row, "ImplicitLineData");
  return result;
}
//---------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkProbeFileReader::MakeOnePointProbe(vtkSmartPointer<vtkVariantArray> row)
{
  for (int i=0; i<4; i++) if (row->GetValue(i).ToString().size()==0) {
    return NULL;
  }
  //
  vtkVariant label = row->GetValue(0);
  double P1[3] = {
    row->GetValue(1).ToDouble(),
    row->GetValue(2).ToDouble(),
    row->GetValue(3).ToDouble() };
  vtkSmartPointer<vtkPointSource> point = vtkSmartPointer<vtkPointSource>::New();
  point->SetCenter(P1);
  point->SetRadius(0);
  point->SetNumberOfPoints(1);
  //
  vtkSmartPointer<vtkPolyData> result = vtkPolyData::SafeDownCast(ESPHI_CopyOutput(point, 0));
  TagDataSet(result, row, "ImplicitPointData");
  return result;
}
//---------------------------------------------------------------------------
int vtkProbeFileReader::RequestData(vtkInformation* request,
                          vtkInformationVector** vtkNotUsed(inputVector),
                          vtkInformationVector* outputVector)
{
  vtkDebugMacro( << "" << "Entering RequestData" )
  //
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataObject* doOutput = outInfo->Get(vtkDataObject::DATA_OBJECT());
  vtkMultiBlockDataSet *group = vtkMultiBlockDataSet::SafeDownCast(doOutput);
  if (!group) {
    vtkErrorMacro("No output data in RequestData");
    return 0;
  }

  if (!outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) ||
      !outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()))
  {
    vtkErrorMacro("Expected information not found in RequestData");
    return 0;
  }

  int updatePiece     = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int updateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  vtkDebugMacro( << "Num pieces " << updateNumPieces << " This piece " << updatePiece )

  //
  // Read probe data
  //
  if (!this->ProbeFileName || !vtksys::SystemTools::FileExists(this->ProbeFileName))
  {
    return 0;
  }

  //----------------------------------------------------------------------------------//
  // Read Probe data
  //----------------------------------------------------------------------------------//
  vtkSmartPointer<vtkTable> table = this->ReadProbeTableData();

  int Id=0, pType = -1;
  for (int i=0; i<table->GetNumberOfRows(); i++) {
    vtkVariantArray *row = table->GetRow(i);
    vtkVariant entry = row->GetValue(0);
    if (entry.ToString()==vtkstd::string("RectProbe")) {
      pType = ESPHI_RECT;
    }
    else if (entry.ToString()==vtkstd::string("DiskProbe")) {
      pType = ESPHI_DISC;
    }
    else if (entry.ToString()==vtkstd::string("LineProbe")) {
      pType = ESPHI_LINE;
    }
    else if (entry.ToString()==vtkstd::string("PointProbe")) {
      pType = ESPHI_POINT;
    }
    // skip comments, lines with # at the start
    else if (entry.ToString().c_str()[0]=='#') continue;
    //
    else if (pType!=-1) {      
      vtkSmartPointer<vtkDataSet> probe;
      switch (pType) {
        case ESPHI_RECT:
          probe = this->MakeOneRectProbe(row);
          // Label X Y Z(Center) Nx Ny Nz(normal) X_size Y_size
          break;
        case ESPHI_DISC:
          probe = this->MakeOneDiskProbe(row);
          break;
        case ESPHI_LINE:
          probe = this->MakeOneLineProbe(row);
          break;
        case ESPHI_POINT:
          probe = this->MakeOnePointProbe(row);
          break;
      }
      if (probe) {
        group->SetBlock(Id, probe);
        vtkStdString name = row->GetValue(0).ToString();
        group->GetMetaData(Id++)->Set(vtkCompositeDataSet::NAME(), name.c_str());
      }
    }
  }

  return 1;
}

//----------------------------------------------------------------------------
void vtkProbeFileReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


