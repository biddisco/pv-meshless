// .NAME vtkProbeFileReader -- 
// .SECTION Description
// vtkProbeFileReader is a polydata filter object 
 
#ifndef __vtkProbeFileReader_h
#define __vtkProbeFileReader_h

#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkSmartPointer.h"
#include <vtkstd/vector>
#include <vtkstd/string>

class vtkTable;
class vtkDataSet;
class vtkPolyData;
class vtkVariantArray;
class vtkImplicitFunction;
class vtkTransform;

class VTK_EXPORT vtkProbeFileReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkProbeFileReader *New();
  vtkTypeRevisionMacro(vtkProbeFileReader, vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
    
  vtkSetStringMacro(ProbeFileName);
  vtkGetStringMacro(ProbeFileName);
      
  vtkSetMacro(Resolution1, int);
  vtkGetMacro(Resolution1, int);
  
  vtkSetMacro(Resolution2, int);
  vtkGetMacro(Resolution2, int);

//BTX
  enum probetypes {
    ESPHI_RECT  = 0,
    ESPHI_DISC  = 1,
    ESPHI_POINT = 2,
    ESPHI_LINE  = 3
  };

  static vtkSmartPointer<vtkDataSet> ESPHI_Copy(vtkDataSet *d);
  static vtkSmartPointer<vtkDataSet> ESPHI_CopyOutput(vtkAlgorithm *a, int port);
  static vtkSmartPointer<vtkPolyData> ESPHI_TransformCopy(vtkAlgorithm *a, int port, vtkTransform *trans);

  static vtkSmartPointer<vtkImplicitFunction> CreateImplicitPlane(
    vtkSmartPointer<vtkVariantArray> row);
  static vtkSmartPointer<vtkImplicitFunction> CreateImplicitBox(
    vtkSmartPointer<vtkVariantArray> row, double bounds[6], double scalefactor);
//ETX

 protected:
   vtkProbeFileReader();
  ~vtkProbeFileReader();

  //
  //
  // Read the data
  //
  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

//BTX
  vtkSmartPointer<vtkTable>    ReadProbeTableData();
  vtkSmartPointer<vtkDataSet>  MakeOneRectProbe(vtkSmartPointer<vtkVariantArray> row);
  vtkSmartPointer<vtkPolyData> MakeOneDiskProbe(vtkSmartPointer<vtkVariantArray> row);
  vtkSmartPointer<vtkPolyData> MakeOneLineProbe(vtkSmartPointer<vtkVariantArray> row);
  vtkSmartPointer<vtkPolyData> MakeOnePointProbe(vtkSmartPointer<vtkVariantArray> row);
//ETX

  //
  char *ProbeFileName; 
  int   numProbes;
  //** Resolution along x direction
  int   Resolution1;
  //** Resolution along y direction, or R
  int   Resolution2;

 private:
  vtkProbeFileReader(const vtkProbeFileReader&); // not implemented
  void operator=(const vtkProbeFileReader&); // not implemented
};

#endif
