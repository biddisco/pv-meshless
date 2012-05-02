/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5SPHReader.h
  Revision of last commit : $Rev: 501 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2008-03-11 20:17:29 +0100 #$

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing 
  1) This copyright notice appears on all copies of source code 
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it 
  must not be reformatted such that the indentation, bracketing or 
  overall style is modified significantly. 

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
// .NAME vtkH5SPHReader - Write H5Part (HDF5) Particle files
// .SECTION Description
// vtkH5SPHReader reads compatible with H5Part : documented here
// http://amas.web.psi.ch/docs/H5Part-doc/h5part.html 

#ifndef __vtkH5SPHReader_h
#define __vtkH5SPHReader_h

#include "vtkH5PartReader.h"
#include <vtkstd/string>
#include <vtkstd/vector>
#include <hdf5.h>
#include "H5Part.h"
#include <map>

class vtkMultiProcessController;

//BTX
struct H5PartFile;
class  FileSeriesFinder;

class compound_info {
  public:
    compound_info() : 
      vtk_type(0), H5DataType(0), offset(0), size(0) {};
      compound_info(const compound_info &c) :
      vtk_type(c.vtk_type), offset(c.offset), size(c.size) 
    {
      H5DataType = H5Tcopy(c.H5DataType);
    }
    compound_info(int v, hid_t h, int o, int s) : 
      vtk_type(v), offset(o), size(s) 
    {
      H5DataType = H5Tcopy(h);
    };
    //
    ~compound_info() 
    { 
      if (H5DataType>0) H5Tclose(H5DataType); 
    }
    int   vtk_type;
    hid_t H5DataType;
    int   offset;
    int   size;
};
//ETX

class VTK_EXPORT vtkH5SPHReader : public vtkH5PartReader
{
public:
//BTX
  typedef vtkstd::vector<vtkstd::string>              stringlist;
  typedef vtkstd::pair<vtkstd::string, compound_info> CompoundType;
  typedef vtkstd::map<vtkstd::string, compound_info>  CompoundInfo;
//ETX

  static vtkH5SPHReader *New();
  vtkTypeRevisionMacro(vtkH5SPHReader,vtkH5PartReader);
  void PrintSelf(ostream& os, vtkIndent indent);   

  // Description:
  // Description:
  // Set/Get the prefix part of the hdf dataset step name
  // by default data is stored as Step#00000 the prefix
  // is "Step" and the width is 5.
  vtkGetStringMacro(StepNamePrefix);
  vtkSetStringMacro(StepNamePrefix);
  
  // Description:
  // Set/Get the width part of the hdf dataset step name
  // by default data is stored as Step#00000 the prefix
  // is "Step" and the width is 5.
  vtkGetMacro(StepNameWidth, int);
  vtkSetMacro(StepNameWidth, int);

  // Description:
  // A convenience function to set both name and width in one call
  void SetStepNameFormat(char *name, int width)
    { this->SetStepNamePrefix(name); this->SetStepNameWidth(width); }

  // In general when fetching multiple time steps from files
  // we need a file pattern, by default we use
  // "PREFIX TIME TEXT0 EXT" which defines the inital part as 
  // <verbatim>
  // prefix - usually path, plus any leadinf characters
  // time   - a numeric part indicating time step
  // text0  - usually some non changing part of the file name
  // ext    - the extension
  // and example using the above pattern might be 
  // 0000000100.dam_breaking_5853510.champs.dat.h5sph
  // </verbatim>
  vtkGetStringMacro(FileNamePattern);
  vtkSetStringMacro(FileNamePattern);

  bool HasStep(int Step);

protected:
   vtkH5SPHReader();
  ~vtkH5SPHReader();
  //
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   OpenFile();
  void  CloseFile();
//  void CopyIntoCoords(int offset, vtkDataArray *source, vtkDataArray *dest);
//BTX
  static int      ScanCompoundType(hid_t loc_id, const char *name, void *opdata);
  int             FindCompoundDataSet(hid_t group_id, const char *group_name, const hid_t type, char * const pattern);
  vtkIdType       GetNumberOfParticles();
//ETX
  //
  // Internal Variables
  //
  char         *StepNamePrefix;
  int           StepNameWidth;
  char         *FileNamePattern;
//
  //BTX
  vtkstd::string                             CompoundName;
  CompoundInfo                               CompoundData;
  int                                        CompoundSize;
  std::vector<int>                           CompoundOffset;
  FileSeriesFinder                          *Finder;
  vtkstd::string                             FileNameInternal;
  vtkstd::string                             FileNameLast;
  //ETX

private:
  vtkH5SPHReader(const vtkH5SPHReader&);  // Not implemented.
  void operator=(const vtkH5SPHReader&);  // Not implemented.
};

#endif
