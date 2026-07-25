/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkH5hutHelper.h

  Copyright (C) CSCS - Swiss National Supercomputing Centre.

=========================================================================*/
#ifndef __vtkH5hutHelper_h
#define __vtkH5hutHelper_h

#include "H5hut.h"

#include "vtkCellType.h"

extern "C" {
hid_t vtkH5hutGetHDF5FileId(h5_file_t f);
hid_t vtkH5hutGetHDF5IterationGroupId(h5_file_t f);
}
#include "vtkCharArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkUnsignedShortArray.h"

#include "vtkMPI.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"
#include "vtkMultiProcessController.h"

#include <string>
#include <vector>
#include <vtksys/SystemTools.hxx>

//----------------------------------------------------------------------------
// Simple random number generator used by the reader.
class Random {
public:
  unsigned int __seed;
  Random(int seed) { __seed = seed; }
  unsigned int getseed() { return __seed; }
  void setseed(int seed) { __seed = seed; }
  double nextNumber() {
    __seed = (__seed * 9301 + 49297) % 233280;
    return __seed / 233280.0;
  }
  int nextNumberInt() {
    __seed = (__seed * 9301 + 49297) % 233280;
    return __seed;
  }
};

//----------------------------------------------------------------------------
// Map a H5hut dataset type constant to a VTK data type.
inline int H5hutTypeToVTKType(h5_int64_t type) {
  if (type == H5_FLOAT64_T) {
    return VTK_DOUBLE;
  } else if (type == H5_FLOAT32_T) {
    return VTK_FLOAT;
  } else if (type == H5_INT8_T) {
    return VTK_CHAR;
  } else if (type == H5_UINT8_T) {
    return VTK_UNSIGNED_CHAR;
  } else if (type == H5_INT16_T) {
    return VTK_SHORT;
  } else if (type == H5_UINT16_T) {
    return VTK_UNSIGNED_SHORT;
  } else if (type == H5_INT32_T) {
    return VTK_INT;
  } else if (type == H5_UINT32_T) {
    return VTK_UNSIGNED_INT;
  } else if (type == H5_INT64_T) {
    return VTK_LONG_LONG;
  } else if (type == H5_UINT64_T) {
    return VTK_UNSIGNED_LONG_LONG;
  }
  return VTK_VOID;
}

//----------------------------------------------------------------------------
// Map a VTK data type to a H5hut type constant.
inline h5_int64_t VTKTypeToH5hutType(int vtkType) {
  switch (vtkType) {
  case VTK_CHAR:
    return H5_INT8_T;
  case VTK_SIGNED_CHAR:
    return H5_INT8_T;
  case VTK_UNSIGNED_CHAR:
    return H5_UINT8_T;
  case VTK_SHORT:
    return H5_INT16_T;
  case VTK_UNSIGNED_SHORT:
    return H5_UINT16_T;
  case VTK_INT:
    return H5_INT32_T;
  case VTK_UNSIGNED_INT:
    return H5_UINT32_T;
  case VTK_LONG:
    return (sizeof(long) == 8 ? H5_INT64_T : H5_INT32_T);
  case VTK_UNSIGNED_LONG:
    return (sizeof(unsigned long) == 8 ? H5_UINT64_T : H5_UINT32_T);
  case VTK_LONG_LONG:
    return H5_INT64_T;
  case VTK_UNSIGNED_LONG_LONG:
    return H5_UINT64_T;
  case VTK_FLOAT:
    return H5_FLOAT32_T;
  case VTK_DOUBLE:
    return H5_FLOAT64_T;
  case VTK_ID_TYPE:
    return (sizeof(vtkIdType) == 8 ? H5_INT64_T : H5_INT32_T);
  default:
    return -1;
  }
}

//----------------------------------------------------------------------------
// Return a vector with the names of all datasets in the current step.
inline std::vector<std::string> H5hutDatasetNames(h5_file_t f) {
  std::vector<std::string> names;
  h5_ssize_t nds = H5PartGetNumDatasets(f);
  if (nds > 0) {
    names.resize(static_cast<size_t>(nds));
    char name[512];
    for (h5_ssize_t i = 0; i < nds; ++i) {
      H5PartGetDatasetName(f, i, name, 512);
      names[i] = std::string(name);
    }
  }
  return names;
}

//----------------------------------------------------------------------------
// Open an H5hut file. In parallel this sets up an MPI communicator property.
inline h5_file_t H5hutOpenFile(const char *fname, h5_int64_t mode,
                               vtkMultiProcessController *controller) {
  h5_file_t f = 0;
  std::string fullPath = vtksys::SystemTools::CollapseFullPath(fname);
#ifdef PARALLEL_IO
  MPI_Comm comm = MPI_COMM_WORLD;
  if (controller) {
    vtkMPICommunicator *vtkComm =
        vtkMPICommunicator::SafeDownCast(controller->GetCommunicator());
    if (vtkComm && vtkComm->GetMPIComm()) {
      comm = *vtkComm->GetMPIComm()->GetHandle();
    }
  }
  h5_prop_t prop = H5CreateFileProp();
  H5SetPropFileMPIOCollective(prop, &comm);
  f = H5OpenFile(fullPath.c_str(), mode, prop);
  H5CloseProp(prop);
#else
  f = H5OpenFile(fullPath.c_str(), mode, H5_PROP_DEFAULT);
#endif
  if (f == static_cast<h5_file_t>(H5_ERR)) {
    f = 0;
  }
  return f;
}

//----------------------------------------------------------------------------
// Close an H5hut file.
inline void H5hutCloseFile(h5_file_t &f) {
  if (f) {
    H5CloseFile(f);
    f = 0;
  }
}

//----------------------------------------------------------------------------
// Write a 1D data array into the current step under the given name.
inline h5_err_t H5hutWriteDataArray(h5_file_t f, const char *name,
                                    vtkDataArray *data) {
  if (!data) {
    return H5_FAILURE;
  }
  h5_int64_t type = VTKTypeToH5hutType(data->GetDataType());
  if (type < 0) {
    return H5_FAILURE;
  }
  return h5u_write_dataset(f, name, data->GetVoidPointer(0),
                           static_cast<h5_types_t>(type));
}

//----------------------------------------------------------------------------
// Map a VTK data type to the corresponding HDF5 native type.
inline hid_t VTKTypeToHDF5NativeType(int vtkType) {
  switch (vtkType) {
  case VTK_CHAR:
  case VTK_SIGNED_CHAR:
    return H5T_NATIVE_CHAR;
  case VTK_UNSIGNED_CHAR:
    return H5T_NATIVE_UCHAR;
  case VTK_SHORT:
    return H5T_NATIVE_SHORT;
  case VTK_UNSIGNED_SHORT:
    return H5T_NATIVE_USHORT;
  case VTK_INT:
    return H5T_NATIVE_INT;
  case VTK_UNSIGNED_INT:
    return H5T_NATIVE_UINT;
  case VTK_LONG:
    return H5T_NATIVE_LONG;
  case VTK_UNSIGNED_LONG:
    return H5T_NATIVE_ULONG;
  case VTK_LONG_LONG:
    return H5T_NATIVE_LLONG;
  case VTK_UNSIGNED_LONG_LONG:
    return H5T_NATIVE_ULLONG;
  case VTK_FLOAT:
    return H5T_NATIVE_FLOAT;
  case VTK_DOUBLE:
    return H5T_NATIVE_DOUBLE;
  case VTK_ID_TYPE:
    return (sizeof(vtkIdType) == 8 ? H5T_NATIVE_INT64 : H5T_NATIVE_INT32);
  default:
    return H5I_INVALID_HID;
  }
}

//----------------------------------------------------------------------------
// Read a 1D data array from the current step/view into the supplied array.
inline h5_err_t H5hutReadDataArray(h5_file_t f, const char *name,
                                   vtkDataArray *data) {
  if (!data) {
    return H5_FAILURE;
  }
  h5_int64_t type = VTKTypeToH5hutType(data->GetDataType());
  if (type < 0) {
    return H5_FAILURE;
  }
  return h5u_read_dataset(f, name, data->GetVoidPointer(0),
                          static_cast<h5_types_t>(type));
}

//----------------------------------------------------------------------------
// Read a 1D data array using an HDF5 strided hyperslab selection.
// This bypasses H5PartSetViewIndices(), which constructs a slow HDF5
// point selection for strided access.
inline h5_err_t H5hutReadDataArrayStrided(h5_file_t f, const char *name,
                                          vtkDataArray *data,
                                          h5_size_t start, h5_size_t stride,
                                          h5_size_t count) {
  if (!f || !name || !data || count == 0) {
    return H5_FAILURE;
  }

  hid_t mem_type = VTKTypeToHDF5NativeType(data->GetDataType());
  if (mem_type == H5I_INVALID_HID) {
    return H5_FAILURE;
  }

  // Use H5hut's current iteration group directly so we do not have to
  // know the step-name formatting used in the file.
  hid_t groupId = vtkH5hutGetHDF5IterationGroupId(f);
  if (groupId < 0) {
    return H5_FAILURE;
  }

  hid_t datasetId = H5Dopen2(groupId, name, H5P_DEFAULT);
  if (datasetId < 0) {
    return H5_FAILURE;
  }

  hid_t fileSpaceId = H5Dget_space(datasetId);
  if (fileSpaceId < 0) {
    H5Dclose(datasetId);
    return H5_FAILURE;
  }

  hsize_t hstart = static_cast<hsize_t>(start);
  hsize_t hstride = static_cast<hsize_t>(stride);
  hsize_t hcount = static_cast<hsize_t>(count);
  herr_t status = H5Sselect_hyperslab(fileSpaceId, H5S_SELECT_SET, &hstart,
                                      &hstride, &hcount, nullptr);
  if (status < 0) {
    H5Sclose(fileSpaceId);
    H5Dclose(datasetId);
    return H5_FAILURE;
  }

  hsize_t memDim = hcount;
  hid_t memSpaceId = H5Screate_simple(1, &memDim, nullptr);
  if (memSpaceId < 0) {
    H5Sclose(fileSpaceId);
    H5Dclose(datasetId);
    return H5_FAILURE;
  }

  status = H5Dread(datasetId, mem_type, memSpaceId, fileSpaceId, H5P_DEFAULT,
                   data->GetVoidPointer(0));

  H5Sclose(memSpaceId);
  H5Sclose(fileSpaceId);
  H5Dclose(datasetId);

  return status >= 0 ? H5_SUCCESS : H5_FAILURE;
}

#endif
