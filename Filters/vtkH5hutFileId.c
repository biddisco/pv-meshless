#include "../h5core/private/h5_types.h"
#include "H5hut.h"

hid_t vtkH5hutGetHDF5FileId(h5_file_t f) { return ((h5_file_p)f)->file; }

hid_t vtkH5hutGetHDF5IterationGroupId(h5_file_t f) {
  return ((h5_file_p)f)->iteration_gid;
}
