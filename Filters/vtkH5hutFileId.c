#include "../h5core/private/h5_types.h"
#include "H5hut.h"

hid_t vtkH5hutGetHDF5FileId(h5_file_t f) { return ((h5_file_p)f)->file; }
