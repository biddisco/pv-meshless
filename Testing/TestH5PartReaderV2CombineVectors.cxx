// SPDX-FileCopyrightText: Copyright (C) CSCS - Swiss National Supercomputing
// Centre SPDX-License-Identifier: LicenseRef-CSCS

#include "vtkH5PartReaderV2.h"

#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkMPI.h"
#include "vtkMPIController.h"
#include "vtkMultiProcessController.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTestUtilities.h"
#include "vtkTesting.h"

#include <hdf5.h>
#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

static const int NumPoints = 10;

//----------------------------------------------------------------------------
void WriteDoubleDataset(hid_t parent, const char *name, const double *data,
                        hsize_t n) {
  hsize_t dims = n;
  hid_t dataspace = H5Screate_simple(1, &dims, nullptr);
  hid_t dataset = H5Dcreate2(parent, name, H5T_NATIVE_DOUBLE, dataspace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
}

//----------------------------------------------------------------------------
// Create a test file with:
//   x, y, z               - bare coordinate arrays
//   E_x, E_y, E_z         - vector with Cartesian suffix
//   B_i, B_j, B_k         - vector with ijk suffix
//   force_u, force_v, force_w - vector with uvw suffix
//   velocity_0, velocity_1, velocity_2 - vector with numeric suffix
//   temperature           - standalone scalar
void CreateTestFile(const std::string &filename) {
  hid_t file =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    std::cerr << "Failed to create test file " << filename << std::endl;
    return;
  }

  hid_t step =
      H5Gcreate2(file, "Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  std::vector<double> x(NumPoints), y(NumPoints), z(NumPoints);
  std::vector<double> Ex(NumPoints), Ey(NumPoints), Ez(NumPoints);
  std::vector<double> Bi(NumPoints), Bj(NumPoints), Bk(NumPoints);
  std::vector<double> fu(NumPoints), fv(NumPoints), fw(NumPoints);
  std::vector<double> v0(NumPoints), v1(NumPoints), v2(NumPoints);
  std::vector<double> temp(NumPoints);

  for (int i = 0; i < NumPoints; ++i) {
    x[i] = i * 1.0;
    y[i] = i * 2.0;
    z[i] = i * 3.0;
    Ex[i] = i * 1.0;
    Ey[i] = i * 2.0;
    Ez[i] = i * 3.0;
    Bi[i] = i * 4.0;
    Bj[i] = i * 5.0;
    Bk[i] = i * 6.0;
    fu[i] = i * 7.0;
    fv[i] = i * 8.0;
    fw[i] = i * 9.0;
    v0[i] = i * 10.0;
    v1[i] = i * 11.0;
    v2[i] = i * 12.0;
    temp[i] = i * 100.0;
  }

  WriteDoubleDataset(step, "x", x.data(), NumPoints);
  WriteDoubleDataset(step, "y", y.data(), NumPoints);
  WriteDoubleDataset(step, "z", z.data(), NumPoints);
  WriteDoubleDataset(step, "E_x", Ex.data(), NumPoints);
  WriteDoubleDataset(step, "E_y", Ey.data(), NumPoints);
  WriteDoubleDataset(step, "E_z", Ez.data(), NumPoints);
  WriteDoubleDataset(step, "B_i", Bi.data(), NumPoints);
  WriteDoubleDataset(step, "B_j", Bj.data(), NumPoints);
  WriteDoubleDataset(step, "B_k", Bk.data(), NumPoints);
  WriteDoubleDataset(step, "force_u", fu.data(), NumPoints);
  WriteDoubleDataset(step, "force_v", fv.data(), NumPoints);
  WriteDoubleDataset(step, "force_w", fw.data(), NumPoints);
  WriteDoubleDataset(step, "velocity_0", v0.data(), NumPoints);
  WriteDoubleDataset(step, "velocity_1", v1.data(), NumPoints);
  WriteDoubleDataset(step, "velocity_2", v2.data(), NumPoints);
  WriteDoubleDataset(step, "temperature", temp.data(), NumPoints);

  H5Gclose(step);
  H5Fclose(file);
}

//----------------------------------------------------------------------------
vtkPolyData *UpdateReader(vtkH5PartReaderV2 *reader) {
  vtkStreamingDemandDrivenPipeline *readerSDDP =
      vtkStreamingDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
  readerSDDP->UpdateDataObject();
  readerSDDP->UpdateInformation();
  vtkInformation *outInfo = readerSDDP->GetOutputInformation(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), 0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), 1);
  outInfo->Set(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 0);
  readerSDDP->Update();
  return reader->GetOutput();
}

//----------------------------------------------------------------------------
int CheckArrayExists(vtkPolyData *output, const char *name, bool shouldExist) {
  vtkDataArray *arr = output->GetPointData()->GetArray(name);
  if (shouldExist && !arr) {
    std::cerr << "ERROR: expected array '" << name << "' but it was not found."
              << std::endl;
    return 1;
  }
  if (!shouldExist && arr) {
    std::cerr << "ERROR: array '" << name
              << "' should not exist but was found with "
              << arr->GetNumberOfComponents() << " components." << std::endl;
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
int CheckVector(vtkPolyData *output, const char *name, int expectedComponents,
                const std::vector<double> &expected) {
  vtkDataArray *arr = output->GetPointData()->GetArray(name);
  if (!arr) {
    std::cerr << "ERROR: vector '" << name << "' not found." << std::endl;
    return 1;
  }
  if (arr->GetNumberOfComponents() != expectedComponents) {
    std::cerr << "ERROR: vector '" << name << "' has "
              << arr->GetNumberOfComponents() << " components, expected "
              << expectedComponents << std::endl;
    return 1;
  }
  if (arr->GetNumberOfTuples() != NumPoints) {
    std::cerr << "ERROR: vector '" << name << "' has "
              << arr->GetNumberOfTuples() << " tuples, expected " << NumPoints
              << std::endl;
    return 1;
  }
  for (int i = 0; i < NumPoints; ++i) {
    for (int c = 0; c < expectedComponents; ++c) {
      double val = arr->GetComponent(i, c);
      double exp = expected[i * expectedComponents + c];
      if (std::abs(val - exp) > 1e-9) {
        std::cerr << "ERROR: vector '" << name << "' tuple " << i
                  << " component " << c << " = " << val << ", expected " << exp
                  << std::endl;
        return 1;
      }
    }
  }
  return 0;
}

//----------------------------------------------------------------------------
int CheckScalar(vtkPolyData *output, const char *name,
                const std::vector<double> &expected) {
  vtkDataArray *arr = output->GetPointData()->GetArray(name);
  if (!arr) {
    std::cerr << "ERROR: scalar '" << name << "' not found." << std::endl;
    return 1;
  }
  if (arr->GetNumberOfComponents() != 1) {
    std::cerr << "ERROR: scalar '" << name << "' has "
              << arr->GetNumberOfComponents() << " components, expected 1"
              << std::endl;
    return 1;
  }
  if (arr->GetNumberOfTuples() != NumPoints) {
    std::cerr << "ERROR: scalar '" << name << "' has "
              << arr->GetNumberOfTuples() << " tuples, expected " << NumPoints
              << std::endl;
    return 1;
  }
  for (int i = 0; i < NumPoints; ++i) {
    double val = arr->GetComponent(i, 0);
    double exp = expected[i];
    if (std::abs(val - exp) > 1e-9) {
      std::cerr << "ERROR: scalar '" << name << "' tuple " << i << " = " << val
                << ", expected " << exp << std::endl;
      return 1;
    }
  }
  return 0;
}

//----------------------------------------------------------------------------
int CheckCoords(vtkPolyData *output, const std::vector<double> &x,
                const std::vector<double> &y, const std::vector<double> &z) {
  vtkPoints *points = output->GetPoints();
  if (!points) {
    std::cerr << "ERROR: no points in output." << std::endl;
    return 1;
  }
  if (points->GetNumberOfPoints() != NumPoints) {
    std::cerr << "ERROR: " << points->GetNumberOfPoints()
              << " points, expected " << NumPoints << std::endl;
    return 1;
  }
  for (int i = 0; i < NumPoints; ++i) {
    double p[3];
    points->GetPoint(i, p);
    if (std::abs(p[0] - x[i]) > 1e-9 || std::abs(p[1] - y[i]) > 1e-9 ||
        std::abs(p[2] - z[i]) > 1e-9) {
      std::cerr << "ERROR: point " << i << " = (" << p[0] << ", " << p[1]
                << ", " << p[2] << "), expected (" << x[i] << ", " << y[i]
                << ", " << z[i] << ")" << std::endl;
      return 1;
    }
  }
  return 0;
}

//----------------------------------------------------------------------------
std::vector<double> ExpectedVector(double scale0, double scale1,
                                   double scale2) {
  std::vector<double> v(NumPoints * 3);
  for (int i = 0; i < NumPoints; ++i) {
    v[i * 3 + 0] = i * scale0;
    v[i * 3 + 1] = i * scale1;
    v[i * 3 + 2] = i * scale2;
  }
  return v;
}

//----------------------------------------------------------------------------
std::vector<double> ExpectedScalar(double scale) {
  std::vector<double> v(NumPoints);
  for (int i = 0; i < NumPoints; ++i) {
    v[i] = i * scale;
  }
  return v;
}

//----------------------------------------------------------------------------
std::vector<double> ExpectedMagnitude(double scale0, double scale1,
                                      double scale2) {
  std::vector<double> v(NumPoints);
  for (int i = 0; i < NumPoints; ++i) {
    v[i] = std::sqrt(i * i *
                     (scale0 * scale0 + scale1 * scale1 + scale2 * scale2));
  }
  return v;
}

} // namespace

//----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
  int initialized = 0;
  if (MPI_Initialized(&initialized) == MPI_SUCCESS && !initialized) {
    MPI_Init(&argc, &argv);
  }

  vtkSmartPointer<vtkMPIController> vtkController =
      vtkSmartPointer<vtkMPIController>::New();
  vtkController->Initialize(&argc, &argv, 1);
  vtkMultiProcessController::SetGlobalController(vtkController);

  int rank = vtkController->GetLocalProcessId();
  (void)rank;

  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c = 1; c < argc; ++c) {
    test->AddArgument(argv[c]);
  }

  char *tempDir = vtkTestUtilities::GetArgOrEnvOrDefault("-T", argc, argv,
                                                         "VTK_TEMP_DIR", ".");
  std::string filename =
      std::string(tempDir) + "/H5PartReaderV2CombineVectors.h5part";
  delete[] tempDir;

  if (vtkController->GetLocalProcessId() == 0) {
    std::remove(filename.c_str());
    CreateTestFile(filename);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int localFail = 0;

  // ------------------------------------------------------------------
  // Test 1: CombineVectorComponents on (default), no magnitude
  // ------------------------------------------------------------------
  {
    vtkSmartPointer<vtkH5PartReaderV2> reader =
        vtkSmartPointer<vtkH5PartReaderV2>::New();
    reader->SetFileName(const_cast<char *>(filename.c_str()));
    vtkPolyData *output = UpdateReader(reader);

    // Coords auto-detected from bare x, y, z
    localFail += CheckCoords(output, ExpectedScalar(1.0), ExpectedScalar(2.0),
                             ExpectedScalar(3.0));

    // Combined vectors
    localFail += CheckVector(output, "E", 3, ExpectedVector(1.0, 2.0, 3.0));
    localFail += CheckVector(output, "B", 3, ExpectedVector(4.0, 5.0, 6.0));
    localFail += CheckVector(output, "force", 3, ExpectedVector(7.0, 8.0, 9.0));
    localFail +=
        CheckVector(output, "velocity", 3, ExpectedVector(10.0, 11.0, 12.0));

    // Standalone scalar
    localFail += CheckScalar(output, "temperature", ExpectedScalar(100.0));

    // Individual component arrays should NOT be present
    const char *notExpected[] = {
        "E_x",        "E_y",        "E_z",     "B_i",     "B_j",
        "B_k",        "force_u",    "force_v", "force_w", "velocity_0",
        "velocity_1", "velocity_2", "x",       "y",       "z"};
    for (const char *name : notExpected) {
      localFail += CheckArrayExists(output, name, false);
    }

    // No magnitude arrays by default
    localFail += CheckArrayExists(output, "E_magnitude", false);

    std::cout << "Test 1 (combine on, no magnitude): "
              << (localFail ? "FAIL" : "PASS") << std::endl;
  }

  // ------------------------------------------------------------------
  // Test 2: CombineVectorComponents on + ExportVectorComponentsMagnitude
  // ------------------------------------------------------------------
  {
    vtkSmartPointer<vtkH5PartReaderV2> reader =
        vtkSmartPointer<vtkH5PartReaderV2>::New();
    reader->SetFileName(const_cast<char *>(filename.c_str()));
    reader->SetExportVectorComponentsMagnitude(1);
    vtkPolyData *output = UpdateReader(reader);

    // Magnitude arrays for combined vectors
    localFail +=
        CheckScalar(output, "E_magnitude", ExpectedMagnitude(1.0, 2.0, 3.0));
    localFail +=
        CheckScalar(output, "B_magnitude", ExpectedMagnitude(4.0, 5.0, 6.0));
    localFail += CheckScalar(output, "force_magnitude",
                             ExpectedMagnitude(7.0, 8.0, 9.0));
    localFail += CheckScalar(output, "velocity_magnitude",
                             ExpectedMagnitude(10.0, 11.0, 12.0));

    // No magnitude for standalone scalar (Nc == 1)
    localFail += CheckArrayExists(output, "temperature_magnitude", false);

    std::cout << "Test 2 (combine on, magnitude on): "
              << (localFail ? "FAIL" : "PASS") << std::endl;
  }

  // ------------------------------------------------------------------
  // Test 3: CombineVectorComponents off — all arrays as separate scalars
  // ------------------------------------------------------------------
  {
    vtkSmartPointer<vtkH5PartReaderV2> reader =
        vtkSmartPointer<vtkH5PartReaderV2>::New();
    reader->SetFileName(const_cast<char *>(filename.c_str()));
    reader->SetCombineVectorComponents(0);
    vtkPolyData *output = UpdateReader(reader);

    // No combined vectors
    localFail += CheckArrayExists(output, "E", false);
    localFail += CheckArrayExists(output, "B", false);
    localFail += CheckArrayExists(output, "force", false);
    localFail += CheckArrayExists(output, "velocity", false);

    // Individual component arrays should be present as 1-component scalars
    localFail += CheckScalar(output, "E_x", ExpectedScalar(1.0));
    localFail += CheckScalar(output, "E_y", ExpectedScalar(2.0));
    localFail += CheckScalar(output, "E_z", ExpectedScalar(3.0));
    localFail += CheckScalar(output, "temperature", ExpectedScalar(100.0));

    std::cout << "Test 3 (combine off): " << (localFail ? "FAIL" : "PASS")
              << std::endl;
  }

  if (vtkController->GetLocalProcessId() == 0) {
    std::remove(filename.c_str());
  }

  int globalFail = 0;
  MPI_Allreduce(&localFail, &globalFail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  vtkController->Finalize();
  return globalFail;
}
