// SPDX-FileCopyrightText: Copyright (C) CSCS - Swiss National Supercomputing
// Centre SPDX-License-Identifier: LicenseRef-CSCS

#include "vtkH5PartReaderV2.h"

#include "vtkDataArray.h"
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

static const int NumPartitions = 4;
static const int PointsPerPartition = 100;
static const int TotalPoints = NumPartitions * PointsPerPartition;

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

void WriteInt64Dataset(hid_t parent, const char *name, const long long *data,
                       hsize_t n) {
  hsize_t dims = n;
  hid_t dataspace = H5Screate_simple(1, &dims, nullptr);
  hid_t dataset = H5Dcreate2(parent, name, H5T_NATIVE_LLONG, dataspace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  H5Sclose(dataspace);
}

void CreateTestFile(const std::string &filename) {
  hid_t file =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    std::cerr << "Failed to create test file " << filename << std::endl;
    return;
  }

  hid_t step =
      H5Gcreate2(file, "Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  std::vector<double> x(TotalPoints);
  std::vector<double> y(TotalPoints);
  std::vector<double> z(TotalPoints);
  std::vector<long long> id(TotalPoints);

  for (int p = 0; p < NumPartitions; ++p) {
    for (int i = 0; i < PointsPerPartition; ++i) {
      int idx = p * PointsPerPartition + i;
      x[idx] = p + i * 0.01;
      y[idx] = i * 0.01;
      z[idx] = i * 0.01;
      id[idx] = idx;
    }
  }

  WriteDoubleDataset(step, "Coords_0", x.data(), TotalPoints);
  WriteDoubleDataset(step, "Coords_1", y.data(), TotalPoints);
  WriteDoubleDataset(step, "Coords_2", z.data(), TotalPoints);
  WriteInt64Dataset(step, "id", id.data(), TotalPoints);

  H5Gclose(step);

  hid_t partGroup =
      H5Gcreate2(file, "Partition#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  std::vector<double> box(NumPartitions * 13, 0.0);
  for (int p = 0; p < NumPartitions; ++p) {
    int off = p * 13;
    box[off + 0] = PointsPerPartition; // count
    box[off + 1] = p;                  // xmin
    box[off + 2] = 0.0;                // ymin
    box[off + 3] = 0.0;                // zmin
    box[off + 4] = p + 1.0;            // xmax
    box[off + 5] = 1.0;                // ymax
    box[off + 6] = 1.0;                // zmax
    box[off + 7] = p;                  // ghost xmin
    box[off + 8] = 0.0;
    box[off + 9] = 0.0;
    box[off + 10] = p + 1.0; // ghost xmax
    box[off + 11] = 1.0;
    box[off + 12] = 1.0;
  }

  hsize_t boxDims = NumPartitions * 13;
  WriteDoubleDataset(partGroup, "Box", box.data(), boxDims);

  H5Gclose(partGroup);
  H5Fclose(file);
}

int ExpectedPointCount(int rank, int size) {
  if (size == 1) {
    return TotalPoints;
  }
  int start = (rank * NumPartitions) / size;
  int end = ((rank + 1) * NumPartitions) / size - 1;
  return (start > end) ? 0 : (end - start + 1) * PointsPerPartition;
}

int ExpectedMinPartition(int rank, int size) {
  int start = (rank * NumPartitions) / size;
  int end = ((rank + 1) * NumPartitions) / size - 1;
  return (start > end) ? -1 : start;
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
  int size = vtkController->GetNumberOfProcesses();

  std::cout << "Rank " << rank << " starting" << std::endl;

  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c = 1; c < argc; ++c) {
    test->AddArgument(argv[c]);
  }

  char *tempDir = vtkTestUtilities::GetArgOrEnvOrDefault("-T", argc, argv,
                                                         "VTK_TEMP_DIR", ".");
  std::string filename =
      std::string(tempDir) + "/H5PartReaderV2Partitioning.h5part";
  delete[] tempDir;

  if (rank == 0) {
    std::remove(filename.c_str());
    CreateTestFile(filename);
    std::cout << "Rank 0 created " << filename << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout << "Rank " << rank << " past barrier" << std::endl;

  int localFail = 0;

  vtkSmartPointer<vtkH5PartReaderV2> reader =
      vtkSmartPointer<vtkH5PartReaderV2>::New();
  reader->SetFileName(const_cast<char *>(filename.c_str()));

  vtkStreamingDemandDrivenPipeline *readerSDDP =
      vtkStreamingDemandDrivenPipeline::SafeDownCast(reader->GetExecutive());
  readerSDDP->UpdateDataObject();
  readerSDDP->UpdateInformation();
  vtkInformation *outInfo = readerSDDP->GetOutputInformation(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), rank);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
               size);
  outInfo->Set(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 0);
  readerSDDP->Update();

  vtkPolyData *output = reader->GetOutput();
  vtkIdType localCount = output->GetNumberOfPoints();
  int expectedCount = ExpectedPointCount(rank, size);

  std::cout << "Rank " << rank << "/" << size << " read " << localCount
            << " points (expected " << expectedCount << ")" << std::endl;

  if (localCount != expectedCount) {
    std::cerr << "ERROR: Rank " << rank << " expected " << expectedCount
              << " points but got " << localCount << std::endl;
    localFail = 1;
  }

  vtkDataArray *ids = output->GetPointData()->GetArray("id");
  if (ids) {
    int minPart = ExpectedMinPartition(rank, size);
    for (vtkIdType i = 0; i < localCount; ++i) {
      long long gid = static_cast<long long>(ids->GetTuple1(i));
      int part = static_cast<int>(gid / PointsPerPartition);
      if (minPart >= 0 &&
          (part < minPart ||
           part >= minPart + expectedCount / PointsPerPartition)) {
        std::cerr << "ERROR: Rank " << rank << " got out-of-range id " << gid
                  << " (partition " << part << ")" << std::endl;
        localFail = 1;
        break;
      }
    }
  }

  long long localCountLL = localCount;
  long long totalCount = 0;
  MPI_Allreduce(&localCountLL, &totalCount, 1, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);

  if (rank == 0 && totalCount != TotalPoints) {
    std::cerr << "ERROR: Total points across ranks " << totalCount
              << " != expected " << TotalPoints << std::endl;
    localFail = 1;
  }

  if (rank == 0) {
    std::remove(filename.c_str());
  }

  int globalFail = 0;
  MPI_Allreduce(&localFail, &globalFail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  vtkController->Finalize();
  return globalFail;
}
