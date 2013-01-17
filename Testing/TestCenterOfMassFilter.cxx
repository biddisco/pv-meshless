#include <cstdlib>
#include <cmath>
#include "vtkCenterOfMassFilter.h"
#include "vtkPoints.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkNew.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New()

const double eps = 1e-4;
inline bool isEqual(const double x, const double y) {
   return std::abs(x - y) <= eps * std::abs(x);
} 

int main() {
  // create a smart pointer to the vtkCenterOfMassFilter
  VTK_CREATE(vtkCenterOfMassFilter, vtkcom);
  // =========================================================================
  // one point data set
  // =========================================================================
  vtkNew<vtkPoints> points;
  points->InsertNextPoint(1.0, 2.0, 3.0);

  vtkNew<vtkDoubleArray> mass;
  mass->InsertNextValue(2.0);

  double com[3];
  vtkcom->ComputeCenterOfMass(points.GetPointer(), mass.GetPointer(), com);

  if (!isEqual(com[0], 1.0)) {
    std::cerr << "com[0] is wrong for one point data set\n";
    return EXIT_FAILURE;
  }
  if (!isEqual(com[1], 2.0)) {
    std::cerr << "com[1] is wrong for one point data set\n";
    return EXIT_FAILURE;
  }
  if (!isEqual(com[2], 3.0)) {
    std::cerr << "com[2] is wrong for one point data set\n";
    return EXIT_FAILURE;
  }
  mass->Initialize();
  points->Initialize();
  // =========================================================================
  // two points data set
  // =========================================================================
  points->InsertNextPoint(0.5, 1.0, 1.5);
  points->InsertNextPoint(-1.0, -2.0, -3.0);
  mass->InsertNextValue(2.0);
  mass->InsertNextValue(1.0);
  vtkcom->ComputeCenterOfMass(points.GetPointer(), mass.GetPointer(), com);
  if (!isEqual(com[0], 0.0)) {
    std::cerr << "com[0] is wrong for two points data set\n";
    return EXIT_FAILURE;
  }
  if (!isEqual(com[1], 0.0)) {
    std::cerr << "com[1] is wrong for two points data set\n";
    return EXIT_FAILURE;
  }
  if (!isEqual(com[2], 0.0)) {
    std::cerr << "com[2] is wrong for two points data set\n";
    return EXIT_FAILURE;
  }
  mass->Initialize();
  points->Initialize();
  return EXIT_SUCCESS;
}
