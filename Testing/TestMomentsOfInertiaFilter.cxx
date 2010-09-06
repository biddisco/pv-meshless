#include <cstdlib>
#include <cmath>
#include "vtkMomentsOfInertiaFilter.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New()

const double eps = 1e-4;
inline bool isEqual(const double x, const double y) {
   return std::abs(x - y) <= eps * std::abs(x);
} 

int main() {
  // create a smart pointer to the vtkCenterOfMassFilter
  VTK_CREATE(vtkMomentsOfInertiaFilter, vtkmi);
  vtkPolyData *vpd = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  points->InsertNextPoint(0.0, 0.0, 0.0);
  vpd->SetPoints(points);
  double ceneterPoint[3] = {0.0, 0.0, 0.0};
  double  inertiaTensor [3][3];
  vtkmi->ComputeInertiaTensor(vpd, vtkstd::string("mass"), ceneterPoint, inertiaTensor);

  return EXIT_SUCCESS;
}
