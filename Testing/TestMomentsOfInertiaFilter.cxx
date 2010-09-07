#include <cstdlib>
#include <cmath>
#include "vtkMomentsOfInertiaFilter.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkTensor.h"
#include "vtkFloatArray.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New()

const double eps = 1e-4;
const int ndim  = 3;

bool isEqual(const double x, const double y) {
   return std::abs(x - y) <= eps * std::abs(x);
} 

bool isEqual(vtkSmartPointer<vtkTensor> A, vtkSmartPointer<vtkTensor> B) {
  /// use max matrix norm 
  double norm = 0.0;
  for (int i=0; i<ndim; i++) {
    for (int j=0; j<ndim; j++) {
      const double delta = std::abs(A->GetComponent(i, j) - B->GetComponent(i, j));
      if (delta>norm) {
	norm = delta;
      }
    }
  }
  return std::abs(norm) <= eps;
}

vtkSmartPointer<vtkTensor> DoubleToTensor(double B[ndim][ndim])  {
  VTK_CREATE(vtkTensor, A);
  for (int i=0; i<ndim; i++) {
    for (int j=0; j<ndim; j++) {
      A->SetComponent(i, j, B[i][j]);
    }
  }  
  return A;
}

vtkSmartPointer<vtkTensor> FillTensor(const double a00, const double a01, const double a02,
				      const double a10, const double a11, const double a12,
				      const double a20, const double a21, const double a22)  {
  VTK_CREATE(vtkTensor, A);
  A->SetComponent(0, 0, a00);
  A->SetComponent(0, 1, a01);
  A->SetComponent(0, 2, a02);
  A->SetComponent(1, 0, a10);
  A->SetComponent(1, 1, a11);
  A->SetComponent(1, 2, a12);
  A->SetComponent(2, 0, a20);
  A->SetComponent(2, 1, a21);
  A->SetComponent(2, 2, a22);
  return A;
}

vtkSmartPointer<vtkTensor> getInertiaTensor(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkFloatArray> dataArray) {
  VTK_CREATE(vtkMomentsOfInertiaFilter, vtkmi);
  VTK_CREATE(vtkPolyData, vpd);
  vpd->SetPoints(points);
  const vtkstd::string arrayname("mass");
  dataArray->SetName(arrayname.c_str());
  vpd->GetPointData()->AddArray(dataArray);
  double ceneterPoint[ndim] = {0.0, 0.0, 0.0};
  double  inertiaTensor [ndim][ndim];
  vtkmi->ComputeInertiaTensor(vpd, arrayname, ceneterPoint, inertiaTensor);
  return DoubleToTensor(inertiaTensor);
}

int main() {
  // create a smart pointer to the vtkCenterOfMassFilter
  VTK_CREATE(vtkPoints, points);
  VTK_CREATE(vtkFloatArray, massArray);

  // set precision for tensor output
  std::ostream& mycerr(std::cerr);
  mycerr.precision(2);
  mycerr.setf(ios::fixed);


  // =========================================================================
  // one point on Y axis
  // =========================================================================
  points->InsertNextPoint(0.0, 1.0, 0.0);
  massArray->InsertNextValue(1.0);
  vtkSmartPointer<vtkTensor>  result =  getInertiaTensor(points, massArray);
  // reference inertia tensor
  vtkSmartPointer<vtkTensor> ref = FillTensor(1.0, 0.0, 0.0,
					      0.0, 0.0, 0.0,
					      0.0, 0.0, 1.0);
  if (!isEqual(result, ref)) {
    std::cerr << "Configuration number 1 faild\n"; 
    std::cerr << "ComputeInertiaTensor returns " ;
    mycerr << *result;
    std::cerr << "Reference tensor is ";
    mycerr << *ref;
    return EXIT_FAILURE;
  }
  points->Initialize();
  massArray->Initialize();

  // =========================================================================
  // one point on X axis
  // =========================================================================
  points->InsertNextPoint(1.0, 0.0, 0.0);
  massArray->InsertNextValue(1.0);
  result =  getInertiaTensor(points, massArray);
  // reference inertia tensor
  ref = FillTensor(0.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
  if (!isEqual(result, ref)) {
    std::cerr << "Configuration number 2 faild\n"; 
    std::cerr << "ComputeInertiaTensor returns: " ;
    mycerr << *result;
    std::cerr << "Reference tensor is ";
    mycerr << *ref;
    return EXIT_FAILURE;
  }
  points->Initialize();
  massArray->Initialize();

  // =========================================================================
  // one point on Z axis
  // =========================================================================
  points->InsertNextPoint(0.0, 0.0, 1.0);
  massArray->InsertNextValue(1.0);
  result =  getInertiaTensor(points, massArray);
  // reference inertia tensor
  ref = FillTensor(1.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 0.0);
  if (!isEqual(result, ref)) {
    std::cerr << "Configuration number 3 faild\n"; 
    std::cerr << "ComputeInertiaTensor returns: " ;
    mycerr << *result;
    std::cerr << "Reference tensor is ";
    mycerr << *ref;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
