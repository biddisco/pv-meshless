#include <iostream>
#include <vector>
#include "KernelBox.h"
#include "KernelCusp.h"
#include "KernelGaussian.h"
#include "KernelQuadratic.h"
#include "KernelSpline3rdOrder.h"
#include "KernelSpline5thOrder.h"
#include "KernelWendland.h"


/// Simpson integradion rule
double SimpsonIntegration(const std::vector<double> &y, 
			  const std::vector<double> &x) {
  const int n = x.size()-1;
  const double aux = 2.*(x[n]-x[0])/(3.*static_cast<double>(n));
  double val = 0.5*(y[0]*x[0]+y[n]*x[n]);
  for(int i=2;i<n;i+=2){
    val+=2*y[i-1]+y[i];
  }
  val+=2*y[n-1];
  return aux*val;
}

/// Integrade a kernel in 3D
double IntegralOfKernel(const Kernel& ker) {
  const double supportlength = ker.maxDistance();
  const int npoint = 5000;
  std::vector<double> x;
  std::vector<double> y;
  for (int i = 0; i<npoint; i++) {
    const double r = static_cast<double>(i) * supportlength / static_cast<double>(npoint);
    x.push_back(r);
    y.push_back(ker.w(r)*r);
  }
  return 2.0*M_PI*SimpsonIntegration(y, x);
}


int main() {
  const double eps = 1e-4;
  const int ndim = 2;
  const double smoothinglength = 1.0/3.0;
  double s;
  //
  s  = IntegralOfKernel(KernelBox(ndim, smoothinglength));
  if ( fabs(s - 1.0) > eps) {
    return EXIT_FAILURE;
  }
  s  = IntegralOfKernel(KernelCusp(ndim, smoothinglength));
  if ( fabs(s - 1.0) > eps) {
    return EXIT_FAILURE;
  }
  s  = IntegralOfKernel(KernelGaussian(ndim, smoothinglength));
  if ( fabs(s - 1.0) > eps) {
    return EXIT_FAILURE;
  }
  s  = IntegralOfKernel(KernelQuadratic(ndim, smoothinglength));
  if ( fabs(s - 1.0) > eps) {
    return EXIT_FAILURE;
  }
  s  = IntegralOfKernel(KernelSpline3rdOrder(ndim, smoothinglength));
  if ( fabs(s - 1.0) > eps) {
    return EXIT_FAILURE;
  }
  s  = IntegralOfKernel(KernelSpline5thOrder(ndim, smoothinglength));
  if ( fabs(s - 1.0) > eps) {
    return EXIT_FAILURE;
  }
  s  = IntegralOfKernel(KernelWendland(ndim, smoothinglength));
  if ( fabs(s - 1.0) > eps) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
