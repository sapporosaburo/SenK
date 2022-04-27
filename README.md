# SenK

SenK is a high-performance linear solver library in C++. The most important feature of SenK is its ILUB preconditioner. 

## ILUB precondtioner

The ILUB preconditioner is one of the ILU preconditioners for the Krylov subspace linear solver. In ILU preconditioning, how to control fill-ins is one of the most significant aspects. If any fill-ins are not permitted, the resultant preconditioner may only slightly improve the difficulty of the problem, but its computational complexity is relatively low. Conversely, if more fill-ins are allowed, the degree of the difficulty may improve more, but the computational effort will also be higher.

In ILUB preconditioning, a rectangle/square block-wise fill-in permission strategy applies. As a result, we can effortlessly and effectively use SIMD operations in forward and backward substitution of the ILUB preconditioner.

## Sample

A sample source code to execute the ILUB preconditioned GMRES(m) solver in shown in `sample.cpp` .

```c++
#include "senk.hpp"
int main(int argc, char **argv) {
  double *val;
  int *cind;
  int *rptr;
  // Get a coefficinet matrix by reading a MatrixMarket file.
  // The coefficinet matrix is stored in CSR format.
  double *tval;
  int *tcind;
  int *trptr;
  // Duplicate the coefficinet matrix for construction of the preconditioner.
  senk::matrix::Duplicate(
    val, cind, rptr, &tval, &tcind, &trptr, N);
  double *bval;
  int *bcind;
  int *brptr;
  const int bnl = 2;
  const int bnw = 1;
  // The block size is bnl * bnw.
  senk::matrix::Csr2Bcsr<bnl, bnw>(
    tval, tcind, trptr, &bval, &bcind, &brptr, N);
  free(tval);
  free(tcind);
  free(trptr);
  senk::matrix::Bcsr2Csr<bnl, bnw>(bval, bcind, brptr, &tval, &tcind, &trptr, N);  
  // ILU(0) factorization with the block structure allocated above.
  senk::matrix::Ilu0(tval, tcind, trptr, N);
  // Execution of the ILUB-GMRES(m)
  int m = 50;
  int max = 6000;
	senk::solver::IlubGmresm<double, bnl, bnw>(
  	val, cind, rptr,
  	blval, blcind, blrptr, buval, bucind, burptr,
  	b, x, nrm_b, max/m, m, N, epsilon
  );
}
```

