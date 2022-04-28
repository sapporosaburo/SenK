# SenK

SenK is a high-performance linear solver library in C++. 

## Sample

A sample source code to execute the ILU(0) preconditioned GMRES(m) solver in shown in `sample.cpp` .

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
  // ILU(0) factorization with the block structure allocated above.
  senk::matrix::Ilu0(tval, tcind, trptr, N);
  double *lval;
  int *lcind;
  int *lrptr;
  double *uval;
  int *ucind;
  int *urptr;
  senk::matrix::Split(
    tval, tcind, trptr,
    &lval, &lcind, &lrptr,
    &uval, &ucind, &urptr,
    nullptr, N, "L-DU", true
  );
  // Execution of the ILU-GMRES(m)
  int m = 50;
  int max = 6000;
  senk::solver::IluGmresm<double, bnl, bnw>(
    val, cind, rptr,
    lval, lcind, lrptr, uval, ucind, urptr,
    b, x, nrm_b, max/m, m, N, epsilon
  );
}
```

