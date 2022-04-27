# SenK

SenK is a high-performance linear solver library in C++. The most important feature of SenK is its ILUB preconditioner. 

## ILUB precondtioner

The ILUB preconditioner is one of the ILU preconditioners for the Krylov subspace linear solver. In ILU preconditioning, how to control fill-ins is one of the most significant aspects. If any fill-ins are not permitted, the resultant preconditioner may only slightly improve the difficulty of the problem, but its computational complexity is relatively low.

## Sample

A sample source code to execute ILUB precondtioned`sample.cpp` 

```c++
#include "senk.hpp"
int main(int argc, char **argv) {
  
}
```

