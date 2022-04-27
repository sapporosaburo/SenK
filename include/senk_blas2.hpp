#ifndef SENK_BLAS2_HPP
#define SENK_BLAS2_HPP

namespace senk {

namespace blas2 {

// Upper triangular solver
template <typename T> inline
void Trsv(T *U, T *b, T *x, int n, int m)
{
    for(int i=m-1; i>=0; i--) {
        T temp = b[i];
        for(int j=m-1; j>i; j--) {
            temp -= U[j*n+i] * x[j];
        }
        x[i] = temp / U[i*n+i];
    }
}

}

}

#endif
