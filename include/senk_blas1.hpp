#ifndef SENK_BLAS1_HPP
#define SENK_BLAS1_HPP

#include <cmath>
#include <cstdio>

#ifdef __GNUC__
#include <omp.h>
#endif

#include "senk_helper.hpp"

namespace senk {

namespace blas1 {

template <typename T> inline
void Copy(T *x, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] = x[i]; }
}
template <typename T> inline
void Scal(T a, T *x, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { x[i] *= a; }
}
template <typename T> inline
void Axpy(T a, T *x, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] += a * x[i]; }
}
template <typename T> inline
void Axpby(T a, T *x, T b, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] = a * x[i] + b * y[i]; }
}
template <typename T> inline
void Axpyz(T a, T *x, T *y, T *z, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { z[i] = a * x[i] + y[i]; }
}
template <typename T> inline
T Dot(T *x, T *y, int N) {
    T res = 0;
    #pragma omp parallel for simd reduction(+: res)
    for(int i=0; i<N; i++) { res += x[i] * y[i]; }
    return res;
}
template <typename T> inline
T Nrm2(T *x, int N) {
    T res = 0;
    #pragma omp parallel for simd reduction(+: res)
    for(int i=0; i<N; i++) { res += x[i] * x[i]; }
    return std::sqrt(res);
}
template <typename T> inline
void HadProd(T *x, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] *= x[i]; }
}
template <typename T> inline
void HadDiv(T *x, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] /= x[i]; }
}
template <typename T> inline
T Ggen(T a, T b, T *c, T *s)
{
    T r;
    r = std::sqrt(a*a + b*b);
    c[0] =  a / r;
    s[0] = -b / r;
    return r;
}
template <typename T> inline
void Grot(T c, T s, T *a, T *b)
{
    T temp = a[0];
    a[0] = c * temp - s * b[0];
    b[0] = s * temp + c * b[0];
}

}

}


#endif
