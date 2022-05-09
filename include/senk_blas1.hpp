/**
 * @file senk_blas1.hpp
 * @brief Level1 BLAS-style functions are written.
 * @author Kengo Suzuki
 * @date 5/8/2021
 */
#ifndef SENK_BLAS1_HPP
#define SENK_BLAS1_HPP

#include <cmath>

#include "senk_helper.hpp"

namespace senk {
/**
 * @brief This namespace contains Level1 BLAS-style functions.
 */
namespace blas1 {
/**
 * @brief Copy x to y.
 * @tparam T The type of vectors.
 * @param x A 1D-array of size N.
 * @param y A 1D-array of size N.
 * @param N The size of vectors.
 */
template <typename T> inline
void Copy(T *x, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] = x[i]; }
}
/**
 * @brief Multiply x by a.
 * @tparam T The type of a vector.
 * @param a A scalar value.
 * @param x A 1D-array of size N.
 * @param N The size of a vector.
 */
template <typename T> inline
void Scal(T a, T *x, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { x[i] *= a; }
}
/**
 * @brief Compute y = a * x + y.
 * @tparam T The type of vectors.
 * @param a A scalar value.
 * @param x A 1D-array of size N.
 * @param y A 1D-array of size N.
 * @param N The size of vectors.
 */
template <typename T> inline
void Axpy(T a, T *x, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] += a * x[i]; }
}
/**
 * @brief Compute y = a * x + b * y.
 * @tparam T The type of vectors.
 * @param a A scalar value.
 * @param x A 1D-array of size N.
 * @param b A scalar value.
 * @param y A 1D-array of size N.
 * @param N The size of vectors.
 */
template <typename T> inline
void Axpby(T a, T *x, T b, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] = a * x[i] + b * y[i]; }
}
/**
 * @brief Compute z = a * x + y.
 * @tparam T The type of vectors.
 * @param a A scalar value.
 * @param x A 1D-array of size N.
 * @param y A 1D-array of size N.
 * @param z A 1D-array of size N.
 * @param N The size of vectors.
 */
template <typename T> inline
void Axpyz(T a, T *x, T *y, T *z, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { z[i] = a * x[i] + y[i]; }
}
/**
 * @brief Compute the dot product of x and y.
 * @tparam T The type of vectors.
 * @param x A 1D-array of size N.
 * @param y A 1D-array of size N.
 * @param N The size of vectors.
 * @return The resulting dot product.
 */
template <typename T> inline
T Dot(T *x, T *y, int N) {
    T res = 0;
    #pragma omp parallel for simd reduction(+: res)
    for(int i=0; i<N; i++) { res += x[i] * y[i]; }
    return res;
}
/**
 * @brief Compute the 2-norm of x.
 * @tparam T The type of the vector.
 * @param x A 1D-array of size N.
 * @param N The size of the vector.
 * @return 2-norm of x.
 */
template <typename T> inline
T Nrm2(T *x, int N) {
    T res = 0;
    #pragma omp parallel for simd reduction(+: res)
    for(int i=0; i<N; i++) { res += x[i] * x[i]; }
    return std::sqrt(res);
}
/**
 * @brief Compute the Hadamard product of x and y.
 * @tparam T The type of the vectors.
 * @param x A 1D-array of size N.
 * @param y A 1D-array of size N.
 * @param N The size of the vectors.
 */
template <typename T> inline
void HadProd(T *x, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] *= x[i]; }
}
/**
 * @brief Compute the element-wise division of x and y.
 * @tparam T The type of the vectors.
 * @param x A 1D-array of size N.
 * @param y A 1D-array of size N.
 * @param N The size of the vectors.
 */
template <typename T> inline
void HadDiv(T *x, T *y, int N) {
    #pragma omp parallel for simd
    for(int i=0; i<N; i++) { y[i] /= x[i]; }
}
/**
 * @brief Generate a Gives rotation matrix.
 * @tparam T The type of the vectors.
 * @param a A scalar value.
 * @param b A scalar value.
 * @param c The resulting value of cos.
 * @param s The resulting value of sin.
 */
template <typename T> inline
T Ggen(T a, T b, T *c, T *s)
{
    T r;
    r = std::sqrt(a*a + b*b);
    c[0] =  a / r;
    s[0] = -b / r;
    return r;
}
/**
 * @brief Compute the Gives rotation.
 * @tparam T The type of the vectors.
 * @param c The value of cos.
 * @param s THe value of sin.
 * @param a A rotated scalar value.
 * @param b A rotated scalar value.
 */
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
