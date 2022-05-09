/**
 * @file senk_gcr.hpp
 * @brief The GCR solvers are defined.
 * @author Kengo Suzuki
 * @date 5/9/2022
 */
#ifndef SENK_GCR_HPP
#define SENK_GCR_HPP

#include "senk_sparse.hpp"
#include "senk_blas1.hpp"

namespace senk {

namespace solver {
/**
 * @brief The Non-preconditioned GCR(m) solver.
 * @tparam T The type of a coefficient matrix and vectors.
 * @param val val array of the CSR storage format.
 * @param cind column index array of the CSR storage format.
 * @param rptr row pointer array of the CSR storage format.
 * @param b A right-hand side vector.
 * @param x An unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param outer The maximum number of iterations of outer loop.
 * @param m The number of the restart period.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T>
void Gcrm(
    T *val, int *cind, int *rptr,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *r      = new T[N];
    T *p      = new T[N*(m+1)];
    T *Ap     = new T[N*(m+1)];
    T *Ar     = new T[N];
    T *dot_Ap = new T[m+1];
    T alpha, beta, nrm_r;
    
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, r, N);
        blas1::Axpby<T>(1, b, -1, r, N);
        blas1::Copy<T>(r, &p[0], N);
        nrm_r = blas1::Nrm2<T>(r, N);
        printf("%d %e\n", i*m, nrm_r/nrm_b);
        if(nrm_r < nrm_b * epsilon) break;
        sparse::SpmvCsr<T>(val, cind, rptr, &p[0], &Ap[0], N);
        for(int j=0; j<m; j++) {
            dot_Ap[j] = blas1::Dot<T>(&Ap[j*N], &Ap[j*N], N);
            alpha = blas1::Dot<T>(&Ap[j*N], r, N) / dot_Ap[j];
            blas1::Axpy<T>(alpha, &p[j*N], x, N);
            blas1::Axpy<T>(-alpha, &Ap[j*N], r, N);
            sparse::SpmvCsr<T>(val, cind, rptr, r, Ar, N);
            beta = -blas1::Dot<T>(&Ap[0], Ar, N) / dot_Ap[0];
            blas1::Axpyz<T>(beta, &p[0], r, &p[(j+1)*N], N);
            blas1::Axpyz<T>(beta, &Ap[0], Ar, &Ap[(j+1)*N], N);
            for(int k=1; k<=j; k++) {
                beta = -blas1::Dot<T>(&Ap[k*N], Ar, N) / dot_Ap[k];
                blas1::Axpy<T>(beta, &p[k*N], &p[(j+1)*N], N);
                blas1::Axpy<T>(beta, &Ap[k*N], &Ap[(j+1)*N], N);
            }
        }
    }
    delete[] r;
    delete[] p;
    delete[] Ap;
    delete[] Ar;
    delete[] dot_Ap;
}
/**
 * @brief The ILU preconditioned GCR(m) solver.
 * @tparam T The type of a coefficient matrix and vectors.
 * @param val val array of the CSR storage format.
 * @param cind column index array of the CSR storage format.
 * @param rptr row pointer array of the CSR storage format.
 * @param lval Same as val, but for the matrix L.
 * @param lcind Same as cind, but for the matrix L.
 * @param lrptr Same as rptr, but for the matrix L.
 * @param uval Same as val, but for the matrix U.
 * @param ucind Same as cind, but for the matrix U.
 * @param urptr Same as rptr, but for the matrix U.
 * @param b A right-hand side vector.
 * @param x An unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param outer The maximum number of iterations of outer loop.
 * @param m The number of the restart period.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T>
void IluGcrm(
    T *val, int *cind, int *rptr,
    T *lval, int *lcind, int *lrptr,
    T *uval, int *ucind, int *urptr,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *r      = new T[N];
    T *Kr     = new T[N];
    T *p      = new T[N*(m+1)];
    T *Ap     = new T[N*(m+1)];
    T *AKr    = new T[N];
    T *dot_Ap = new T[m+1];
    T alpha, beta, nrm_r;
    
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, r, N);
        blas1::Axpby<T>(1, b, -1, r, N);
        sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, r, &p[0], N);
        sparse::SptrsvCsr_u<T>(uval, ucind, urptr, &p[0], &p[0], N);
        nrm_r = blas1::Nrm2<T>(r, N);
        printf("%d %e\n", i*m, nrm_r/nrm_b);
        if(nrm_r < nrm_b * epsilon) break;
        sparse::SpmvCsr<T>(val, cind, rptr, &p[0], &Ap[0], N);
        for(int j=0; j<m; j++) {
            dot_Ap[j] = blas1::Dot<T>(&Ap[j*N], &Ap[j*N], N);
            alpha = blas1::Dot<T>(&Ap[j*N], r, N) / dot_Ap[j];
            blas1::Axpy<T>(alpha, &p[j*N], x, N);
            blas1::Axpy<T>(-alpha, &Ap[j*N], r, N);
            sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, r, Kr, N);
            sparse::SptrsvCsr_u<T>(uval, ucind, urptr, Kr, Kr, N);
            sparse::SpmvCsr<T>(val, cind, rptr, Kr, AKr, N);
            beta = -blas1::Dot<T>(&Ap[0], AKr, N) / dot_Ap[0];
            blas1::Axpyz<T>(beta, &p[0], Kr, &p[(j+1)*N], N);
            blas1::Axpyz<T>(beta, &Ap[0], AKr, &Ap[(j+1)*N], N);
            for(int k=1; k<=j; k++) {
                beta = -blas1::Dot<T>(&Ap[k*N], AKr, N) / dot_Ap[k];
                blas1::Axpy<T>(beta, &p[k*N], &p[(j+1)*N], N);
                blas1::Axpy<T>(beta, &Ap[k*N], &Ap[(j+1)*N], N);
            }
        }
    }
    delete[] r;
    delete[] Kr;
    delete[] p;
    delete[] Ap;
    delete[] AKr;
    delete[] dot_Ap;
}
/**
 * @brief The ILUB preconditioned GCR(m) solver.
 * @tparam T The type of a coefficient matrix and vectors.
 * @param val val array of the CSR storage format.
 * @param cind column index array of the CSR storage format.
 * @param rptr row pointer array of the CSR storage format.
 * @param blval values of L in the BCSR format.
 * @param blcind colum positions of blocks of L in the BCSR format.
 * @param blrptr starting positions of row blocks of L in the BCSR format.
 * @param buval values of U in the BCSR format.
 * @param bucind colum positions of blocks of U in the BCSR format.
 * @param burptr starting positions of row blocks of U in the BCSR format.
 * @param b A right-hand side vector.
 * @param x An unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param outer The maximum number of iterations of outer loop.
 * @param m The number of the restart period.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T, int bnl, int bnw>
void IlubGcrm(
    T *val, int *cind, int *rptr,
    T *blval, int *blcind, int *blrptr,
    T *buval, int *bucind, int *burptr,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *r      = new T[N];
    T *Kr     = new T[N];
    T *p      = new T[N*(m+1)];
    T *Ap     = new T[N*(m+1)];
    T *AKr    = new T[N];
    T *dot_Ap = new T[m+1];
    T alpha, beta, nrm_r;
    
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, r, N);
        blas1::Axpby<T>(1, b, -1, r, N);
        sparse::SptrsvBcsr_l<T, bnl, bnw>(blval, blcind, blrptr, r, &p[0], N);
        sparse::SptrsvBcsr_u<T, bnl, bnw>(buval, bucind, burptr, &p[0], &p[0], N);
        nrm_r = blas1::Nrm2<T>(r, N);
        printf("%d %e\n", i*m, nrm_r/nrm_b);
        if(nrm_r < nrm_b * epsilon) break;
        sparse::SpmvCsr<T>(val, cind, rptr, &p[0], &Ap[0], N);
        for(int j=0; j<m; j++) {
            dot_Ap[j] = blas1::Dot<T>(&Ap[j*N], &Ap[j*N], N);
            alpha = blas1::Dot<T>(&Ap[j*N], r, N) / dot_Ap[j];
            blas1::Axpy<T>(alpha, &p[j*N], x, N);
            blas1::Axpy<T>(-alpha, &Ap[j*N], r, N);
            sparse::SptrsvBcsr_l<T, bnl, bnw>(blval, blcind, blrptr, r, Kr, N);
            sparse::SptrsvBcsr_u<T, bnl, bnw>(buval, bucind, burptr, Kr, Kr, N);
            sparse::SpmvCsr<T>(val, cind, rptr, Kr, AKr, N);
            beta = -blas1::Dot<T>(&Ap[0], AKr, N) / dot_Ap[0];
            blas1::Axpyz<T>(beta, &p[0], Kr, &p[(j+1)*N], N);
            blas1::Axpyz<T>(beta, &Ap[0], AKr, &Ap[(j+1)*N], N);
            for(int k=1; k<=j; k++) {
                beta = -blas1::Dot<T>(&Ap[k*N], AKr, N) / dot_Ap[k];
                blas1::Axpy<T>(beta, &p[k*N], &p[(j+1)*N], N);
                blas1::Axpy<T>(beta, &Ap[k*N], &Ap[(j+1)*N], N);
            }
        }
    }
    delete[] r;
    delete[] Kr;
    delete[] p;
    delete[] Ap;
    delete[] AKr;
    delete[] dot_Ap;
}

} // namespace solver

} // namespace senk

#endif
