/**
 * @file senk_bicgstab.hpp
 * @brief The BiCGStab method is defined.
 * @author Kengo Suzuki
 * @date 5/8/2021
 */
#ifndef SENK_BICGSTAB_HPP
#define SENK_BICGSTAB_HPP

#include "senk_sparse.hpp"
#include "senk_blas1.hpp"

namespace senk {

namespace solver {
/**
 * @brief Non-preconditioned BiCGStab solver
 * @tparam T The type of a coefficient matrix and vectors.
 * @param val val array of the CSR storage format.
 * @param cind col-index array of the CSR storage format.
 * @param rptr row-ptr array of the CSR storage format.
 * @param b The right-hand side vector.
 * @param x The unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param max_iter The maximum number of iterations.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T>
void Bicgstab(
    T *val, int *cind, int *rptr,
    T *b, T *x, T nrm_b,
    int max_iter, int N, T epsilon)
{
    int i;
    int flag = 0;
    T *r    = new T[N];
    T *rstr = new T[N];
    T *p    = new T[N];
    T *Ap   = new T[N];
    T *s    = new T[N];
    T *As   = new T[N];
    T *temp = new T[N];
    T alpha, beta, omega;
    T r_rstr, prev;
    T nrm_r = nrm_b;

    sparse::SpmvCsr<T>(val, cind, rptr, x, r, N);
    blas1::Axpby<T>(1, b, -1, r, N);
    blas1::Copy<T>(r, rstr, N);
    blas1::Copy<T>(r, p, N);
    r_rstr = blas1::Dot<T>(r, rstr, N);
    for(i=0; i<max_iter; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, p, Ap, N);
        alpha = r_rstr / blas1::Dot<T>(Ap, rstr, N);
        blas1::Axpyz<T>(-alpha, Ap, r, s, N);
        sparse::SpmvCsr<T>(val, cind, rptr, s, As, N);
        omega = blas1::Dot<T>(As, s, N) / blas1::Dot<T>(As, As, N);
        blas1::Axpy<T>(alpha, p, x, N);
        blas1::Axpy<T>(omega, s, x, N);
        blas1::Axpyz<T>(-omega, As, s, r, N);
        nrm_r = blas1::Nrm2<T>(r, N);
        printf("%d %e\n", i+1, nrm_r/nrm_b);
        if(nrm_r < epsilon * nrm_b) {
            printf("# iter %d\n", i+1);
            printf("# res %e\n", nrm_r/nrm_b);
            flag = 1;
            break;
        }
        prev = r_rstr;
        r_rstr = blas1::Dot<T>(r, rstr, N);
        beta = alpha / omega * r_rstr / prev;
        blas1::Axpyz<T>(-omega, Ap, p, temp, N);
        blas1::Axpyz<T>(beta, temp, r, p, N);
    }
    if(!flag) {
        printf("# iter %d (max)\n", i);
        printf("# res %e\n", nrm_r/nrm_b);
    }
    delete[] r;
    delete[] rstr;
    delete[] p;
    delete[] Ap;
    delete[] s;
    delete[] As;
    delete[] temp;
}
/**
 * @brief ILU preconditioned BiCGStab solver
 * @tparam T The type of a coefficient matrix and vectors.
 * @param val val array of the CSR storage format.
 * @param cind col-index array of the CSR storage format.
 * @param rptr row-ptr array of the CSR storage format.
 * @param lval Same as val, but for the matrix L.
 * @param lcind Same as cind, but for the matrix L.
 * @param lrptr Same as rptr, but for the matrix L.
 * @param uval Same as val, but for the matrix U.
 * @param ucind Same as cind, but for the matrix U.
 * @param urptr Same as rptr, but for the matrix U.
 * @param b The right-hand side vector.
 * @param x The unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param max_iter The maximum number of iterations.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T>
void IluBicgstab(
    T *val, int *cind, int *rptr,
    T *lval, int *lcind, int *lrptr,
    T *uval, int *ucind, int *urptr,
    T *b, T *x, T nrm_b,
    int max_iter, int N, T epsilon)
{
    int i;
    int flag = 0;
    T *r    = new T[N];
    T *rstr = new T[N];
    T *p    = new T[N];
    T *Kp   = new T[N];
    T *AKp  = new T[N];
    T *s    = new T[N];
    T *Ks   = new T[N];
    T *AKs  = new T[N];
    T *temp = new T[N];

    T alpha, beta, omega;
    T r_rstr, prev;
    T nrm_r = nrm_b;

    sparse::SpmvCsr<T>(val, cind, rptr, x, r, N);
    blas1::Axpby<T>(1, b, -1, r, N);
    blas1::Copy<T>(r, rstr, N);
    blas1::Copy<T>(r, p, N);
    r_rstr = blas1::Dot<T>(r, rstr, N);
    for(i=0; i<max_iter; i++) {
        sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, p, Kp, N);
        sparse::SptrsvCsr_u<T>(uval, ucind, urptr, Kp, Kp, N);
        sparse::SpmvCsr<T>(val, cind, rptr, Kp, AKp, N);
        alpha = r_rstr / blas1::Dot<T>(AKp, rstr, N);
        blas1::Axpyz(-alpha, AKp, r, s, N);
        sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, s, Ks, N);
        sparse::SptrsvCsr_u<T>(uval, ucind, urptr, Ks, Ks, N);
        sparse::SpmvCsr<T>(val, cind, rptr, Ks, AKs, N);
        omega = blas1::Dot<T>(AKs, s, N) / blas1::Dot<T>(AKs, AKs, N);
        blas1::Axpy<T>(alpha, Kp, x, N);
        blas1::Axpy<T>(omega, Ks, x, N);
        blas1::Axpyz<T>(-omega, AKs, s, r, N);
        nrm_r = blas1::Nrm2<T>(r, N);
        printf("%d %e\n", i+1, nrm_r/nrm_b);
        if(nrm_r < epsilon * nrm_b) {
            printf("# iter %d\n", i+1);
            printf("# res %e\n", nrm_r/nrm_b);
            flag = 1;
            break;
        }
        prev = r_rstr;
        r_rstr = blas1::Dot<T>(r, rstr, N);
        beta = alpha / omega * r_rstr / prev;
        blas1::Axpyz<T>(-omega, AKp, p, temp, N);
        blas1::Axpyz<T>(beta, temp, r, p, N);
    }
    if(!flag) {
        printf("# iter %d (max)\n", i);
        printf("# res %e\n", nrm_r/nrm_b);
    }
    delete[] r;
    delete[] rstr;
    delete[] p;
    delete[] Kp;
    delete[] AKp;
    delete[] s;
    delete[] Ks;
    delete[] AKs;
    delete[] temp;
}
/**
 * @brief ILUB preconditioned BiCGStab solver
 * @tparam T The type of a coefficient matrix and vectors.
 * @tparam bnl The number of rows of the block.
 * @tparam bnw The number of columns of the block.
 * @param val val array of the CSR storage format.
 * @param cind col-index array of the CSR storage format.
 * @param rptr row-ptr array of the CSR storage format.
 * @param blval values of L in the BCSR format.
 * @param blcind colum positions of blocks of L in the BCSR format.
 * @param blrptr starting positions of row blocks of L in the BCSR format.
 * @param buval values of U in the BCSR format.
 * @param bucind colum positions of blocks of U in the BCSR format.
 * @param burptr starting positions of row blocks of U in the BCSR format.
 * @param b The right-hand side vector.
 * @param x The unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param max_iter The maximum number of iterations.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T, int bnl, int bnw>
void IlubBicgstab(
    T *val, int *cind, int *rptr,
    T *blval, int *blcind, int *blrptr,
    T *buval, int *bucind, int *burptr,
    T *b, T *x, T nrm_b,
    int max_iter, int N, T epsilon)
{
    int i;
    int flag = 0;
    T *r    = new T[N];
    T *rstr = new T[N];
    T *p    = new T[N];
    T *Kp   = new T[N];
    T *AKp  = new T[N];
    T *s    = new T[N];
    T *Ks   = new T[N];
    T *AKs  = new T[N];
    T *temp = new T[N];

    T alpha, beta, omega;
    T r_rstr, prev;
    T nrm_r = nrm_b;

    sparse::SpmvCsr<T>(val, cind, rptr, x, r, N);
    blas1::Axpby<T>(1, b, -1, r, N);
    blas1::Copy<T>(r, rstr, N);
    blas1::Copy<T>(r, p, N);
    r_rstr = blas1::Dot<T>(r, rstr, N);
    for(i=0; i<max_iter; i++) {
        sparse::SptrsvBcsr_l<T, bnl, bnw>(blval, blcind, blrptr, p, Kp, N);
        sparse::SptrsvBcsr_u<T, bnl, bnw>(buval, bucind, burptr, Kp, Kp, N);
        sparse::SpmvCsr<T>(val, cind, rptr, Kp, AKp, N);
        alpha = r_rstr / blas1::Dot(AKp, rstr, N);
        blas1::Axpyz(-alpha, AKp, r, s, N);
        sparse::SptrsvBcsr_l<T, bnl, bnw>(blval, blcind, blrptr, s, Ks, N);
        sparse::SptrsvBcsr_u<T, bnl, bnw>(buval, bucind, burptr, Ks, Ks, N);
        sparse::SpmvCsr<T>(val, cind, rptr, Ks, AKs, N);
        omega = blas1::Dot<T>(AKs, s, N) / blas1::Dot<T>(AKs, AKs, N);
        blas1::Axpy<T>(alpha, Kp, x, N);
        blas1::Axpy<T>(omega, Ks, x, N);
        blas1::Axpyz<T>(-omega, AKs, s, r, N);
        nrm_r = blas1::Nrm2<T>(r, N);
        printf("%d %e\n", i+1, nrm_r/nrm_b);
        if(nrm_r < epsilon * nrm_b) {
            printf("# iter %d\n", i+1);
            printf("# res %e\n", nrm_r/nrm_b);
            flag = 1;
            break;
        }
        prev = r_rstr;
        r_rstr = blas1::Dot<T>(r, rstr, N);
        beta = alpha / omega * r_rstr / prev;
        blas1::Axpyz<T>(-omega, AKp, p, temp, N);
        blas1::Axpyz<T>(beta, temp, r, p, N);
    }
    if(!flag) {
        printf("# iter %d (max)\n", i);
        printf("# res %e\n", nrm_r/nrm_b);
    }
    delete[] r;
    delete[] rstr;
    delete[] p;
    delete[] Kp;
    delete[] AKp;
    delete[] s;
    delete[] Ks;
    delete[] AKs;
    delete[] temp;
}

//
//void Bicgstab_IR(
//    double *val, int *fval, int *cind, int *rptr,
//    double *b, double *x, double nrm_b,
//    int max_iter, int N, double epsilon);
//void Bicgstab_IR(
//    double *val, short *fval, int *cind, int *rptr,
//    double *b, double *x, double nrm_b,
//    int max_iter, int N, double epsilon);
//
//void IluBicgstab(
//    float *val, int *cind, int *rptr,
//    float *lval, int *lcind, int *lrptr,
//    float *uval, int *ucind, int *urptr,
//    float *b, float *x, double nrm_b,
//    int max_iter, int N, double epsilon);
//void IluBicgstab_IR(
//    double *val, float *fval, int *cind, int *rptr,
//    float *lval, int *lcind, int *lrptr,
//    float *uval, int *ucind, int *urptr,
//    double *b, double *x, double nrm_b,
//    int max_iter, int N, double epsilon);
//void IluBicgstab_IR(
//    double *val, int *fval, int *cind, int *rptr,
//    int *lval, int *lcind, int *lrptr,
//    int *uval, int *ucind, int *urptr,
//    double *b, double *x, double nrm_b,
//    int max_iter, int N, double epsilon);

} // namespace solver

} // namespace senk

#endif
