/**
 * @file senk_gmres.hpp
 * @brief The GMRES solvers are defined.
 * @author Kengo Suzuki
 * @date 5/9/2022
 */
#ifndef SENK_GMRES_HPP
#define SENK_GMRES_HPP

#include "senk_sparse.hpp"
#include "senk_blas1.hpp"
#include "senk_blas2.hpp"

namespace senk {
/**
 * @brief Contains solvers.
 */
namespace solver {
/**
 * @brief The Non-preconditioned GMRES(m) solver.
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
void Gmresm(
    T *val, int *cind, int *rptr,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *c = new T[m];
    T *s = new T[m];
    T *e = new T[m+1];
    T *H = new T[(m+1)*m];
    T *V = new T[N*(m+1)];
    T *y = new T[m];

    int flag = 0;
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, &V[0], N);
        senk::blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        senk::blas1::Scal<T>(1/e[0], &V[0], N);
        int j;
        for(j=0; j<m; j++) {
            sparse::SpmvCsr<T>(val, cind, rptr, &V[j*N], &V[(j+1)*N], N);
            for(int k=0; k<=j; k++) {
                H[j*(m+1)+k] = blas1::Dot<T>(&V[k*N], &V[(j+1)*N], N);
                blas1::Axpy<T>(-H[j*(m+1)+k], &V[k*N], &V[(j+1)*N], N);
            }
            H[j*(m+1)+j+1] = blas1::Nrm2<T>(&V[(j+1)*N], N);
            blas1::Scal<T>(1/H[j*(m+1)+j+1], &V[(j+1)*N], N);
            for(int k=0; k<j; k++) {
                blas1::Grot<T>(c[k], s[k], &H[j*(m+1)+k], &H[j*(m+1)+k+1]);
            }
            H[j*(m+1)+j] = blas1::Ggen<T>(H[j*(m+1)+j], H[j*(m+1)+j+1], &c[j], &s[j]);
            H[j*(m+1)+j+1] = 0;
            e[j+1] = s[j] * e[j];
            e[j] = c[j] * e[j];
#if PRINT_RES
            printf("# e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
#endif
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("%s iter %d\n", ITER_SYMBOL, i*m+j+1);
                printf("%s res %e\n", RES_SYMBOL, std::abs(e[j+1])/nrm_b);
                j++;
                flag = 1;
                break;
            }
        }
        blas2::Trsv<T>(H, e, y, m+1, j);
        for(int k=0; k<j; k++) {
            blas1::Axpy<T>(y[k], &V[k*N], x, N);
        }
        if(flag == 1) break;
    }
    if(!flag) {
        printf("# iter %d\n", outer*m);
        printf("# res : Check by using senk::test\n");
    }
    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
}

/**
 * @brief The Non-preconditioned GMRES(m) solver.
 * @tparam T The type of a coefficient matrix and vectors.
 * @tparam bnl The number of rows of the block.
 * @tparam bnw The number of columns of the block.
 * @param bval val array of the BCSR storage format.
 * @param bcind column index array of the BCSR storage format.
 * @param brptr row pointer array of the BCSR storage format.
 * @param b A right-hand side vector.
 * @param x An unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param outer The maximum number of iterations of outer loop.
 * @param m The number of the restart period.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T, int bnl, int bnw>
void Gmresm(
    T *bval, int *bcind, int *brptr,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *c = new T[m];
    T *s = new T[m];
    T *e = new T[m+1];
    T *H = new T[(m+1)*m];
    T *V = new T[N*(m+1)];
    T *y = new T[m];

    int flag = 0;
    for(int i=0; i<outer; i++) {
        sparse::SpmvBcsr<T, bnl, bnw>(bval, bcind, brptr, x, &V[0], N);
        senk::blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        senk::blas1::Scal<T>(1/e[0], &V[0], N);
        int j;
        for(j=0; j<m; j++) {
            sparse::SpmvBcsr<T, bnl, bnw>(bval, bcind, brptr, &V[j*N], &V[(j+1)*N], N);
            for(int k=0; k<=j; k++) {
                H[j*(m+1)+k] = blas1::Dot<T>(&V[k*N], &V[(j+1)*N], N);
                blas1::Axpy<T>(-H[j*(m+1)+k], &V[k*N], &V[(j+1)*N], N);
            }
            H[j*(m+1)+j+1] = blas1::Nrm2<T>(&V[(j+1)*N], N);
            blas1::Scal<T>(1/H[j*(m+1)+j+1], &V[(j+1)*N], N);
            for(int k=0; k<j; k++) {
                blas1::Grot<T>(c[k], s[k], &H[j*(m+1)+k], &H[j*(m+1)+k+1]);
            }
            H[j*(m+1)+j] = blas1::Ggen<T>(H[j*(m+1)+j], H[j*(m+1)+j+1], &c[j], &s[j]);
            H[j*(m+1)+j+1] = 0;
            e[j+1] = s[j] * e[j];
            e[j] = c[j] * e[j];
#if PRINT_RES
            printf("# e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
#endif
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("%s iter %d\n", ITER_SYMBOL, i*m+j+1);
                printf("%s res %e\n", RES_SYMBOL, std::abs(e[j+1])/nrm_b);
                j++;
                flag = 1;
                break;
            }
        }
        blas2::Trsv<T>(H, e, y, m+1, j);
        for(int k=0; k<j; k++) {
            blas1::Axpy<T>(y[k], &V[k*N], x, N);
        }
        if(flag == 1) break;
    }
    if(!flag) {
        printf("# iter %d\n", outer*m);
        printf("# res : Check by using senk::test\n");
    }
    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
}

/**
 * @brief The ILU preconditioned GMRES(m) solver.
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
void IluGmresm(
    T *val, int *cind, int *rptr,
    T *lval, int *lcind, int *lrptr,
    T *uval, int *ucind, int *urptr,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *c = new T[m];
    T *s = new T[m];
    T *e = new T[m+1];
    T *H = new T[(m+1)*m];
    T *V = new T[N*(m+1)];
    T *y = new T[m];
    T *t = new T[N];

    int flag = 0;
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, &V[0], N);
        blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        blas1::Scal<T>(1/e[0], &V[0], N);
        int j;
        for(j=0; j<m; j++) {
            sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, &V[j*N], t, N);
            sparse::SptrsvCsr_u<T>(uval, ucind, urptr, t, t, N);
            sparse::SpmvCsr<T>(val, cind, rptr, t, &V[(j+1)*N], N);
            for(int k=0; k<=j; k++) {
                H[j*(m+1)+k] = blas1::Dot<T>(&V[k*N], &V[(j+1)*N], N);
                blas1::Axpy<T>(-H[j*(m+1)+k], &V[k*N], &V[(j+1)*N], N);
            }
            H[j*(m+1)+j+1] = blas1::Nrm2<T>(&V[(j+1)*N], N);
            blas1::Scal<T>(1/H[j*(m+1)+j+1], &V[(j+1)*N], N);
            for(int k=0; k<j; k++) {
                blas1::Grot<T>(c[k], s[k], &H[j*(m+1)+k], &H[j*(m+1)+k+1]);
            }
            H[j*(m+1)+j] = blas1::Ggen<T>(H[j*(m+1)+j], H[j*(m+1)+j+1], &c[j], &s[j]);
            H[j*(m+1)+j+1] = 0;
            e[j+1] = s[j] * e[j];
            e[j] = c[j] * e[j];
#if PRINT_RES
            printf("# e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
#endif
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("%s iter %d\n", ITER_SYMBOL, i*m+j+1);
                printf("%s res %e\n", RES_SYMBOL, std::abs(e[j+1])/nrm_b);
                j++;
                flag = 1;
                break;
            }
        }
        blas2::Trsv<T>(H, e, y, m+1, j);
        blas1::Scal<T>(y[0], &V[0], N);
        for(int k=1; k<j; k++) {
            blas1::Axpy<T>(y[k], &V[k*N], &V[0], N);
        }
        sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, &V[0], t, N);
        sparse::SptrsvCsr_u<T>(uval, ucind, urptr, t, t, N);
        blas1::Axpy<T>(1, t, x, N);

        if(flag == 1) break;
    }
    if(!flag) {
        printf("# iter %d\n", outer*m);
        printf("# res : Check by using senk::test\n");
    }
    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
    delete[] t;
}
/**
 * @brief The ILU preconditioned GMRES(m) solver parallelized by AMC ordering.
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
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 * @param b A right-hand side vector.
 * @param x An unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param outer The maximum number of iterations of outer loop.
 * @param m The number of the restart period.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T>
void AmcIluGmresm(
    T *val, int *cind, int *rptr,
    T *lval, int *lcind, int *lrptr,
    T *uval, int *ucind, int *urptr,
    int *cptr, int cnum,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *c = new T[m];
    T *s = new T[m];
    T *e = new T[m+1];
    T *H = new T[(m+1)*m];
    T *V = new T[N*(m+1)];
    T *y = new T[m];
    T *t = new T[N];

    int flag = 0;
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, &V[0], N);
        blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        blas1::Scal<T>(1/e[0], &V[0], N);
        int j;
        for(j=0; j<m; j++) {
            sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, &V[j*N], t, N, cptr, cnum);
            sparse::SptrsvCsr_u<T>(uval, ucind, urptr, t, t, N, cptr, cnum);
            sparse::SpmvCsr<T>(val, cind, rptr, t, &V[(j+1)*N], N);
            for(int k=0; k<=j; k++) {
                H[j*(m+1)+k] = blas1::Dot<T>(&V[k*N], &V[(j+1)*N], N);
                blas1::Axpy<T>(-H[j*(m+1)+k], &V[k*N], &V[(j+1)*N], N);
            }
            H[j*(m+1)+j+1] = blas1::Nrm2<T>(&V[(j+1)*N], N);
            blas1::Scal<T>(1/H[j*(m+1)+j+1], &V[(j+1)*N], N);
            for(int k=0; k<j; k++) {
                blas1::Grot<T>(c[k], s[k], &H[j*(m+1)+k], &H[j*(m+1)+k+1]);
            }
            H[j*(m+1)+j] = blas1::Ggen<T>(H[j*(m+1)+j], H[j*(m+1)+j+1], &c[j], &s[j]);
            H[j*(m+1)+j+1] = 0;
            e[j+1] = s[j] * e[j];
            e[j] = c[j] * e[j];
#if PRINT_RES
            printf("# e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
#endif
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("%s iter %d\n", ITER_SYMBOL, i*m+j+1);
                printf("%s res %e\n", RES_SYMBOL, std::abs(e[j+1])/nrm_b);
                j++;
                flag = 1;
                break;
            }
        }
        blas2::Trsv<T>(H, e, y, m+1, j);
        blas1::Scal<T>(y[0], &V[0], N);
        for(int k=1; k<j; k++) {
            blas1::Axpy<T>(y[k], &V[k*N], &V[0], N);
        }
        sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, &V[0], t, N, cptr, cnum);
        sparse::SptrsvCsr_u<T>(uval, ucind, urptr, t, t, N, cptr, cnum);
        blas1::Axpy<T>(1, t, x, N);

        if(flag == 1) break;
    }
    if(!flag) {
        printf("# iter %d\n", outer*m);
        printf("# res : Check by using senk::test\n");
    }
    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
    delete[] t;
}
/**
 * @brief The ILU preconditioned GMRES(m) solver parallelized by ABMC ordering.
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
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 * @param bsize The number of rows/columns of the blocks used in ABMC.
 * @param b A right-hand side vector.
 * @param x An unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param outer The maximum number of iterations of outer loop.
 * @param m The number of the restart period.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T>
void AbmcIluGmresm(
    T *val, int *cind, int *rptr,
    T *lval, int *lcind, int *lrptr,
    T *uval, int *ucind, int *urptr,
    int *cptr, int cnum, int bsize,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *c = new T[m];
    T *s = new T[m];
    T *e = new T[m+1];
    T *H = new T[(m+1)*m];
    T *V = new T[N*(m+1)];
    T *y = new T[m];
    T *t = new T[N];

    int flag = 0;
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, &V[0], N);
        blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        blas1::Scal<T>(1/e[0], &V[0], N);
        int j;
        for(j=0; j<m; j++) {
            sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, &V[j*N], t, N, cptr, cnum, bsize);
            sparse::SptrsvCsr_u<T>(uval, ucind, urptr, t, t, N, cptr, cnum, bsize);
            sparse::SpmvCsr<T>(val, cind, rptr, t, &V[(j+1)*N], N);
            for(int k=0; k<=j; k++) {
                H[j*(m+1)+k] = blas1::Dot<T>(&V[k*N], &V[(j+1)*N], N);
                blas1::Axpy<T>(-H[j*(m+1)+k], &V[k*N], &V[(j+1)*N], N);
            }
            H[j*(m+1)+j+1] = blas1::Nrm2<T>(&V[(j+1)*N], N);
            blas1::Scal<T>(1/H[j*(m+1)+j+1], &V[(j+1)*N], N);
            for(int k=0; k<j; k++) {
                blas1::Grot<T>(c[k], s[k], &H[j*(m+1)+k], &H[j*(m+1)+k+1]);
            }
            H[j*(m+1)+j] = blas1::Ggen<T>(H[j*(m+1)+j], H[j*(m+1)+j+1], &c[j], &s[j]);
            H[j*(m+1)+j+1] = 0;
            e[j+1] = s[j] * e[j];
            e[j] = c[j] * e[j];
#if PRINT_RES
            printf("# e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
#endif
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("%s iter %d\n", ITER_SYMBOL, i*m+j+1);
                printf("%s res %e\n", RES_SYMBOL, std::abs(e[j+1])/nrm_b);
                j++;
                flag = 1;
                break;
            }
        }
        blas2::Trsv<T>(H, e, y, m+1, j);
        blas1::Scal<T>(y[0], &V[0], N);
        for(int k=1; k<j; k++) {
            blas1::Axpy<T>(y[k], &V[k*N], &V[0], N);
        }
        sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, &V[0], t, N, cptr, cnum, bsize);
        sparse::SptrsvCsr_u<T>(uval, ucind, urptr, t, t, N, cptr, cnum, bsize);
        blas1::Axpy<T>(1, t, x, N);

        if(flag == 1) break;
    }
    if(!flag) {
        printf("# iter %d\n", outer*m);
        printf("# res : Check by using senk::test\n");
    }
    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
    delete[] t;
}
/**
 * @brief The ILU preconditioned GMRES(m) solver parallelized by the block Jacobi method.
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
 * @param bnum The number of blocks.
 * @param b A right-hand side vector.
 * @param x An unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param outer The maximum number of iterations of outer loop.
 * @param m The number of the restart period.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T>
void BjIluGmresm(
    T *val, int *cind, int *rptr,
    T *lval, int *lcind, int *lrptr,
    T *uval, int *ucind, int *urptr,
    int bnum,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *c = new T[m];
    T *s = new T[m];
    T *e = new T[m+1];
    T *H = new T[(m+1)*m];
    T *V = new T[N*(m+1)];
    T *y = new T[m];
    T *t = new T[N];

    int flag = 0;
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, &V[0], N);
        blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        blas1::Scal<T>(1/e[0], &V[0], N);
        int j;
        for(j=0; j<m; j++) {
            sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, &V[j*N], t, N, bnum);
            sparse::SptrsvCsr_u<T>(uval, ucind, urptr, t, t, N, bnum);
            sparse::SpmvCsr<T>(val, cind, rptr, t, &V[(j+1)*N], N);
            for(int k=0; k<=j; k++) {
                H[j*(m+1)+k] = blas1::Dot<T>(&V[k*N], &V[(j+1)*N], N);
                blas1::Axpy<T>(-H[j*(m+1)+k], &V[k*N], &V[(j+1)*N], N);
            }
            H[j*(m+1)+j+1] = blas1::Nrm2<T>(&V[(j+1)*N], N);
            blas1::Scal<T>(1/H[j*(m+1)+j+1], &V[(j+1)*N], N);
            for(int k=0; k<j; k++) {
                blas1::Grot<T>(c[k], s[k], &H[j*(m+1)+k], &H[j*(m+1)+k+1]);
            }
            H[j*(m+1)+j] = blas1::Ggen<T>(H[j*(m+1)+j], H[j*(m+1)+j+1], &c[j], &s[j]);
            H[j*(m+1)+j+1] = 0;
            e[j+1] = s[j] * e[j];
            e[j] = c[j] * e[j];
#if PRINT_RES
            printf("# e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
#endif
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("%s iter %d\n", ITER_SYMBOL, i*m+j+1);
                printf("%s res %e\n", RES_SYMBOL, std::abs(e[j+1])/nrm_b);
                j++;
                flag = 1;
                break;
            }
        }
        blas2::Trsv<T>(H, e, y, m+1, j);
        blas1::Scal<T>(y[0], &V[0], N);
        for(int k=1; k<j; k++) {
            blas1::Axpy<T>(y[k], &V[k*N], &V[0], N);
        }
        sparse::SptrsvCsr_l<T>(lval, lcind, lrptr, &V[0], t, N, bnum);
        sparse::SptrsvCsr_u<T>(uval, ucind, urptr, t, t, N, bnum);
        blas1::Axpy<T>(1, t, x, N);

        if(flag == 1) break;
    }
    if(!flag) {
        printf("# iter %d\n", outer*m);
        printf("# res : Check by using senk::test\n");
    }
    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
    delete[] t;
}
/**
 * @brief The ILUB preconditioned GMRES(m) solver.
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
void IlubGmresm(
    T *val, int *cind, int *rptr,
    T *blval, int *blcind, int *blrptr,
    T *buval, int *bucind, int *burptr,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *c = new T[m];
    T *s = new T[m];
    T *e = new T[m+1];
    T *H = new T[(m+1)*m];
    T *V = new T[N*(m+1)];
    T *y = new T[m];
    T *t = new T[N];

    int flag = 0;
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, &V[0], N);
        blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        blas1::Scal(1/e[0], &V[0], N);
        int j;
        for(j=0; j<m; j++) {
            sparse::SptrsvBcsr_l<T, bnl, bnw>(blval, blcind, blrptr, &V[j*N], t, N);
            sparse::SptrsvBcsr_u<T, bnl, bnw>(buval, bucind, burptr, t, t, N);
            sparse::SpmvCsr<T>(val, cind, rptr, t, &V[(j+1)*N], N);
            for(int k=0; k<=j; k++) {
                H[j*(m+1)+k] = blas1::Dot<T>(&V[k*N], &V[(j+1)*N], N);
                blas1::Axpy<T>(-H[j*(m+1)+k], &V[k*N], &V[(j+1)*N], N);
            }
            H[j*(m+1)+j+1] = blas1::Nrm2<T>(&V[(j+1)*N], N);
            blas1::Scal<T>(1/H[j*(m+1)+j+1], &V[(j+1)*N], N);
            for(int k=0; k<j; k++) {
                blas1::Grot<T>(c[k], s[k], &H[j*(m+1)+k], &H[j*(m+1)+k+1]);
            }
            H[j*(m+1)+j] = blas1::Ggen<T>(H[j*(m+1)+j], H[j*(m+1)+j+1], &c[j], &s[j]);
            H[j*(m+1)+j+1] = 0;
            e[j+1] = s[j] * e[j];
            e[j] = c[j] * e[j];
#if PRINT_RES
            printf("# e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
#endif
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("%s iter %d\n", ITER_SYMBOL, i*m+j+1);
                printf("%s res %e\n", RES_SYMBOL, std::abs(e[j+1])/nrm_b);
                j++;
                flag = 1;
                break;
            }
        }
        blas2::Trsv<T>(H, e, y, m+1, j);
        blas1::Scal<T>(y[0], &V[0], N);
        for(int k=1; k<j; k++) {
            blas1::Axpy<T>(y[k], &V[k*N], &V[0], N);
        }
        sparse::SptrsvBcsr_l<T, bnl, bnw>(blval, blcind, blrptr, &V[0], t, N);
        sparse::SptrsvBcsr_u<T, bnl, bnw>(buval, bucind, burptr, t, t, N);
        blas1::Axpy<T>(1, t, x, N);

        if(flag == 1) break;
    }
    if(!flag) {
        printf("# iter %d\n", outer*m);
        printf("# res : Check by using senk::test\n");
    }
    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
    delete[] t;
}
/**
 * @brief The ILUB preconditioned GMRES(m) solver parallelized by ABMC ordering.
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
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 * @param bsize The number of rows/columns of the blocks used in ABMC.
 * @param b A right-hand side vector.
 * @param x An unknown vector.
 * @param nrm_b The 2-norm of b.
 * @param outer The maximum number of iterations of outer loop.
 * @param m The number of the restart period.
 * @param N The size of the matrix and the vectors.
 * @param epsilon The convergence criterion.
 */
template <typename T, int bnl, int bnw>
void AbmcIlubGmresm(
    T *val, int *cind, int *rptr,
    T *blval, int *blcind, int *blrptr,
    T *buval, int *bucind, int *burptr,
    int *cptr, int cnum, int bsize,
    T *b, T *x, T nrm_b,
    int outer, int m, int N, T epsilon)
{
    T *c = new T[m];
    T *s = new T[m];
    T *e = new T[m+1];
    T *H = new T[(m+1)*m];
    T *V = new T[N*(m+1)];
    T *y = new T[m];
    T *t = new T[N];

    int flag = 0;
    for(int i=0; i<outer; i++) {
        sparse::SpmvCsr<T>(val, cind, rptr, x, &V[0], N);
        blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        blas1::Scal(1/e[0], &V[0], N);
        int j;
        for(j=0; j<m; j++) {
            sparse::SptrsvBcsr_l<T, bnl, bnw>(blval, blcind, blrptr, &V[j*N], t, N, cptr, cnum, bsize);
            sparse::SptrsvBcsr_u<T, bnl, bnw>(buval, bucind, burptr, t, t, N, cptr, cnum, bsize);
            sparse::SpmvCsr<T>(val, cind, rptr, t, &V[(j+1)*N], N);
            for(int k=0; k<=j; k++) {
                H[j*(m+1)+k] = blas1::Dot<T>(&V[k*N], &V[(j+1)*N], N);
                blas1::Axpy<T>(-H[j*(m+1)+k], &V[k*N], &V[(j+1)*N], N);
            }
            H[j*(m+1)+j+1] = blas1::Nrm2<T>(&V[(j+1)*N], N);
            blas1::Scal<T>(1/H[j*(m+1)+j+1], &V[(j+1)*N], N);
            for(int k=0; k<j; k++) {
                blas1::Grot<T>(c[k], s[k], &H[j*(m+1)+k], &H[j*(m+1)+k+1]);
            }
            H[j*(m+1)+j] = blas1::Ggen<T>(H[j*(m+1)+j], H[j*(m+1)+j+1], &c[j], &s[j]);
            H[j*(m+1)+j+1] = 0;
            e[j+1] = s[j] * e[j];
            e[j] = c[j] * e[j];
#if PRINT_RES
            printf("# e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
#endif
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("%s iter %d\n", ITER_SYMBOL, i*m+j+1);
                printf("%s res %e\n", RES_SYMBOL, std::abs(e[j+1])/nrm_b);
                j++;
                flag = 1;
                break;
            }
        }
        blas2::Trsv<T>(H, e, y, m+1, j);
        blas1::Scal<T>(y[0], &V[0], N);
        for(int k=1; k<j; k++) {
            blas1::Axpy<T>(y[k], &V[k*N], &V[0], N);
        }
        sparse::SptrsvBcsr_l<T, bnl, bnw>(blval, blcind, blrptr, &V[0], t, N, cptr, cnum, bsize);
        sparse::SptrsvBcsr_u<T, bnl, bnw>(buval, bucind, burptr, t, t, N, cptr, cnum, bsize);
        blas1::Axpy<T>(1, t, x, N);

        if(flag == 1) break;
    }
    if(!flag) {
        printf("# iter %d\n", outer*m);
        printf("# res : Check by using senk::test\n");
    }
    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
    delete[] t;
}

} // namespace solver

} // namespace senk

#endif
