#ifndef SENK_GMRES_HPP
#define SENK_GMRES_HPP

#include "senk_sparse.hpp"
#include "senk_blas1.hpp"
#include "senk_blas2.hpp"

namespace senk {

namespace solver {

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
        blas1::Axpby<T>(1, b, -1, &V[0], N);
        e[0] = blas1::Nrm2<T>(&V[0], N);
        blas1::Scal<T>(1/e[0], &V[0], N);
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
            printf("e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("# iter %d\n", i*m+j+1);
                printf("# res %e\n", std::abs(e[j+1])/nrm_b);
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

    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
}

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
            printf("e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("# iter %d\n", i*m+j+1);
                printf("# res %e\n", std::abs(e[j+1])/nrm_b);
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

    delete[] c;
    delete[] s;
    delete[] e;
    delete[] H;
    delete[] V;
    delete[] y;
    delete[] t;
}

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
            printf("e[%d] = %e\n", j+1, std::abs(e[j+1]/nrm_b));
            if(std::abs(e[j+1]) <= nrm_b*epsilon) {
                printf("# iter %d\n", i*m+j+1);
                printf("# res %e\n", std::abs(e[j+1])/nrm_b);
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
