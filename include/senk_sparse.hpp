#ifndef SENK_SPARSE_HPP
#define SENK_SPARSE_HPP

#include <vector>
#include <iostream>
//#include <omp.h>

namespace senk {

namespace sparse {

enum Triangle {
    Upper,
    Lower,
    UnitLower
};

template <typename T> inline
void SpmvCsr(T *val, int *cind, int *rptr, T *x, T *y, int N) {
    #pragma omp parallel for
    for(int i=0; i<N; i++) {
        T temp = 0;
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            temp += val[j] * x[cind[j]];
        }
        y[i] = temp;
    }
}

template <typename T> inline
void SpmvCsr(T *val, int *cind, int *rptr, T *diag, T *x, T *y, int N) {
    #pragma omp parallel for
    for(int i=0; i<N; i++) {
        T temp = x[i] * diag[i];
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            temp += val[j] * x[cind[j]];
        }
        y[i] = temp;
    }
}

template <typename T> inline
void SpmvSell(T *val, int *cind, int *wid, int len, T *x, T *y, int N)
{
    int block = (N+len-1)/len;
    #pragma omp parallel for
    for(int i=0; i<block; i++) {
        int start = wid[i] * len;
        int temp = (i==len-1 && N%len!=0) ? N % len : len;
        for(int k=0; k<temp; k++) {
            y[i*len+k] = val[start+k] * x[cind[start+k]];
        }
        for(int j=1; j<wid[i+1]-wid[i]; j++) {
            int off = start+j*len;
            for(int k=0; k<temp; k++) {
                y[i*len+k] += val[off+k] * x[cind[off+k]];
            }
        }
    }
}

template <typename T> inline
void SptrsvCsr_l(T *val, int *cind, int *rptr, T *x, T *y, int N)
{
    // L is assumed to be unit lower triangular.
    for(int i=0; i<N; i++) {
        T temp = x[i];
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            temp -= val[j] * y[cind[j]];
        }
        y[i] = temp;
    }
}

template <typename T> inline
void SptrsvCsr_u(T *val, int *cind, int *rptr, T *x, T *y, int N)
{
    // L is assumed to be general upper triangular.
    // Diagonal has been inverted.
    for(int i=N-1; i>=0; i--) {
        T temp = x[i];
        int j;
        for(j=rptr[i+1]-1; j>=rptr[i]+1; j--) {
            temp -= val[j] * y[cind[j]];    
        }
        y[i] = temp * val[j];
    }
}

template <typename T, int bnl, int bnw> inline
void SptrsvBcsr_l(
    T *bval, int *bcind, int *brptr, 
    T *x, T *y, int N)
{
    // L is assumed to be unit lower triangular.
    int b_size = bnl * bnw;
    for(int i=0; i<N; i+=bnl) {
        int bidx = i / bnl;
        #pragma omp simd simdlen(bnl)
        for(int j=0; j<bnl; j++) {
            y[i+j] = x[i+j];
        }
        //最後の2つが対角ブロック
        for(int j=brptr[bidx]; j<brptr[bidx+1]; j++) {
            int x_ind = bcind[j]*bnw;
            for(int l=0; l<bnw; l++) {
                int off = j*b_size+l*bnl;
                #pragma omp simd simdlen(bnl)
                for(int k=0; k<bnl; k++) {
                    y[i+k] -= bval[off+k] * y[x_ind+l];
                }
            }
        }
    }
}

template <typename T, int bnl, int bnw> inline
void SptrsvBcsr_u(
    T *bval, int *bcind, int *brptr,
    T *x, T *y, int N)
{
    int b_size = bnl * bnw;
    int b_rem = bnl / bnw;
    for(int i=N-bnl; i>=0; i-=bnl) {
        int bidx = i / bnl;
        #pragma omp simd simdlen(bnl)
        for(int j=0; j<bnl; j++) {
            y[i+j] = x[i+j];
        }
        for(int j=brptr[bidx+1]-1; j>=brptr[bidx]+b_rem; j--) {
            int x_ind = bcind[j]*bnw;
            for(int l=0; l<bnw; l++) {
                int off = j*b_size+l*bnl;
                #pragma omp simd simdlen(bnl)
                for(int k=0; k<bnl; k++) {
                    y[i+k] -= bval[off+k] * y[x_ind+l];
                }
            }
        }
        int pos = brptr[bidx]+b_rem-1;
        for(int k=b_rem-1; k>=0; k--) {
            for(int j=bnw-1; j>=0; j--) {
                int off = pos*b_size+j*bnl;
                int idx = k*bnw+j;
                y[i+idx] *= bval[off+idx];
                for(int l=k*bnw+j-1; l>=0; l--) {
                    y[i+l] -= bval[off+l] * y[i+idx];
                }
            }
            pos--;
        }
    }
}

void SpmmCscCsc(
    double *l_val, int *l_rind, int *l_cptr,
    double *r_val, int *r_rind, int *r_cptr,
    double **val, int **rind, int **cptr,
    int L, int M, int R); // -> L x R matrix

// ---- integer ---- //

template <int bit>
void SpmvCsr(
    int *val, int *cind, int *rptr,
    int *x, int *y, int N)
{
    #pragma omp parallel for
    for(int i=0; i<N; i++) {
        long temp = 0;
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            temp += (long)val[j] * (long)x[cind[j]];
        }
        y[i] = (int)(temp >> bit);
    }
}

template <int bit>
void SpmvCsr(
    short *val, int *cind, int *rptr,
    int *x, int *y, int N)
{
    #pragma omp parallel for
    for(int i=0; i<N; i++) {
        long temp = 0;
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            temp += (long)val[j] * (long)x[cind[j]];
        }
        y[i] = (int)(temp >> bit);
    }
}

template <int bit>
void SptrsvCsr(
    int *val, int *cind, int *rptr,
    int *x, int *y, Triangle type, int N)
{
    switch (type) {
        case Upper:
            for(int i=N-1; i>=0; i--) {
                long temp = (long)x[i] << bit;
                int j;
                for(j=rptr[i+1]-1; j>=rptr[i]+1; j--) {
                    temp -= (long)val[j] * (long)y[cind[j]];    
                }
                y[i] = (int)((temp >> bit) * (long)val[j] >> bit);
            }
            break;
        case Lower:
            for(int i=0; i<N; i++) {
                long temp = (long)x[i] << bit;
                int j;
                for(j=rptr[i]; j<rptr[i+1]-1; j++) {
                    temp -= (long)val[j] * (long)y[cind[j]];
                }
                y[i] = (int)((temp >> bit) * (long)val[j] >> bit);
            }
            break;
        case UnitLower:
            for(int i=0; i<N; i++) {
                long temp = (long)x[i] << bit;
                int j;
                for(j=rptr[i]; j<rptr[i+1]; j++) {
                    temp -= (long)val[j] * (long)y[cind[j]];
                }
                y[i] = (int)(temp >> bit);
            }
            break;
        default:
            std::cerr << "SptrsvCsr: type is not valit." << std::endl;
            std::exit(1);
    }
}

}

}

#endif
