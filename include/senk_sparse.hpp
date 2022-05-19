/**
 * @file senk_sparse.hpp
 * @brief Functions related to sparse matrices are written.
 * @author Kengo Suzuki
 * @date 5/9/2022
 */
#ifndef SENK_SPARSE_HPP
#define SENK_SPARSE_HPP

namespace senk {
/**
 * @brief Functions related to sparse matrices and sparse vectors are defined.
 */
namespace sparse {
/**
 * @brief Perform SpMV using the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the CSR format.
 * @param cind A col-index array in the CSR format.
 * @param rptr A row-pointer array in the CSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 */
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
/**
 * @brief Perform SpMV using the CSR format, which stores diagonal elements separately.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the CSR format.
 * @param cind A col-index array in the CSR format.
 * @param rptr A row-pointer array in the CSR format.
 * @param diag An array that stores diagonal elements.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 */
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
/**
 * @brief Perform SpMV using the sliced-ELLPACK (SELL-c) format.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the SELL-c format.
 * @param cind A col-index array in the SELL-c format.
 * @param wid An array that indicates the starting position of the slices.
 * @param len The size of the slices.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 */
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
/**
 * @brief Perform the sparse lower triangular solve on a matrix stored in the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the CSR format.
 * @param cind A col-index array in the CSR format.
 * @param rptr A row-pointer array in the CSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 */
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
/**
 * @brief Perform the sparse upper triangular solve on a matrix stored in the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the CSR format.
 * @param cind A col-index array in the CSR format.
 * @param rptr A row-pointer array in the CSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 */
template <typename T> inline
void SptrsvCsr_u(T *val, int *cind, int *rptr, T *x, T *y, int N)
{
    // U is assumed to be general upper triangular.
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
/**
 * @brief Perform the sparse lower triangular solve in parallel on a AMC ordered matrix stored in the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the CSR format.
 * @param cind A col-index array in the CSR format.
 * @param rptr A row-pointer array in the CSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 */
template <typename T> inline
void SptrsvCsr_l(
    T *val, int *cind, int *rptr, T *x, T *y,
    int N, int *cptr, int cnum)
{
    // L is assumed to be unit lower triangular.
    #pragma omp parallel
    {
        for(int k=0; k<cnum; k++) {
            int start = cptr[k];
            int end = cptr[k+1];
            #pragma omp for
            for(int i=start; i<end; i++) {
                T temp = x[i];
                for(int j=rptr[i]; j<rptr[i+1]; j++) {
                    temp -= val[j] * y[cind[j]];
                }
                y[i] = temp;
            }
        } 
    }
}
/**
 * @brief Perform the sparse upper triangular solve in parallel on a AMC ordered matrix stored in the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the CSR format.
 * @param cind A col-index array in the CSR format.
 * @param rptr A row-pointer array in the CSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 */
template <typename T> inline
void SptrsvCsr_u(T *val, int *cind, int *rptr, T *x, T *y,
    int N, int *cptr, int cnum)
{
    // U is assumed to be general upper triangular.
    // Diagonal has been inverted.
    #pragma omp parallel
    {
        for(int k=cnum-1; k>=0; k--) {
            int start = cptr[k];
            int end = cptr[k+1];
            #pragma omp for
            for(int i=end-1; i>=start; i--) {
                T temp = x[i];
                int j;
                for(j=rptr[i+1]-1; j>=rptr[i]+1; j--) {
                    temp -= val[j] * y[cind[j]];    
                }
                y[i] = temp * val[j];
            }
        }
    }
}
/**
 * @brief Perform the sparse lower triangular solve in parallel on a ABMC ordered matrix stored in the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the CSR format.
 * @param cind A col-index array in the CSR format.
 * @param rptr A row-pointer array in the CSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 * @param bsize The number of rows/columns of the blocks used in ABMC.
 */
template <typename T> inline
void SptrsvCsr_l(
    T *val, int *cind, int *rptr, T *x, T *y,
    int N, int *cptr, int cnum, int bsize)
{
    // L is assumed to be unit lower triangular.
    #pragma omp parallel
    {
        for(int k=0; k<cnum; k++) {
            int start = cptr[k];
            int end = cptr[k+1];
            #pragma omp for
            for(int i=start; i<end; i++) {
                int base = i*bsize;
                for(int l=0; l<bsize; l++) {
                    int idx = base+l;
                    T temp = x[idx];
                    for(int j=rptr[idx]; j<rptr[idx+1]; j++) {
                        temp -= val[j] * y[cind[j]];
                    }
                    y[idx] = temp;
                }
            }
        } 
    }
}
/**
 * @brief Perform the sparse upper triangular solve in parallel on a ABMC ordered matrix stored in the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @param val A val array in the CSR format.
 * @param cind A col-index array in the CSR format.
 * @param rptr A row-pointer array in the CSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 * @param bsize The number of rows/columns of the blocks used in ABMC.
 */
template <typename T> inline
void SptrsvCsr_u(T *val, int *cind, int *rptr, T *x, T *y,
    int N, int *cptr, int cnum, int bsize)
{
    // U is assumed to be general upper triangular.
    // Diagonal has been inverted.
    #pragma omp parallel
    {
        for(int k=cnum-1; k>=0; k--) {
            int start = cptr[k];
            int end = cptr[k+1];
            #pragma omp for
            for(int i=end-1; i>=start; i--) {
                int base = i*bsize;
                for(int l=bsize-1; l>=0; l--) {
                    int idx = base+l;
                    T temp = x[idx];
                    int j;
                    for(j=rptr[idx+1]-1; j>=rptr[idx]+1; j--) {
                        temp -= val[j] * y[cind[j]];    
                    }
                    y[idx] = temp * val[j];
                }
            }
        }
    }
}
/**
 * @brief Perform the sparse lower triangular solve for a matrix stored in the BCSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @tparam bnl The number of rows of the block.
 * @tparam bnw The number of columns of the block.
 * @param bval A val array in the BCSR format.
 * @param bcind A block col-index array in the BCSR format.
 * @param brptr A block row-pointer array in the BCSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 */
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
/**
 * @brief Perform the sparse upper triangular solve for a matrix stored in the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @tparam bnl The number of rows of the block.
 * @tparam bnw The number of columns of the block.
 * @param bval A val array in the BCSR format.
 * @param bcind A block col-index array in the BCSR format.
 * @param brptr A block row-pointer array in the BCSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 */
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
/**
 * @brief Perform the sparse lower triangular solve for a ABMC reordered matrix stored in the BCSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @tparam bnl The number of rows of the block.
 * @tparam bnw The number of columns of the block.
 * @param bval A val array in the BCSR format.
 * @param bcind A block col-index array in the BCSR format.
 * @param brptr A block row-pointer array in the BCSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 * @param bsize The number of rows/columns of the blocks used in ABMC.
 */
template <typename T, int bnl, int bnw> inline
void SptrsvBcsr_l(
    T *bval, int *bcind, int *brptr, T *x, T *y,
    int N, int *cptr, int cnum, int bsize)
{
    // L is assumed to be unit lower triangular.
    int b_size = bnl * bnw;
    #pragma omp parallel
    {
        for(int k=0; k<cnum; k++) {
            int start = cptr[k];
            int end = cptr[k+1];
            #pragma omp for
            for(int i=start; i<end; i++) {
                int base = i*bsize;
                for(int l=0; l<bsize; l+=bnl) {
                    int idx = base+l;
                    int bidx = idx / bnl;
                    #pragma omp simd simdlen(bnl)
                    for(int j=0; j<bnl; j++) {
                        y[idx+j] = x[idx+j];
                    }
                    for(int j=brptr[bidx]; j<brptr[bidx+1]; j++) {
                        int x_ind = bcind[j]*bnw;
                        for(int m=0; m<bnw; m++) {
                            int off = j*b_size+m*bnl;
                            #pragma omp simd simdlen(bnl)
                            for(int n=0; n<bnl; n++) {
                                y[idx+n] -= bval[off+n] * y[x_ind+m];
                            }
                        }
                    }
                }
            }
        }
    }
}
/**
 * @brief Perform the sparse upper triangular solve for a ABMC reordered matrix stored in the CSR format.
 * @tparam T The Type of the matrix and the vectors.
 * @tparam bnl The number of rows of the block.
 * @tparam bnw The number of columns of the block.
 * @param bval A val array in the BCSR format.
 * @param bcind A block col-index array in the BCSR format.
 * @param brptr A block row-pointer array in the BCSR format.
 * @param x Input vector of size N.
 * @param y Output vector of size N.
 * @param N The size of vectors.
 * @param cptr The starting index of each color is stored.
 * @param cnum The number of colors.
 * @param bsize The number of rows/columns of the blocks used in ABMC.
 */
template <typename T, int bnl, int bnw> inline
void SptrsvBcsr_u(
    T *bval, int *bcind, int *brptr, T *x, T *y,
    int N, int *cptr, int cnum, int bsize)
{
    int b_size = bnl * bnw;
    int b_rem = bnl / bnw;
    #pragma omp parallel
    {
        for(int k=cnum-1; k>=0; k--) {
            int start = cptr[k];
            int end = cptr[k+1];
            #pragma omp for
            for(int i=end-1; i>=start; i--) {
                int base = i*bsize;
                for(int l=bsize-bnl; l>=0; l-=bnl) {
                    int idx = base+l;
                    int bidx = idx / bnl;
                    #pragma omp simd simdlen(bnl)
                    for(int j=0; j<bnl; j++) {
                        y[idx+j] = x[idx+j];
                    }
                    for(int j=brptr[bidx+1]-1; j>=brptr[bidx]+b_rem; j--) {
                        int x_ind = bcind[j]*bnw;
                        for(int n=0; n<bnw; n++) {
                            int off = j*b_size+n*bnl;
                            #pragma omp simd simdlen(bnl)
                            for(int m=0; m<bnl; m++) {
                                y[idx+m] -= bval[off+m] * y[x_ind+n];
                            }
                        }
                    }
                    int pos = brptr[bidx]+b_rem-1;
                    for(int m=b_rem-1; m>=0; m--) {
                        for(int j=bnw-1; j>=0; j--) {
                            int off = pos*b_size+j*bnl;
                            int ind = m*bnw+j;
                            y[idx+ind] *= bval[off+ind];
                            for(int n=m*bnw+j-1; n>=0; n--) {
                                y[idx+n] -= bval[off+n] * y[idx+ind];
                            }
                        }
                        pos--;
                    }
                }
            }
        }
    }
}
// ---- experimental ---- //
/*
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
*/
}

}

#endif
