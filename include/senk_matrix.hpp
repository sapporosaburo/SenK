#ifndef SENK_MATRIX_HPP
#define SENK_MATRIX_HPP

#include "senk_utils.hpp"

namespace senk {

namespace matrix {

void RemoveZeros(double **val, int **cind, int *rptr, int N);
bool CheckStructure(double *val, int *cind, int *rptr, int N);
void Padding(double **val, int **cind, int **rptr, int size, int *N);
void Split(
    double *val, int *cind, int *rptr,
    double **lval, int **lcind, int **lrptr,
    double **uval, int **ucind, int **urptr,
    double **diag, int N, const char *key, bool invDiag);
void Expand(
    double *lval, int *lcind, int *lrptr,
    double **val, int **cind, int **rptr,
    int N);
void Duplicate(
    double *tval, int *tcind, int *trptr,
    double **val, int **cind, int **rptr,
    int N);
void Scaling(double *val, int *cind, int *rptr, double *b, int N);
int Csr2PaddedCsr(
    double **val, int **cind, int **rptr,
    int size, int N);
int Csr2Sell(
    double *val, int *cind, int *rptr,
    double **s_val, int **s_cind, int **s_len,
    int size, int N);
void Csr2Csc(
    double *val, int *cind, int *rptr,
    double **cval, int **crind, int **ccptr,
    int N, int M);
void Ilu0(double *val, int *cind, int *rptr, int N);
void Ilup(double **val, int **cind, int **rptr, int N, int p);

template <int bnl, int bnw>
int Csr2Bcsr(
    double *val, int *cind, int *rptr,
    double **bval, int **bcind, int **brptr, int N)
{
    if(N % bnl != 0 || N % bnw != 0) {
        printf("Error: Csr2Bcsr\n");
        exit(EXIT_FAILURE);
    }
    int i, j;
    int *t_rptr = senk::utils::SafeMalloc<int>(bnl);
    int num_block = 0;
    *brptr = senk::utils::SafeMalloc<int>(N/bnl+1);
    (*brptr)[0] = 0;
    // Search the number of block
    for(int i=0; i<N; i+=bnl) {
        for(j=0; j<bnl; j++) {
            t_rptr[j] = (rptr[i+j] != rptr[i+j+1])? rptr[i+j] : -1;
        }
        int now_bcol = -1;
        int flag = 0;
        while(flag == 0) {
            int min_bcol = N;
            for(j=0; j<bnl; j++) {
                if(t_rptr[j] == -1) continue;
                if(cind[t_rptr[j]]/bnw < min_bcol) {
                    min_bcol = cind[t_rptr[j]]/bnw;
                }
            }
            if(min_bcol > now_bcol) {
                now_bcol = min_bcol; num_block++;
            }
            for(j=0; j<bnl; j++) {
                if( t_rptr[j] == -1 ) continue;
                if( cind[t_rptr[j]]/bnw == min_bcol ) {
                    t_rptr[j]++;
                    if( t_rptr[j] == rptr[i+j+1] ) {
                        t_rptr[j] = -1;
                    }
                }
            }
            flag = 1;
            for(j=0; j<bnl; j++) {
                if(t_rptr[j] != -1) flag = 0;
            }
        }
        (*brptr)[i/bnl+1] = num_block;
    }
    *bcind = senk::utils::SafeMalloc<int>(num_block);
    *bval = senk::utils::SafeCalloc<double>(num_block*bnl*bnw);
    num_block = 0;
    // Assign val to bval
    int bval_pos = -bnl*bnw;
    for(i=0; i<N; i+=bnl) {
        for(j=0; j<bnl; j++) {
            t_rptr[j] = rptr[i+j];
            if( t_rptr[j] == rptr[i+j+1] ) {
                t_rptr[j] = -1;
            }
        }
        int now_bcol = -1;
        int flag = 0;
        while(flag == 0) {
            int min_col = N;
            for(j=0; j<bnl; j++) {
                if( t_rptr[j] == -1 ) continue;
                if( cind[t_rptr[j]] < min_col ) {
                    min_col = cind[t_rptr[j]];
                }
            }
            if(min_col/bnw > now_bcol) {
                bval_pos += bnl * bnw;
                now_bcol = min_col/bnw;
                (*bcind)[num_block] = now_bcol;
                num_block++;
            }
            int off = min_col % bnw;
            for(j=0; j<bnl; j++) {
                if( t_rptr[j] == -1 ) continue;
                if( cind[t_rptr[j]] == min_col ) {
                    (*bval)[bval_pos+off*bnl+j] = val[t_rptr[j]];
                    t_rptr[j]++;
                    if( t_rptr[j] == rptr[i+j+1] ) {
                        t_rptr[j] = -1;
                    }
                }
            }
            flag = 1;
            for(j=0; j<bnl; j++) {
                if(t_rptr[j] != -1) flag = 0;
            }
        }
    }
    free(t_rptr);
    return num_block;
}

template <int bnl, int bnw>
int Bcsr2Csr(
    double *bval, int *bcind, int *brptr,
    double **val, int **cind, int **rptr, int N)
{
    int bsize = bnl * bnw;
    int num_block = brptr[N/bnl];
    int nnz;
    *val = senk::utils::SafeMalloc<double>(num_block*bnl*bnw);
    *cind = senk::utils::SafeMalloc<int>(num_block*bnl*bnw);
    *rptr = senk::utils::SafeMalloc<int>(N+1);
    (*rptr)[0] = 0;
    int count = 0;
    for(int i=0; i<N; i++) {
        int bid = i / bnl;
        int id = i % bnl; // 0 to bnl-1
        for(int bj=brptr[bid]; bj<brptr[bid+1]; bj++) {
            for(int j=0; j<bnw; j++) {
                (*cind)[count] = (bcind[bj])*bnw+j;
                (*val)[count] = bval[bj*bsize+j*bnl+id];
                count++;
            }
        }
        (*rptr)[i+1] = count;
    }
    nnz = (*rptr)[N];
    return nnz;
}

}

}

#endif
