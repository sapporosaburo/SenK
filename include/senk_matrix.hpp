/**
 * @file senk_matrix.hpp
 * @brief Functions related (sparse) matrices are defined.
 * @author Kengo Suzuki
 * @date 5/9/2022
 */
#ifndef SENK_MATRIX_HPP
#define SENK_MATRIX_HPP

#include <cstdio>
#include <cstring>

#include "senk_utils.hpp"
#include "senk_class.hpp"

namespace senk {
/**
 * @brief Functions related to matrix are located in namespace matrix.
 */
namespace matrix {
/**
 * @brief Remove zero elements form an input matrix.
 * @tparam The Type of the input matrix.
 * @param val The pointer to the array storing val in the CSR format.
 * @param cind The pointer to the array storing column index in the CSR format.
 * @param rptr The pointer to the array storing row pointer in the CSR format.
 * @param N The size of the input matrix.
 */
template <typename T>
void RemoveZeros(T **val, int **cind, int *rptr, int N)
{
    int off = 0;
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]+off; j<rptr[i+1]; j++) {
            if((*val)[j] == 0) { off++; continue; }
            (*val)[j-off] = (*val)[j];
            (*cind)[j-off] = (*cind)[j];
        }
        rptr[i+1] -= off;
    }
    *val = utils::SafeRealloc<T>(*val, rptr[N]);
    *cind = utils::SafeRealloc<int>(*cind, rptr[N]);
}

template <typename T>
bool CheckStructure(T *val, int *cind, int *rptr, int N)
{
    for(int i=0; i<N; i++) {
        bool flag = false;
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            if(cind[j] == i) { flag = true; break; }
        }
        if(!flag) { return false; }
    }
    return true;
}

template <typename T>
void Csr2Csc(
    T *val, int *cind, int *rptr,
    T **cval, int **crind, int **ccptr,
    int N, int M)
{
    // The size of the input matrix is N * M
    int nnz = rptr[N];
    *cval  = utils::SafeMalloc<T>(nnz);
    *crind = utils::SafeMalloc<int>(nnz);
    *ccptr = utils::SafeMalloc<int>(M+1);
    int *num = utils::SafeCalloc<int>(M);
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]; j<rptr[i+1]; j++) { num[cind[j]]++; }
    }
    (*ccptr)[0] = 0;
    for(int i=0; i<M; i++) {
        (*ccptr)[i+1] = (*ccptr)[i] + num[i];
        num[i] = 0;
    }
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            int off = (*ccptr)[cind[j]];
            int pos = num[cind[j]];
            (*crind)[off+pos] = i;
            (*cval)[off+pos] = val[j];
            num[cind[j]]++;
        }
    }
    free(num);
}

template <typename T>
int Csr2Sell(
    T *val, int *cind, int *rptr,
    T **s_val, int **s_cind, int **s_len,
    int size, int N)
{
    int num_slice = (N+size-1)/size;
    *s_len = utils::SafeMalloc<int>(num_slice+1);
    int temp_slice;
    int row_id;
    int row_max;
    int row_count;
    int nnz = 0;
    
    (*s_len)[0] = 0;
    for(int i=0; i<num_slice; i++) {
        temp_slice = size;
        if(i==num_slice-1 && N%size!=0) temp_slice = N % size;
        row_max = 0;
        for(int j=0; j<temp_slice; j++) {
            row_id = i*size+j;
            if(row_max < rptr[row_id+1] - rptr[row_id]) {
                row_max = rptr[row_id+1] - rptr[row_id];
            }
        }
        (*s_len)[i+1] = (*s_len)[i] + row_max;
        nnz += row_max * temp_slice;
    }
    *s_val  = utils::SafeMalloc<T>(nnz);
    *s_cind = utils::SafeMalloc<int>(nnz);
    for(int i=0; i<num_slice; i++) {
        temp_slice = size;
        if(i==num_slice-1 && N%size!=0) temp_slice = N % size;
        for(int j=0; j<temp_slice; j++) {
            row_id = i*size+j;
            row_count = 0;
            for(int k=rptr[row_id]; k<rptr[row_id+1]; k++) {
                (*s_val)[(*s_len)[i]*size+row_count*temp_slice+j] = val[k];
                (*s_cind)[(*s_len)[i]*size+row_count*temp_slice+j] = cind[k];
                row_count++;
            }
            if(row_count < (*s_len)[i+1] - (*s_len)[i]) {
                int temp = (*s_len)[i+1] - (*s_len)[i] - row_count;
                for(int k=0; k<temp; k++) {
                    (*s_val)[(*s_len)[i]*size+row_count*temp_slice+j] = 0;
                    (*s_cind)[(*s_len)[i]*size+row_count*temp_slice+j] = 0;
                    row_count++;
                }
            }
        }
    }
    return nnz;
}

template <typename T>
void Csr2Bcsr(
    T *val, int *cind, int *rptr,
    T **bval, int **bcind, int **brptr,
    int N, int bnl, int bnw)
{
    if(N % bnl || N % bnw) {
        printf("Error: Csr2Bcsr\n");
        exit(EXIT_FAILURE);
    }
    *brptr = utils::SafeMalloc<int>(N/bnl+1);
    int cnt = 0;
    (*brptr)[0] = 0;
    int *ptr = utils::SafeMalloc<int>(bnl);
// Count the number of block
    for(int i=0; i<N; i+=bnl) {
        // Initialize "ptr"
        for(int j=0; j<bnl; j++) { ptr[j] = rptr[i+j]; }
        while(true) {
            // Find minimum col value
            int min = N;
            for(int j=0; j<bnl; j++) {
                if(ptr[j] != -1 && cind[ptr[j]] < min)
                    min = cind[ptr[j]];
            }
            if(min == N) break;
            // Increment ptr[j] whose col-idx is in min block
            for(int j=0; j<bnl; j++) {
                if(ptr[j] == -1) continue;
                while(cind[ptr[j]]/bnw == min/bnw) {
                    ptr[j]++;
                    if(ptr[j] >= rptr[i+j+1]) {ptr[j] = -1; break;}
                }
            }
            cnt++;
        }
        (*brptr)[i/bnl+1] = cnt;
    }
    *bcind = utils::SafeMalloc<int>(cnt);
    *bval  = utils::SafeCalloc<T>(cnt*bnl*bnw);
// Assign val to bval
    cnt = 0;
    for(int i=0; i<N; i+=bnl) {
        // Initialize "ptr"
        for(int j=0; j<bnl; j++) { ptr[j] = rptr[i+j]; }
        while(true) {
            int min = N;
            for(int j=0; j<bnl; j++) {
                if(ptr[j] != -1 && cind[ptr[j]] < min)
                    min = cind[ptr[j]];
            }
            if(min == N) break;
            for(int j=0; j<bnl; j++) {
                if(ptr[j] == -1) continue;
                while(cind[ptr[j]]/bnw == min/bnw) {
                    int off = cind[ptr[j]] % bnw;
                    (*bval)[cnt*bnl*bnw+off*bnl+j] = val[ptr[j]];
                    //printf("%e\n", (*bval)[cnt*bnl*bnw+off*bnl+j]);
                    ptr[j]++;
                    if(ptr[j] >= rptr[i+j+1]) {ptr[j] = -1; break;}
                }
            }
            (*bcind)[cnt] = min/bnw;
            cnt++;
        }
    }
    free(ptr);
}

template <typename T>
int Bcsr2Csr(
    T *bval, int *bcind, int *brptr,
    T **val, int **cind, int **rptr,
    int N, int bnl, int bnw)
{
    int bsize = bnl * bnw;
    int num_block = brptr[N/bnl];
    int nnz;
    *val = senk::utils::SafeMalloc<T>(num_block*bnl*bnw);
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

template <typename T>
void Padding(T **val, int **cind, int **rptr, int size, int *N)
{
    int remain = (N[0] % size == 0) ? 0 : size - N[0] % size;
    int NNZ = rptr[0][N[0]];
    *val  = utils::SafeRealloc<T>(*val, NNZ+remain);
    *cind = utils::SafeRealloc<int>(*cind, NNZ+remain);
    *rptr = utils::SafeRealloc<int>(*rptr, N[0]+1+remain);
    for(int i=0; i<remain; i++) val[0][NNZ+i] = 1;
    for(int i=0; i<remain; i++) cind[0][NNZ+i] = N[0]+i;
    for(int i=0; i<remain; i++) rptr[0][N[0]+1+i] = rptr[0][N[0]+i]+1;
    N[0] += remain;
}

template <typename T>
void Split(
    T *val, int *cind, int *rptr,
    T **lval, int **lcind, int **lrptr,
    T **uval, int **ucind, int **urptr,
    T **diag, int N, const char *key, bool invDiag)
{
    int lNNZ = 0;
    int uNNZ = 0;
    int dNNZ = 0;
    int *l_ptr, *u_ptr, *d_ptr;
    T *lv_ptr, *uv_ptr, *dv_ptr;
    int *lc_ptr, *uc_ptr, *dc_ptr;
    int *lr_ptr, *ur_ptr;
    if(std::strcmp(key, "LD-U") == 0) {
        l_ptr=&lNNZ; u_ptr=&uNNZ; d_ptr=&lNNZ;
    }else if(std::strcmp(key, "L-DU") == 0) {
        l_ptr=&lNNZ; u_ptr=&uNNZ; d_ptr=&uNNZ;
    }else if(std::strcmp(key, "L-D-U") == 0) {
        l_ptr=&lNNZ; u_ptr=&uNNZ; d_ptr=&dNNZ;
    }else if(std::strcmp(key, "LU-D") == 0) {
        l_ptr=&lNNZ; u_ptr=&lNNZ; d_ptr=&dNNZ;
    }else { printf("Split: Keyword is not valid."); exit(1); }
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            if(cind[j] < i) { (*l_ptr)++; }
            else if(cind[j] > i) { (*u_ptr)++; }
            else { (*d_ptr)++; }
        }
    }
    *lval  = utils::SafeMalloc<T>(lNNZ);
    *lcind = utils::SafeMalloc<int>(lNNZ);
    *lrptr = utils::SafeMalloc<int>(N+1);
    if(uNNZ != 0) {
        *uval  = utils::SafeMalloc<T>(uNNZ);
        *ucind = utils::SafeMalloc<int>(uNNZ);
        *urptr = utils::SafeMalloc<int>(N+1);
    }
    if(dNNZ == N) { *diag = utils::SafeMalloc<T>(dNNZ); }
    lNNZ = 0; uNNZ = 0; dNNZ = 0;
    if(std::strcmp(key, "LD-U") == 0) {
        lv_ptr = *lval;  uv_ptr = *uval;  dv_ptr = *lval;
        lc_ptr = *lcind; uc_ptr = *ucind; dc_ptr = *lcind;
        lr_ptr = *lrptr; ur_ptr = *urptr;
    }else if(std::strcmp(key, "L-DU") == 0) {
        lv_ptr = *lval;  uv_ptr = *uval;  dv_ptr = *uval;
        lc_ptr = *lcind; uc_ptr = *ucind; dc_ptr = *ucind;
        lr_ptr = *lrptr; ur_ptr = *urptr;
    }else if(std::strcmp(key, "L-D-U") == 0) {
        lv_ptr = *lval;  uv_ptr = *uval;  dv_ptr = *diag;
        lc_ptr = *lcind; uc_ptr = *ucind; dc_ptr = nullptr;
        lr_ptr = *lrptr; ur_ptr = *urptr;
    }else if(std::strcmp(key, "LU-D") == 0) {
        lv_ptr = *lval;  uv_ptr = *lval;  dv_ptr = *diag;
        lc_ptr = *lcind; uc_ptr = *lcind; dc_ptr = nullptr;
        lr_ptr = *lrptr; ur_ptr = *lrptr;
    }else { return; }
    lr_ptr[0] = 0;
    ur_ptr[0] = 0;
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            if(cind[j] < i) {
                lv_ptr[*l_ptr] = val[j]; lc_ptr[*l_ptr] = cind[j];
                (*l_ptr)++;
            }else if(cind[j] > i) {
                uv_ptr[*u_ptr] = val[j]; uc_ptr[*u_ptr] = cind[j];
                (*u_ptr)++;
            }else {
                dv_ptr[*d_ptr] = (invDiag) ? 1/val[j] : val[j];
                if(dc_ptr) { dc_ptr[*d_ptr] = cind[j]; }
                (*d_ptr)++;
            }
        }
        lr_ptr[i+1] = *l_ptr;
        ur_ptr[i+1] = *u_ptr;
    }
}

template <typename T>
void Expand(
    T *tval, int *tcind, int *trptr,
    T **val, int **cind, int **rptr,
    int N, const char *key)
{
    T   *lval,  *uval,  *tval2;
    int *lcind, *ucind, *tcind2;
    int *lrptr, *urptr, *trptr2;
    int nnz = trptr[N]*2 - N;
    Csr2Csc<T>(tval, tcind, trptr, &tval2, &tcind2, &trptr2, N, N);
    if(std::strcmp(key, "L") == 0) {
        lval = tval;  lcind = tcind;  lrptr = trptr;
        uval = tval2; ucind = tcind2; urptr = trptr2;
    }else if(std::strcmp(key, "U") == 0) {
        lval = tval2; lcind = tcind2; lrptr = trptr2;
        uval = tval;  ucind = tcind;  urptr = trptr;
    }else { printf("Expand: Keyword is not valid."); exit(1); }
    *val  = utils::SafeMalloc<T>(nnz);
    *cind = utils::SafeMalloc<int>(nnz);
    *rptr = utils::SafeMalloc<int>(N+1);
    int cnt = 0; (*rptr)[0] = cnt;
    for(int i=0; i<N; i++) {
        for(int j=lrptr[i]; j<lrptr[i+1]; j++) {
            (*val)[cnt] = lval[j]; (*cind)[cnt] = lcind[j];
            cnt++;
        }
        for(int j=urptr[i]+1; j<urptr[i+1]; j++) {
            (*val)[cnt] = uval[j]; (*cind)[cnt] = ucind[j];
            cnt++;
        }
        (*rptr)[i+1] = cnt;
    }
    free(tval2);
    free(tcind2);
    free(trptr2);
}

template <typename T>
void Duplicate(
    T *tval, int *tcind, int *trptr,
    T **val, int **cind, int **rptr,
    int N)
{
    *val  = utils::SafeMalloc<T>(trptr[N]);
    *cind = utils::SafeMalloc<int>(trptr[N]);
    *rptr = utils::SafeMalloc<int>(N+1);
    (*rptr)[0] = trptr[0];
    for(int i=0; i<N; i++) {
        (*rptr)[i+1] = trptr[i+1];
        for(int j=trptr[i]; j<trptr[i+1]; j++) {
            (*val)[j] = tval[j]; (*cind)[j] = tcind[j];
        }
    }
}

template <typename T>
void GetDiag(T *val, int *cind, int *rptr, T **diag, int N)
{
    *diag = utils::SafeMalloc<T>(N);
    #pragma omp parallel for
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            if(cind[j] == i) {
                (*diag)[i] = val[j];
                break;
            }
        }
    }
}

template <typename T>
void Scaling(T *val, int *cind, int *rptr, T *ld, T *rd, int N)
{
    #pragma omp parallel for
    for(int i=0; i<N; i++) {
        T left = ld[i];
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            val[j] = left * val[j] * rd[cind[j]];
        }
    }
}

template <typename T>
void Ilu0(T *val, int *cind, int *rptr, int N)
{
    for(int i=1; i<N; i++) {
        for(int k=rptr[i]; k<rptr[i+1]; k++) {
            if(cind[k] >= i) break;
            for(int l=rptr[cind[k]]; l<rptr[cind[k]+1]; l++) {
                if(cind[l] == cind[k]) {
                    if(val[l] == 0) {
                        printf("Error: Ilu0, 0 pivot\n");
                        exit(1);
                    }
                    val[k] = val[k] / val[l];
                    break;
                }
            }
            int pos = rptr[cind[k]];
            for(int j=k+1; j<rptr[i+1]; j++) {
                for(int l=pos; l<rptr[cind[k]+1]; l++) {
                    if(cind[l] < cind[j]) continue;
                    pos = l;
                    if(cind[l] == cind[j]) {
                        val[j] -= val[k] * val[l];
                        pos++;
                    }
                    break;
                }
            }
        }
    }
}

template <typename T>
void Ilup(T **val, int **cind, int **rptr, int N, int p)
{
    int i, j;
    //int NNZ = (*rptr)[N];
    sparse::SpVec<T, int> temp;
    sparse::SpVec<T, int> temp2;
    int *zeros = utils::SafeCalloc<int>(N);

    T pivot = 0;
    int pivot_lev = 0;

    int first_len = (*rptr)[1] - (*rptr)[0];
    T   *new_val  = utils::SafeMalloc<T>(first_len);
    int *new_cind = utils::SafeMalloc<int>(first_len);
    int *new_lev  = utils::SafeCalloc<int>(first_len);
    int *new_rptr = utils::SafeMalloc<int>(N+1);
    int new_len = first_len;
    //一行目
    memcpy(new_val, *val, sizeof(T)*first_len);
    memcpy(new_cind, *cind, sizeof(int)*first_len);
    new_rptr[0] = 0;
    new_rptr[1] = first_len;

    for(i=1; i<N; i++) {
        int now_len = (*rptr)[i+1] - (*rptr)[i];
        temp.Clear();
        temp.Append(&(*val)[(*rptr)[i]], zeros, &(*cind)[(*rptr)[i]], now_len);
        int count = 0;
        while(count < temp.GetLen() && temp.GetIdx(count) < i) {
            int k_ptr = count;
            if(temp.GetLev(k_ptr) > p) { count++; continue; }
            int k = temp.GetIdx(k_ptr);
            temp2.Clear();
            for(j=0; j<count; j++) {
                temp2.Append(temp.GetVal(j), temp.GetLev(j), temp.GetIdx(j));
            }
            for(j=new_rptr[k]; j<new_rptr[k+1]; j++) {
                if(new_cind[j] == k) {
                    pivot = temp.GetVal(k_ptr) / new_val[j];
                    pivot_lev = temp.GetLev(k_ptr);
                    break;
                }
            }
            j++;
            k_ptr++;
            temp2.Append(pivot, 0, k);
            while(k_ptr<temp.GetLen() || j<new_rptr[k+1]) {
                int t_col1 = N+1;
                int t_col2 = N+1;
                if(k_ptr < temp.GetLen()) t_col1 = temp.GetIdx(k_ptr);
                if(j < new_rptr[k+1]) t_col2 = new_cind[j];
                if(t_col1<t_col2) {
                    T t_val = temp.GetVal(k_ptr);
                    int t_lev    = temp.GetLev(k_ptr);
                    int t_col    = temp.GetIdx(k_ptr);
                    temp2.Append(t_val, t_lev, t_col);
                    k_ptr++;
                }else if(t_col1==t_col2) {
                    T t_val = temp.GetVal(k_ptr) - pivot*new_val[j];
                    int t_lev    = temp.GetLev(k_ptr);
                    int t_col    = temp.GetIdx(k_ptr);
                    if(t_lev > pivot_lev+new_lev[j]+1) {
                        t_lev = pivot_lev+new_lev[j]+1;
                    }
                    temp2.Append(t_val, t_lev, t_col);
                    k_ptr++;
                    j++;
                }else { // (k>j) fill-in
                    T t_val = -pivot*new_val[j];
                    int t_lev    = pivot_lev+new_lev[j]+1;
                    int t_col    = new_cind[j];
                    //printf("lev %d\n", t_lev);
                    temp2.Append(t_val, t_lev, t_col);
                    j++;
                }
            }
            count++;
            temp.Clear();
            for(j=0; j<temp2.GetLen(); j++) {
                temp.Append(temp2.GetVal(j), temp2.GetLev(j), temp2.GetIdx(j));
            }
        }
        new_val  = utils::SafeRealloc<T>(new_val, new_len+temp.GetLen());
        new_cind = utils::SafeRealloc<int>(new_cind, new_len+temp.GetLen());
        new_lev  = utils::SafeRealloc<int>(new_lev, new_len+temp.GetLen());
        for(j=0; j<temp.GetLen(); j++) {
            if(temp.GetLev(j) <= p) {
                new_val[new_len]  = temp.GetVal(j);
                new_cind[new_len] = temp.GetIdx(j);
                new_lev[new_len]  = temp.GetLev(j);
                new_len++;
            }
        }
        new_rptr[i+1] = new_len;
    }
    (*val)  = utils::SafeRealloc<T>(new_val, new_len);
    (*cind) = utils::SafeRealloc<int>(new_cind, new_len);
    (*rptr) = utils::SafeRealloc<int>(new_rptr, N+1);
}

template <typename T>
void AllocLevelZero(T **val, int **cind, int **rptr, int N, int p)
{
    int i, j;
    //int NNZ = (*rptr)[N];
    sparse::SpVec<T, int> temp;
    sparse::SpVec<T, int> temp2;
    int *zeros = utils::SafeCalloc<int>(N);

    T pivot = 0;
    int pivot_lev = 0;

    int first_len = (*rptr)[1] - (*rptr)[0];
    T   *new_val  = utils::SafeMalloc<T>(first_len);
    int *new_cind = utils::SafeMalloc<int>(first_len);
    int *new_lev  = utils::SafeCalloc<int>(first_len);
    int *new_rptr = utils::SafeMalloc<int>(N+1);
    int new_len = first_len;
    //一行目
    memcpy(new_val, *val, sizeof(T)*first_len);
    memcpy(new_cind, *cind, sizeof(int)*first_len);
    new_rptr[0] = 0;
    new_rptr[1] = first_len;

    for(i=1; i<N; i++) {
        int now_len = (*rptr)[i+1] - (*rptr)[i];
        temp.Clear();
        temp.Append(&(*val)[(*rptr)[i]], zeros, &(*cind)[(*rptr)[i]], now_len);
        int count = 0;
        while(count < temp.GetLen() && temp.GetIdx(count) < i) {
            int k_ptr = count;
            if(temp.GetLev(k_ptr) > p) { count++; continue; }
            int k = temp.GetIdx(k_ptr);
            temp2.Clear();
            for(j=0; j<count; j++) {
                temp2.Append(temp.GetVal(j), temp.GetLev(j), temp.GetIdx(j));
            }
            for(j=new_rptr[k]; j<new_rptr[k+1]; j++) {
                if(new_cind[j] == k) {
                    pivot = temp.GetVal(k_ptr);// / new_val[j];
                    pivot_lev = temp.GetLev(k_ptr);
                    break;
                }
            }
            j++;
            k_ptr++;
            temp2.Append(pivot, 0, k);
            while(k_ptr<temp.GetLen() || j<new_rptr[k+1]) {
                int t_col1 = N+1;
                int t_col2 = N+1;
                if(k_ptr < temp.GetLen()) t_col1 = temp.GetIdx(k_ptr);
                if(j < new_rptr[k+1]) t_col2 = new_cind[j];
                if(t_col1<t_col2) {
                    T t_val = temp.GetVal(k_ptr);
                    int t_lev    = temp.GetLev(k_ptr);
                    int t_col    = temp.GetIdx(k_ptr);
                    temp2.Append(t_val, t_lev, t_col);
                    k_ptr++;
                }else if(t_col1==t_col2) {
                    T t_val = temp.GetVal(k_ptr);
                    int t_lev    = temp.GetLev(k_ptr);
                    int t_col    = temp.GetIdx(k_ptr);
                    if(t_lev > pivot_lev+new_lev[j]+1) {
                        t_lev = pivot_lev+new_lev[j]+1;
                    }
                    temp2.Append(t_val, t_lev, t_col);
                    k_ptr++;
                    j++;
                }else { // (k>j) fill-in
                    T t_val = 0;
                    int t_lev    = pivot_lev+new_lev[j]+1;
                    int t_col    = new_cind[j];
                    //printf("lev %d\n", t_lev);
                    temp2.Append(t_val, t_lev, t_col);
                    j++;
                }
            }
            count++;
            temp.Clear();
            for(j=0; j<temp2.GetLen(); j++) {
                temp.Append(temp2.GetVal(j), temp2.GetLev(j), temp2.GetIdx(j));
            }
        }
        new_val  = utils::SafeRealloc<T>(new_val, new_len+temp.GetLen());
        new_cind = utils::SafeRealloc<int>(new_cind, new_len+temp.GetLen());
        new_lev  = utils::SafeRealloc<int>(new_lev, new_len+temp.GetLen());
        for(j=0; j<temp.GetLen(); j++) {
            if(temp.GetLev(j) <= p) {
                new_val[new_len]  = temp.GetVal(j);
                new_cind[new_len] = temp.GetIdx(j);
                new_lev[new_len]  = temp.GetLev(j);
                new_len++;
            }
        }
        new_rptr[i+1] = new_len;
    }
    (*val)  = utils::SafeRealloc<T>(new_val, new_len);
    (*cind) = utils::SafeRealloc<int>(new_cind, new_len);
    (*rptr) = utils::SafeRealloc<int>(new_rptr, N+1);
}

template <typename T>
void AllocBlockZero(T **val, int **cind, int **rptr, int N, int bnl, int bnw)
{
    T *bval;
    int *bcind;
    int *brptr;
    Csr2Bcsr<T>(*val, *cind, *rptr, &bval, &bcind, &brptr, N, bnl, bnw);
    free(*val);
    free(*cind);
    free(*rptr);
    Bcsr2Csr<T>(bval, bcind, brptr, val, cind, rptr, N, bnl, bnw);
    free(bval);
    free(bcind);
    free(brptr);
}

template <typename T>
void RemoveOffDiagonal(T **val, int **cind, int **rptr, int N, int bnum)
{
    int bsize = N / bnum;
    int off = 0;
    for(int i=0; i<N; i++) {
        int start = (*rptr)[i]+off;
        int min = i / bsize * bsize;
        int max = min + bsize;
        for(int j=start; j<(*rptr)[i+1]; j++) {
            if((*cind)[j] < min || max <= (*cind)[j]) {
                off++;
                continue;
            }
            (*val)[j-off] = (*val)[j];
            (*cind)[j-off] = (*cind)[j];
        }
        (*rptr)[i+1] -= off;
    }
    *val = utils::SafeRealloc<double>(*val, (*rptr)[N]);
    *cind = utils::SafeRealloc<int>(*cind, (*rptr)[N]);
}

template <typename T>
void Reordering(T *val, int *cind, int *rptr, int *RP, int N)
{
    #pragma omp parallel for
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            cind[j] = RP[cind[j]];
        }
        helper::QuickSort<int, double>(cind, val, rptr[i], rptr[i+1]-1);
    }
}

template <typename T>
void Reordering(T *val, int *cind, int *rptr, int *LP, int *RP, int N)
{
    int NNZ = rptr[N];
    double *tval = utils::SafeMalloc<double>(NNZ);
    int *tcind = utils::SafeMalloc<int>(NNZ);
    int *trptr = utils::SafeMalloc<int>(N+1);
    trptr[0] = 0;
    for(int i=0; i<N; i++) {
        int id = LP[i];
        trptr[i+1] = trptr[i] + (rptr[id+1] - rptr[id]);
    }
    #pragma omp parallel for
    for(int i=0; i<N; i++) {
        int id = LP[i];
        for(int j=0; j<rptr[id+1]-rptr[id]; j++) {
            tval[trptr[i]+j] = val[rptr[id]+j];
            tcind[trptr[i]+j] = RP[cind[rptr[id]+j]];
        }
        helper::QuickSort<int, double>(tcind, tval, trptr[i], trptr[i+1]-1);
    }
    utils::Copy<double>(tval, val, NNZ);
    utils::Copy<int>(tcind, cind, NNZ);
    utils::Copy<int>(trptr, rptr, N+1);
    free(tval);
    free(tcind);
    free(trptr);
}

/*
int Csr2PaddedCsr(
    double **val, int **cind, int **rptr,
    int size, int N)
{
    int nnz = 0;
    int *t_rptr = new int[N+1];
    t_rptr[0] = 0;
    for(int i=0; i<N; i++) {
        int num = rptr[0][i+1] - rptr[0][i];
        nnz += (num+size-1)/size * size;
        t_rptr[i+1] = nnz;
    }
    double *t_val = new double[nnz];
    int *t_cind = new int[nnz];
    for(int i=0; i<N; i++) {
        int num = rptr[0][i+1] - rptr[0][i];
        int cnt = 0;
        for(int j=t_rptr[i]; j<t_rptr[i+1]; j++) {
            if(cnt < num) {
                t_val[j] = val[0][rptr[0][i]+cnt];
                t_cind[j] = cind[0][rptr[0][i]+cnt];
            }else {
                t_val[j] = 0;
                t_cind[j] = t_cind[j-1];
            }
            cnt++;
        }
    }
    *val = (double*)realloc(*val, sizeof(double)*nnz);
    *cind = (int*)realloc(*cind, sizeof(int)*nnz);
    for(int i=0; i<nnz; i++) {val[0][i] = t_val[i];}
    for(int i=0; i<nnz; i++) {cind[0][i] = t_cind[i];}
    for(int i=0; i<N+1; i++) {rptr[0][i] = t_rptr[i];}
    
    return nnz;
}
*/
} // namespace matrix

} // namespace senk

#endif
