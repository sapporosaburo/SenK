#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>

#ifdef __GNUC__
#include <omp.h>
#endif

#include "senk_matrix.hpp"
#include "senk_utils.hpp"

namespace senk {

namespace matrix {

using senk::utils::SafeMalloc;
using senk::utils::SafeRealloc;

void RemoveZeros(double **val, int **cind, int *rptr, int N)
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
    *val = SafeRealloc<double>(*val, rptr[N]);
    *cind = SafeRealloc<int>(*cind, rptr[N]);
}

bool CheckStructure(double *val, int *cind, int *rptr, int N)
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

void Padding(double **val, int **cind, int **rptr, int size, int *N)
{
    int remain = (N[0] % size == 0) ? 0 : size - N[0] % size;
    int NNZ = rptr[0][N[0]];
    *val  = SafeRealloc<double>(*val, NNZ+remain);
    *cind = SafeRealloc<int>(*cind, NNZ+remain);
    *rptr = SafeRealloc<int>(*rptr, N[0]+1+remain);
    for(int i=0; i<remain; i++) val[0][NNZ+i] = 1;
    for(int i=0; i<remain; i++) cind[0][NNZ+i] = N[0]+i;
    for(int i=0; i<remain; i++) rptr[0][N[0]+1+i] = rptr[0][N[0]+i]+1;
    N[0] += remain;
}

void Split(
    double *val, int *cind, int *rptr,
    double **lval, int **lcind, int **lrptr,
    double **uval, int **ucind, int **urptr,
    double **diag, int N, const char *key, bool invDiag)
{
    int lNNZ = 0;
    int uNNZ = 0;
    int dNNZ = 0;
    int *l_ptr, *u_ptr, *d_ptr;
    double *lv_ptr, *uv_ptr, *dv_ptr;
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
    }else { printf("Split: Keyword is not valid."); return; }
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            if(cind[j] < i) { (*l_ptr)++; }
            else if(cind[j] > i) { (*u_ptr)++; }
            else { (*d_ptr)++; }
        }
    }
    *lval  = SafeMalloc<double>(lNNZ);
    *lcind = SafeMalloc<int>(lNNZ);
    *lrptr = SafeMalloc<int>(N+1);
    if(uNNZ != 0) {
        *uval  = SafeMalloc<double>(uNNZ);
        *ucind = SafeMalloc<int>(uNNZ);
        *urptr = SafeMalloc<int>(N+1);
    }
    if(dNNZ == N) { *diag = SafeMalloc<double>(dNNZ); }
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
                lv_ptr[*l_ptr] = val[j];
                lc_ptr[*l_ptr] = cind[j];
                (*l_ptr)++;
            }else if(cind[j] > i) {
                uv_ptr[*u_ptr] = val[j];
                uc_ptr[*u_ptr] = cind[j];
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

void Expand(
    double *lval, int *lcind, int *lrptr,
    double **val, int **cind, int **rptr,
    int N)
{
    double *uval;
    int *ucind;
    int *urptr;
    int lnnz = lrptr[N];
    int nnz = lnnz*2 - N;
    Csr2Csc(lval, lcind, lrptr, &uval, &ucind, &urptr, N, N);
    *val = new double[nnz];
    *cind = new int[nnz];
    *rptr = new int[N+1];
    int count = 0;
    (*rptr)[0] = count;
    for(int i=0; i<N; i++) {
        for(int j=lrptr[i]; j<lrptr[i+1]; j++) {
            (*val)[count] = lval[j];
            (*cind)[count] = lcind[j];
            count++;
        }
        for(int j=urptr[i]+1; j<urptr[i+1]; j++) {
            (*val)[count] = uval[j];
            (*cind)[count] = ucind[j];
            count++;
        }
        (*rptr)[i+1] = count;
    }
    delete[] uval;
    delete[] ucind;
    delete[] urptr;
}

void Duplicate(
    double *tval, int *tcind, int *trptr,
    double **val, int **cind, int **rptr,
    int N)
{
    *val = new double[trptr[N]];
    *cind = new int[trptr[N]];
    *rptr = new int[N+1];
    (*rptr)[0] = trptr[0];
    for(int i=0; i<N; i++) {
        (*rptr)[i+1] = trptr[i+1];
        for(int j=trptr[i]; j<trptr[i+1]; j++) {
            (*val)[j] = tval[j];
            (*cind)[j] = tcind[j];
        }
    }
}

void Scaling(double *val, int *cind, int *rptr, double *b, int N)
{
    for(int i=0; i<N; i++) {
        double max = std::abs(val[rptr[i]]);
        for(int j=rptr[i]+1; j<rptr[i+1]; j++) {
            if(std::abs(val[j]) > max) max = std::abs(val[j]);
        }
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            val[j] /= max;
        }
        b[i] /= max;
    }
}

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

int Csr2Sell(
    double *val, int *cind, int *rptr,
    double **s_val, int **s_cind, int **s_len,
    int size, int N)
{
    int num_slice = (N+size-1)/size;
    *s_len = new int[num_slice+1];
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
    *s_val = new double[nnz];
    *s_cind = new int[nnz];
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

void Csr2Csc(
    double *val, int *cind, int *rptr,
    double **cval, int **crind, int **ccptr,
    int N, int M)
{
    // The size of the input matrix is N * M
    int nnz = rptr[N];
    *cval = new double[nnz];
    *crind = new int[nnz];
    *ccptr = new int[M+1];

    int *num = (int*)calloc(M, sizeof(int));
    for(int i=0; i<N; i++) {
        for(int j=rptr[i]; j<rptr[i+1]; j++) {
            num[cind[j]]++;
        }
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

void Ilu0(double *val, int *cind, int *rptr, int N)
{
    for(int i=1; i<N; i++) {
        for(int k=rptr[i]; k<rptr[i+1]; k++) {
            if(cind[k] >= i) break;
            for(int l=rptr[cind[k]]; l<rptr[cind[k]+1]; l++) {
                if(cind[l] == cind[k]) {
                    if(val[l] == 0) {
                        printf("Error: Ilu0, 0 pivot\n");
                        std::exit(1);
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

//---- experimental ----//

static void append_dii(
    double **d, int **i1, int **i2, int *len, int *mlen,
    double d_val, int i1_val, int i2_val) {
    if(*len == *mlen) {
        *mlen *= 2;
        *d = utils::SafeRealloc<double>(*d, *mlen);
        *i1 = utils::SafeRealloc<int>(*i1, *mlen);
        *i2 = utils::SafeRealloc<int>(*i2, *mlen);
    }
    (*d)[*len] = d_val;
    (*i1)[*len] = i1_val;
    (*i2)[*len] = i2_val;
    *len += 1;
}

void Ilup(double **val, int **cind, int **rptr, int N, int p)
{
    int i, j;
    //int NNZ = (*rptr)[N];
    int LEN = 8;
    double *temp_val = utils::SafeMalloc<double>(LEN);
    int *temp_col    = utils::SafeMalloc<int>(LEN);
    int *temp_lev    = utils::SafeMalloc<int>(LEN);
    int temp_len     = 0;
    int temp_mlen    = LEN;
    double *temp2_val = utils::SafeMalloc<double>(LEN);
    int *temp2_col    = utils::SafeMalloc<int>(LEN);
    int *temp2_lev    = utils::SafeMalloc<int>(LEN);
    int temp2_len     = 0;
    int temp2_mlen    = LEN;
    double pivot = 0;
    int pivot_lev = 0;

    int first_len = (*rptr)[1] - (*rptr)[0];
    double *new_val = utils::SafeMalloc<double>(first_len);
    int *new_cind = utils::SafeMalloc<int>(first_len);
    int *new_lev = utils::SafeCalloc<int>(first_len);
    int new_len = first_len;
    int *new_rptr = SafeMalloc<int>(N+1);
    //一行目
    memcpy(new_val, *val, sizeof(double)*first_len);
    memcpy(new_cind, *cind, sizeof(int)*first_len);
    new_rptr[0] = 0;
    new_rptr[1] = first_len;

    for(i=1; i<N; i++) {
        int now_len = (*rptr)[i+1] - (*rptr)[i];
        while(temp_mlen < now_len) {
            temp_mlen *= 2;
            temp_val = utils::SafeRealloc<double>(temp_val, temp_mlen);
            temp_col = utils::SafeRealloc<int>(temp_col, temp_mlen);
            temp_lev = utils::SafeRealloc<int>(temp_lev, temp_mlen);
        }
        temp_len = now_len;
        memcpy(temp_val, &(*val)[(*rptr)[i]], sizeof(double)*now_len);
        memcpy(temp_col, &(*cind)[(*rptr)[i]], sizeof(int)*now_len);
        for(j=0; j<now_len; j++) temp_lev[j] = 0;
        int count = 0;
        while(count < temp_len && temp_col[count] < i) {
            int k_ptr = count;
            if(temp_lev[k_ptr] > p) {
                count++;
                continue;
            }
            int k = temp_col[k_ptr];
            temp2_len = 0;
            for(j=0; j<count; j++) {
                append_dii(&temp2_val, &temp2_col, &temp2_lev, &temp2_len, &temp2_mlen, temp_val[j], temp_col[j], temp_lev[j]);
            }
            for(j=new_rptr[k]; j<new_rptr[k+1]; j++) {
                if(new_cind[j] == k) {
                    pivot = temp_val[k_ptr] / new_val[j];
                    pivot_lev = temp_lev[k_ptr];
                    break;
                }
            }
            j++;
            k_ptr++;
            append_dii(&temp2_val, &temp2_col, &temp2_lev, &temp2_len, &temp2_mlen, pivot, k, 0);
            while(k_ptr<temp_len || j<new_rptr[k+1]) {
                int t_col1 = INT_MAX;
                int t_col2 = INT_MAX;
                if(k_ptr < temp_len) t_col1 = temp_col[k_ptr];
                if(j < new_rptr[k+1]) t_col2 = new_cind[j];
                if(t_col1<t_col2) {
                    double t_val = temp_val[k_ptr];
                    int t_col = temp_col[k_ptr];
                    int t_lev = temp_lev[k_ptr];
                    append_dii(&temp2_val, &temp2_col, &temp2_lev, &temp2_len, &temp2_mlen, t_val, t_col, t_lev);
                    k_ptr++;
                }else if(t_col1==t_col2) {
                    double t_val = temp_val[k_ptr] - pivot*new_val[j];
                    int t_col = temp_col[k_ptr];
                    int t_lev = temp_lev[k_ptr];
                    if(t_lev > pivot_lev+new_lev[j]+1) {
                        t_lev = pivot_lev+new_lev[j]+1;
                    }
                    append_dii(&temp2_val, &temp2_col, &temp2_lev, &temp2_len, &temp2_mlen, t_val, t_col, t_lev);
                    k_ptr++;
                    j++;
                }else { // (k>j) fill-in
                    double t_val = -pivot*new_val[j];
                    int t_col = new_cind[j];
                    int t_lev = pivot_lev+new_lev[j]+1;
                    //printf("lev %d\n", t_lev);
                    append_dii(&temp2_val, &temp2_col, &temp2_lev, &temp2_len, &temp2_mlen, t_val, t_col, t_lev);
                    j++;
                }
            }
            count++;
            temp_len = temp2_len;
            while(temp_mlen < temp2_len) {
                temp_mlen *= 2;
                temp_val = utils::SafeRealloc<double>(temp_val, temp_mlen);
                temp_col = utils::SafeRealloc<int>(temp_col, temp_mlen);
                temp_lev = utils::SafeRealloc<int>(temp_lev, temp_mlen);
            }
            memcpy(temp_val, temp2_val, sizeof(double)*temp2_len);
            memcpy(temp_col, temp2_col, sizeof(int)*temp2_len);
            memcpy(temp_lev, temp2_lev, sizeof(int)*temp2_len);
        }
        new_val  = utils::SafeRealloc<double>(new_val, new_len+temp_len);
        new_cind = utils::SafeRealloc<int>(new_cind, new_len+temp_len);
        new_lev  = utils::SafeRealloc<int>(new_lev, new_len+temp_len);
        for(j=0; j<temp_len; j++) {
            if(temp_lev[j] <= p) {
                new_val[new_len] = temp_val[j];
                new_cind[new_len] = temp_col[j];
                new_lev[new_len] = temp_lev[j];
                new_len++;
            }
        }
        new_rptr[i+1] = new_len;
    }
    (*val)  = utils::SafeRealloc<double>(new_val, new_len);
    (*cind) = utils::SafeRealloc<int>(new_cind, new_len);
    (*rptr) = utils::SafeRealloc<int>(new_rptr, N+1);
    free(temp_val);
    free(temp_col);
    free(temp_lev);
    free(temp2_val);
    free(temp2_col);
    free(temp2_lev); 
}

}

}
