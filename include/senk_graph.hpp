#ifndef SENK_GRAPH_HPP
#define SENK_GRAPH_HPP

#include "senk_utils.hpp"

namespace senk {

namespace graph {

void Blocking( // Simple Blocking
    int *cind, int *rptr,
    int **b_cind, int **b_rptr, int **block,
    int N, int bsize)
{
    *block = utils::SafeMalloc<int>(N);
    for(int i=0; i<N; i++) { (*block)[i] = i; }
    
    *b_rptr = utils::SafeMalloc<int>(N/bsize+1);
    *b_cind = utils::SafeMalloc<int>(rptr[N]); // At most rptr[N]
    int cnt = 0;
    (*b_rptr)[0] = 0;
    int *ptr = utils::SafeMalloc<int>(bsize);
    for(int i=0; i<N; i+=bsize) {
        // Initialize "ptr"
        for(int j=0; j<bsize; j++) { ptr[j] = rptr[(*block)[i+j]]; }
        while(true) {
            // Find minimum col index
            int min = N;
            for(int j=0; j<bsize; j++) {
                if(ptr[j] != -1 && cind[ptr[j]] < min) {
                    min = cind[ptr[j]];
                }
            }
            if(min == N) break;
            (*b_cind)[cnt] = min/bsize;
            for(int j=0; j<bsize; j++) {
                if(ptr[j] == -1) continue;
                while(cind[ptr[j]]/bsize == min/bsize) {
                    ptr[j]++;
                    if(ptr[j] >= rptr[(*block)[i+j]+1]) {
                        ptr[j] = -1; break;
                    }
                }
            }
            cnt++;
        }
        (*b_rptr)[i/bsize+1] = cnt;
    }
}

void GetAdjacency(
    int *cind, int *rptr,
    int **am_cind, int **am_rptr,
    int N, bool isSym)
{
    if(isSym) {
        *am_cind = utils::SafeMalloc<int>(rptr[N]);
        *am_rptr = utils::SafeMalloc<int>(N+1);
        utils::Copy<int>(cind, *am_cind, rptr[N]);
        utils::Copy<int>(rptr, *am_rptr, N+1);
    }else {
        int *rind = utils::SafeMalloc<int>(rptr[N]);
        int *cptr = utils::SafeMalloc<int>(N+1);
        int *len = utils::SafeCalloc<int>(N);
        for(int i=0; i<N; i++) {
            for(int j=rptr[i]; j<rptr[i+1]; j++) {
                len[cind[j]]++;
            }
        }
        cptr[0] = 0;
        for(int i=0; i<N; i++) {
            cptr[i+1] = cptr[i] + len[i];
            len[i] = 0;
        }
        for(int i=0; i<N; i++) {
            for(int j=rptr[i]; j<rptr[i+1]; j++) {
                rind[cptr[cind[j]]+len[cind[j]]] = i;
                len[cind[j]]++;
            }
        }
        *am_cind = utils::SafeMalloc<int>(rptr[N]*2);
        *am_rptr = utils::SafeMalloc<int>(N+1);
        (*am_rptr)[0] = 0;
        int cnt = 0;
        for(int i=0; i<N; i++) {
            int csr_ptr = rptr[i];
            int csc_ptr = cptr[i];
            while(csr_ptr < rptr[i+1] || csc_ptr < cptr[i+1]) {
                int csr_idx = (csr_ptr < rptr[i+1])? cind[csr_ptr] : N;
                int csc_idx = (csc_ptr < cptr[i+1])? rind[csc_ptr] : N;
                if(csr_idx < csc_idx) {
                    (*am_cind)[cnt] = csr_idx;
                    csr_ptr++;
                }else if(csr_idx == csc_idx){
                    (*am_cind)[cnt] = csr_idx;
                    csr_ptr++;
                    csc_ptr++;
                }else {
                    (*am_cind)[cnt] = csc_idx;
                    csc_ptr++;
                }
                cnt++;
            }
            (*am_rptr)[i+1] = cnt;
        }
        free(rind);
        free(cptr);
        free(len);
    }
}

void Coloring( // Greedy Coloring
    int *am_cind, int *am_rptr, int **color, int *num_color, int N)
{
    *color = utils::SafeCalloc<int>(N);
    int *visit = utils::SafeCalloc<int>(N);
    int now = 1;
    int cnt = 0;
    while(cnt < N) {
        for(int i=0; i<N; i++) {
            if(visit[i] == now || (*color)[i] != 0) continue;
            (*color)[i] = now;
            cnt++;
            for(int j=am_rptr[i]; j<am_rptr[i+1]; j++) {
                visit[am_cind[j]] = now;
            }
        }
        now++;
    }
    *num_color = now - 1;
    free(visit);
}

void GetAMCPermutation(
    int *cind, int *rptr,
    int *num_color, int **size_color, int **LP, int **RP,
    int N, bool isSym)
{
    int *am_cind;
    int *am_rptr;
    GetAdjacency(cind, rptr, &am_cind, &am_rptr, N, isSym);
    int *color;
    // Greedy Coloring
    Coloring(am_cind, am_rptr, &color, num_color, N);
    *LP = utils::SafeMalloc<int>(N);
    *size_color = utils::SafeMalloc<int>(*num_color+1);
    (*size_color)[0] = 0;
    int cnt = 0;
    for(int i=0; i<*num_color; i++) {
        int temp_cnt = 0;
        for(int j=0; j<N; j++) {
            if(color[j]-1 != i) continue;
            (*LP)[cnt] = j;
            cnt++;
            temp_cnt++;
        }
        (*size_color)[i+1] = (*size_color)[i] + temp_cnt;
    }
    *RP = utils::SafeMalloc<int>(N);
    for(int i=0; i<N; i++) {(*RP)[(*LP)[i]] = i;}
}

void GetABMCPermutation(
    int *cind, int *rptr,
    int *num_color, int **size_color, int **LP, int **RP,
    int N, int bsize, bool isSym)
{
    if(N%bsize) {
        printf("Error: GetABMCPermutation\n");
        exit(1);
    }
    int *b_cind;
    int *b_rptr;
    int *block;
    // Simple Blocking (or Grouping)
    Blocking(cind, rptr, &b_cind, &b_rptr, &block, N, bsize);
    int *am_cind;
    int *am_rptr;
    GetAdjacency(b_cind, b_rptr, &am_cind, &am_rptr, N/bsize, isSym);
    int *color;
    // Greedy Coloring
    Coloring(am_cind, am_rptr, &color, num_color, N/bsize);
    *LP = utils::SafeMalloc<int>(N);
    *size_color = utils::SafeMalloc<int>(*num_color+1);
    (*size_color)[0] = 0;
    int cnt = 0;
    for(int i=0; i<*num_color; i++) {
        for(int j=0; j<N/bsize; j++) {
            if(color[j]-1 != i) continue;
            for(int k=0; k<bsize; k++) {
                (*LP)[cnt] = block[j*bsize+k];
                cnt++;
            }
        }
        (*size_color)[i+1] = cnt/bsize;
    }
    *RP = utils::SafeMalloc<int>(N);
    for(int i=0; i<N; i++) {(*RP)[(*LP)[i]] = i;}
}


} // namespace graph

} // namespace senk

#endif
