/**
 * @file senk_graph.hpp
 * @brief Functions related to the graph theory are written.
 * @author Kengo Suzuki
 * @date 5/9/2022
 */
#ifndef SENK_GRAPH_HPP
#define SENK_GRAPH_HPP

#include "senk_utils.hpp"

namespace senk {
/**
 * @brief Contains functions related to the graph theory.
 */
namespace graph {
/**
 * @brief Partition a input graph into blocks by using simple blocking.
 * @param cind An array that stores column indices.
 * @param rptr An array that stores row pointer.
 * @param b_cind A pointer to an array that stores column indices of the resulting blocked matrix.
 * @param b_rptr A pointer to an array that stores row pointer of the resulting blocked matrix.
 * @param block A pointer to an array that maps the original index to the block.
 * @param N The size of a input graph (matrix)
 * @param bsize The size of the block.
 */
static void SimpleBlocking( // Simple Blocking
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
/**
 * @brief Partition a input graph into blocks by using connected blocking.
 * @param cind An array that stores column indices.
 * @param rptr An array that stores row pointer.
 * @param b_cind A pointer to an array that stores column indices of the resulting blocked matrix.
 * @param b_rptr A pointer to an array that stores row pointer of the resulting blocked matrix.
 * @param block A pointer to an array that maps the original index to the block.
 * @param N The size of a input graph (matrix)
 * @param bsize The size of the block.
 */
static void ConnectedBlocking( // Connected Blocking
    int *cind, int *rptr,
    int **b_cind, int **b_rptr, int **block,
    int N, int bsize)
{
    *block = utils::SafeCalloc<int>(N);
    int *inv_block = utils::SafeCalloc<int>(N);
    int *visited = utils::SafeCalloc<int>(N);
    int *queue = utils::SafeCalloc<int>(bsize+1);
    int qs, qe;
    int count = 0, tcnt = 0;
    for(int i=0; i<N; i++) {
        if(visited[i]) continue;
        qs = 0; qe = 0; // Clear
        (*block)[count] = i; count++; tcnt++;
        inv_block[i] = (count-1) / bsize;
        visited[i] = 1;
        if(tcnt == bsize) { tcnt = 0; continue; }
        queue[qe] = i; qe++; // Enqueue
        bool isFull = false;
        while(!isFull) {
            int q = queue[qs]; qs++; // Dequeue
            for(int j=rptr[q]; j<rptr[q+1]; j++) {
                if(visited[cind[j]]) continue;
                (*block)[count] = cind[j]; count++; tcnt++;
                inv_block[cind[j]] = (count-1) / bsize;
                visited[cind[j]] = 1;
                if(tcnt == bsize) { tcnt = 0; isFull=true; break; }
                queue[qe] = cind[j]; qe++; // Enqueue
            }
            if(qs > qe) { exit(1); }
            if(qs == qe) { break; }
        }
    }
    for(int i=0; i<N; i+=bsize) senk::helper::QuickSort<int>(*block, i, i+bsize-1);

    *b_rptr = utils::SafeMalloc<int>(N/bsize+1);
    *b_cind = utils::SafeMalloc<int>(rptr[N]); // At most rptr[N]
    int *temp = utils::SafeMalloc<int>(rptr[N]); // At most rptr[N]
    int cnt = 0;
    (*b_rptr)[0] = 0;
    for(int i=0; i<N; i+=bsize) {
        int len = 0;
        for(int j=0; j<bsize; j++) {
            for(int k=rptr[(*block)[i+j]]; k<rptr[(*block)[i+j]+1]; k++) {
                temp[len] = inv_block[cind[k]]; len++;
            }
        }
        helper::QuickSort<int>(temp, 0, len-1);
        (*b_cind)[cnt] = temp[0]; cnt++;
        int prev = temp[0];
        for(int j=1; j<len; j++) {
            if(temp[j] == prev) continue;
            (*b_cind)[cnt] = temp[j]; cnt++;
            prev = temp[j];
        }
        (*b_rptr)[i/bsize+1] = cnt;
    }
    free(visited);
    free(queue);
    free(inv_block);
    free(temp);
}

/**
 * @brief Create an adjacency matrix from a input matrix.
 * @param cind An array that stores column indices.
 * @param rptr An array that stores row pointer.
 * @param am_cind A pointer to an array that stores column indices of the resulting adjacency matrix.
 * @param am_rptr A pointer to an array that stores row pointer of the resulting adjacency matrix.
 * @param N The size of a input graph (matrix)
 * @param isSym It true, a procedure for symmetric matrices is performed.
 */
static void GetAdjacency(
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
/**
 * @brief Coloring an adjacency matrix from a input matrix.
 * @param cind An array that stores column indices.
 * @param rptr An array that stores row pointer.
 * @param color A pointer to an array that maps indices to colors.
 * @param num_color An array that stores the number of indices in each color. 
 * @param N The size of a input graph (matrix)
 */
static void Coloring( // Greedy Coloring
    int *cind, int *rptr, int **color, int *num_color, int N)
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
            for(int j=rptr[i]; j<rptr[i+1]; j++) {
                visit[cind[j]] = now;
            }
        }
        now++;
    }
    *num_color = now - 1;
    free(visit);
}
/**
 * @brief Create an permutation matrix based on the AMC ordering technique@cite iwashita2002AMC.
 * @param cind An array that stores column indices.
 * @param rptr An array that stores row pointer.
 * @param num_color A variable to receive the number of colors.
 * @param size_color A pointer to receive the size of each color.
 * @param LP A pointer to the resulting vector representing the left permutation matrix.
 * @param RP A pointer to the resulting vector representing the right permutation matrix.
 * @param N The size of the input graph (matrix)
 * @param isSym Whether the input matrix is symmetric or not.
 */
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
/**
 * @brief Create an permutation matrix based on the ABMC ordering technique@cite iwashita2012ABMC, @cite li2019enhancement.
 * @param cind An array that stores column indices.
 * @param rptr An array that stores row pointer.
 * @param num_color A variable to receive the number of colors.
 * @param size_color A pointer to receive the size of each color.
 * @param LP A pointer to the resulting vector representing the left permutation matrix.
 * @param RP A pointer to the resulting vector representing the right permutation matrix.
 * @param N The size of the input graph (matrix)
 * @param bsize The size of the block.
 * @param isSym Whether the input matrix is symmetric or not.
 * @param bmethod Specify the blocking method: "simple" or "connect".
 */
void GetABMCPermutation(
    int *cind, int *rptr,
    int *num_color, int **size_color, int **LP, int **RP,
    int N, int bsize, bool isSym, const char *bmethod)
{
    if(N%bsize) {
        printf("Error: GetABMCPermutation\n");
        exit(1);
    }
    int *b_cind;
    int *b_rptr;
    int *block;
    if(std::strcmp(bmethod, "simple") == 0) {
        SimpleBlocking(cind, rptr, &b_cind, &b_rptr, &block, N, bsize);
    }else if(std::strcmp(bmethod, "connect") == 0) {
        ConnectedBlocking(cind, rptr, &b_cind, &b_rptr, &block, N, bsize);
    }else {
        printf("Error: GetABMCPermutation\n");
        printf("The blocking method is invalid\n");
        exit(1);
    }
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
