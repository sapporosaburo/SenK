#ifndef SENK_UTILS_HPP
#define SENK_UTILS_HPP

#include <cstdlib>
//#include <omp.h>

namespace senk {

namespace utils {

template <typename T>
T *SafeMalloc(int size)
{
    T *res = (T*)std::malloc(sizeof(T)*size);
    if(!res) { printf("Error: SafeMalloc\n"); exit(1); }
    else { return res; }
}

template <typename T>
T *SafeCalloc(int size)
{
    T *res = (T*)std::calloc(size, sizeof(T));
    if(!res) { printf("Error: SafeCalloc\n"); exit(1); }
    else { return res; }
}

template <typename T>
T *SafeRealloc(T *old, int size)
{
    T *res = (T*)std::realloc(old, sizeof(T)*size);
    if(!res) { printf("Error: SafeRealloc\n"); exit(1); }
    else { return res; }
}

template <typename T>
void SafeFree(T **ptr)
{
    if(*ptr) {free(*ptr);}
    *ptr = nullptr;
}

template <typename T>
void Copy(T *in, T *out, int size)
{
    #pragma omp parallel for
    for(int i=0; i<size; i++) {
        out[i] = in[i];
    }
}

template <typename T>
void Set(T val, T *out, int size)
{
    #pragma omp parallel for
    for(int i=0; i<size; i++) {
        out[i] = val;
    }
}

template <typename T1, typename T2>
void Convert(T1 *in, T2 *out, int size)
{
    #pragma omp parallel for
    for(int i=0; i<size; i++) {
        out[i] = (T2)in[i];
    }
}

template <typename T1, typename T2, int bit>
void Convert(T1 *in, T2 *out, int size)
{
    #pragma omp parallel for
    for(int i=0; i<size; i++) {
        out[i] = (T2)(in[i] * (1<<bit));
    }
}

}

}

#endif
