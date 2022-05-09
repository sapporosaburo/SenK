/**
 * @file senk_utils.hpp
 * @brief Utility functions are defined.
 * @author Kengo Suzuki
 * @date 5/9/2022
 */
#ifndef SENK_UTILS_HPP
#define SENK_UTILS_HPP

#include <cstdlib>

namespace senk {
/**
 * @brief Contains utility functions.
 */
namespace utils {
/**
 * @brief Allocate memory
 * @tparam An unit type
 * @param size The size of memory to be allocated.
 */
template <typename T>
T *SafeMalloc(int size)
{
    T *res = (T*)std::malloc(sizeof(T)*size);
    if(!res) { printf("Error: SafeMalloc\n"); exit(1); }
    else { return res; }
}
/**
 * @brief Allocate memory and clear it.
 * @tparam An unit type
 * @param size The size of memory to be allocated.
 */
template <typename T>
T *SafeCalloc(int size)
{
    T *res = (T*)std::calloc(size, sizeof(T));
    if(!res) { printf("Error: SafeCalloc\n"); exit(1); }
    else { return res; }
}
/**
 * @brief Reallocate memory.
 * @tparam An unit type
 * @param old The pointer to the memory to be reallocated.
 * @param size The size of memory to be allocated.
 */
template <typename T>
T *SafeRealloc(T *old, int size)
{
    T *res = (T*)std::realloc(old, sizeof(T)*size);
    if(!res) { printf("Error: SafeRealloc\n"); exit(1); }
    else { return res; }
}
/**
 * @brief Free allocated memory.
 * @tparam T An unit type
 * @param ptr The pointer to the memory to be freed.
 */
template <typename T>
void SafeFree(T **ptr)
{
    if(*ptr) {free(*ptr);}
    *ptr = nullptr;
}
/**
 * @brief Copy in to out.
 * @tparam T The type of in and out.
 * @param in The input array.
 * @param out The output array. 
 * @param size The size of the arrays.
 */
template <typename T>
void Copy(T *in, T *out, int size)
{
    #pragma omp parallel for
    for(int i=0; i<size; i++) {
        out[i] = in[i];
    }
}
/**
 * @brief Set values of out to val.
 * @tparam T The type of the value.
 * @param val The input value.
 * @param out The output array. 
 * @param size The size of the array.
 */
template <typename T>
void Set(T val, T *out, int size)
{
    #pragma omp parallel for
    for(int i=0; i<size; i++) {
        out[i] = val;
    }
}
/**
 * @brief Convert type of arrays.
 * @tparam T1 The type of the original array.
 * @tparam T2 The type of the converted array.
 * @param in The input array.
 * @param out The output array.
 * @param size The size of the array.
 */
template <typename T1, typename T2>
void Convert(T1 *in, T2 *out, int size)
{
    #pragma omp parallel for
    for(int i=0; i<size; i++) {
        out[i] = (T2)in[i];
    }
}
/*
template <typename T1, typename T2, int bit>
void Convert(T1 *in, T2 *out, int size)
{
    #pragma omp parallel for
    for(int i=0; i<size; i++) {
        out[i] = (T2)(in[i] * (1<<bit));
    }
}
*/
}

}

#endif
