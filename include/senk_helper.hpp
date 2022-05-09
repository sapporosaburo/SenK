/**
 * @file senk_helper.hpp
 * @brief Some supplemental functions are defined.
 * @author Kengo Suzuki
 * @date 5/9/2022
 */
#ifndef SENK_HELPER_HPP
#define SENK_HELPER_HPP

namespace senk {
/**
 * @brief Contains supplemental functions.
 */
namespace helper {
/**
 * @brief Swap two variables.
 * @tparam T Type of the variables.
 * @param a The variable.
 * @param b The variable.
 */
template <typename T>
void Swap(T *a, T *b)
{
    T temp = a[0];
    a[0] = b[0];
    b[0] = temp;
}
/**
 * @brief Sort an array by the quick sort.
 * @tparam T Type of the array.
 * @param key The variable.
 * @param left Starting position of the range to be sorted.
 * @param right End position of the range to be sorted.
 */
template <typename T>
void QuickSort(T *key, int left, int right)
{
    int Left, Right;
    T pivot;
    Left = left; Right = right;
    pivot = key[(left + right) / 2];
    while (1) {
        while (key[Left] < pivot) Left++;
        while (pivot < key[Right]) Right--; 
        if (Left >= Right) break;
        Swap<T>(&key[Left], &key[Right]);
        Left++; Right--;
    }
    if (left < Left-1) QuickSort<T>(key, left, Left-1);
    if (Right+1 < right) QuickSort<T>(key, Right+1, right);
}
/**
 * @brief Sort arrays by the quick sort.
 * @tparam T Type of the key array.
 * @tparam T Type of the dependent array.
 * @param key A key array.
 * @param sub A subordinate array.
 * @param left Starting position of the range to be sorted.
 * @param right End position of the range to be sorted.
 */
template <typename T, typename T2>
void QuickSort(T *key, T2 *sub, int left, int right)
{
    int Left, Right;
    T pivot;
    Left = left; Right = right;
    pivot = key[(left + right) / 2];
    while (1) {
        while (key[Left] < pivot) Left++;
        while (pivot < key[Right]) Right--; 
        if (Left >= Right) break;
        Swap<T>(&key[Left], &key[Right]);
        Swap<T2>(&sub[Left], &sub[Right]);
        Left++; Right--;
    }
    if (left < Left-1) QuickSort<T, T2>(key, sub, left, Left-1);
    if (Right+1 < right) QuickSort<T, T2>(key, sub, Right+1, right);
}
/**
 * @brief Sort arrays by the quick sort.
 * @tparam T Type of the key array.
 * @tparam T2 Type of the 1st dependent array.
 * @tparam T3 Type of the 2nd dependent array.
 * @param key A key array.
 * @param sub A dependent array.
 * @param sub2 A dependent array.
 * @param left Starting position of the range to be sorted.
 * @param right End position of the range to be sorted.
 */
template <typename T, typename T2, typename T3>
void QuickSort(T *key, T2 *sub, T3 *sub2, int left, int right)
{
    int Left, Right;
    T pivot;
    Left = left; Right = right;
    pivot = key[(left + right) / 2];
    while (1) {
        while (key[Left] < pivot) Left++;
        while (pivot < key[Right]) Right--; 
        if (Left >= Right) break;
        Swap<T>(&key[Left], &key[Right]);
        Swap<T2>(&sub[Left], &sub[Right]);
        Swap<T3>(&sub2[Left], &sub2[Right]);
        Left++; Right--;
    }
    if (left < Left-1) QuickSort<T, T2, T3>(key, sub, sub2, left, Left-1);
    if (Right+1 < right) QuickSort<T, T2, T3>(key, sub, sub2, Right+1, right);
}
/**
 * @brief Sort an array in descending order by the quick sort.
 * @tparam T Type of the array.
 * @param key The variable.
 * @param left Starting position of the range to be sorted.
 * @param right End position of the range to be sorted.
 */
template <typename T>
void QuickSortDesc(T *key, int left, int right)
{
    int Left, Right;
    T pivot;
    Left = left; Right = right;
    pivot = key[(left + right) / 2];
    while (1) {
        while (key[Left] > pivot) Left++;
        while (pivot > key[Right]) Right--; 
        if (Left >= Right) break;
        Swap<T>(&key[Left], &key[Right]);
        Left++; Right--;
    }
    if (left < Left-1) QuickSortDesc<T>(key, left, Left-1);
    if (Right+1 < right) QuickSortDesc<T>(key, Right+1, right);
}
/**
 * @brief Sort arrays in descending order by the quick sort.
 * @tparam T Type of the key array.
 * @tparam T Type of the dependent array.
 * @param key A key array.
 * @param sub A subordinate array.
 * @param left Starting position of the range to be sorted.
 * @param right End position of the range to be sorted.
 */
template <typename T, typename T2>
void QuickSortDesc(T *key, T2 *sub, int left, int right)
{
    int Left, Right;
    T pivot;
    Left = left; Right = right;
    pivot = key[(left + right) / 2];
    while (1) {
        while (key[Left] > pivot) Left++;
        while (pivot > key[Right]) Right--; 
        if (Left >= Right) break;
        Swap<T>(&key[Left], &key[Right]);
        Swap<T2>(&sub[Left], &sub[Right]);
        Left++; Right--;
    }
    if (left < Left-1) QuickSortDesc<T, T2>(key, sub, left, Left-1);
    if (Right+1 < right) QuickSortDesc<T, T2>(key, sub, Right+1, right);
}
/**
 * @brief Sort arrays in descending order by the quick sort.
 * @tparam T Type of the key array.
 * @tparam T2 Type of the 1st dependent array.
 * @tparam T3 Type of the 2nd dependent array.
 * @param key A key array.
 * @param sub A dependent array.
 * @param sub2 A dependent array.
 * @param left Starting position of the range to be sorted.
 * @param right End position of the range to be sorted.
 */
template <typename T, typename T2, typename T3>
void QuickSortDesc(T *key, T2 *sub, T3 *sub2, int left, int right)
{
    int Left, Right;
    T pivot;
    Left = left; Right = right;
    pivot = key[(left + right) / 2];
    while (1) {
        while (key[Left] > pivot) Left++;
        while (pivot > key[Right]) Right--; 
        if (Left >= Right) break;
        Swap<T>(&key[Left], &key[Right]);
        Swap<T2>(&sub[Left], &sub[Right]);
        Swap<T3>(&sub2[Left], &sub2[Right]);
        Left++; Right--;
    }
    if (left < Left-1) QuickSortDesc<T, T2, T3>(key, sub, sub2, left, Left-1);
    if (Right+1 < right) QuickSortDesc<T, T2, T3>(key, sub, sub2, Right+1, right);
}
/*
template <typename T>
T Sqrt(T x)
{
    T s, t;
    if(x <= 0) return 0;
    s = 1; t = x;
    while(s < t) { s<<=1; t>>=1; }
    do { t = s; s=(x/s+s)>>1; } while (s < t);
    return t;
}
*/
} // namespace helper

} // namespace senk

#endif
