#ifndef SENK_HELPER_HPP
#define SENK_HELPER_HPP

namespace senk {

namespace helper {

template <typename T>
void Swap(T *a, T *b)
{
    T temp = a[0];
    a[0] = b[0];
    b[0] = temp;
}

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

} // namespace helper

} // namespace senk

#endif
