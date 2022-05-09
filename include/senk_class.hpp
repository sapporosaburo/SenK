/**
 * @file senk_class.hpp
 * @brief Classes used in SenK are written.
 * @author Kengo Suzuki
 * @date 5/9/2022
 */
#ifndef SENK_CLASS_HPP
#define SENK_CLASS_HPP

#include <tuple>
#include "senk_utils.hpp"

#define SP_LEN 8

namespace senk {

namespace sparse {
/**
 * @brief Sparse vector class
 * @tparam T1 Type of the val array.
 * @tparam T2 Type of the lev array.
 */
template <typename T1, typename T2>
class SpVec {
private:
    //! An array that stores the values of the nonzero elements.
    T1 *val;
    //! An array that stores the supplemental values of the nonzero elements.
    T2 *lev;
    //! An array that stores the indices of the nonzero elements.
    int *idx;
    //! The number of the nonzero elements.
    int len;
    //! The size of allocated memories.
    int mlen;
public:
    /**
     * @brief Constructor.
     * @details Allocate memory of size SP_LEN.
     */
    SpVec() {
        val = utils::SafeMalloc<T1>(SP_LEN);
        lev = utils::SafeMalloc<T2>(SP_LEN);
        idx = utils::SafeMalloc<int>(SP_LEN);
        len = 0;
        mlen = SP_LEN;
    }
    /**
     * @brief Destructor.
     * @details Free all memories for the arrays. 
     */
    ~SpVec() {
        utils::SafeFree(&val);
        utils::SafeFree(&lev);
        utils::SafeFree(&idx);
    }
    /**
     * @brief Append argument values to the arrays.
     * @param t_val A value to be appended to val.
     * @param t_lev A value to be appended to lev.
     * @param t_idx A value to be appended to idx.
     */
    inline void Append(T1 t_val, T2 t_lev, int t_idx) {
        if(len == mlen) {
            mlen *= 2;
            val = utils::SafeRealloc(val, mlen);
            lev = utils::SafeRealloc(lev, mlen);
            idx = utils::SafeRealloc(idx, mlen);
        }
        val[len] = t_val;
        lev[len] = t_lev;
        idx[len] = t_idx;
        len++;
    }
    /**
     * @brief Append argument values to the arrays.
     * @param t_val Values to be appended to val.
     * @param t_lev Values to be appended to lev.
     * @param t_idx Values to be appended to idx.
     * @param t_len Length of these argument values.
     */
    inline void Append(T1 *t_val, T2 *t_lev, int *t_idx, int t_len) {
        if(len+t_len > mlen) {
            while(len+t_len > mlen) mlen *= 2;
            val = utils::SafeRealloc(val, mlen);
            lev = utils::SafeRealloc(lev, mlen);
            idx = utils::SafeRealloc(idx, mlen);
        }
        for(int i=0; i<t_len; i++) {
            val[len+i] = t_val[i];
            lev[len+i] = t_lev[i];
            idx[len+i] = t_idx[i];
        }
        len += t_len;
    }
    /**
     * @brief Return len.
     */
    inline int GetLen() { return len; }
    /**
     * @brief Return i-th value of val.
     */
    inline T1 GetVal(int i) { return val[i]; }
    /**
     * @brief Return i-th value of lev.
     */
    inline T2 GetLev(int i) { return lev[i]; }
    /**
     * @brief Return i-th value of idx.
     */
    inline int GetIdx(int i) { return idx[i]; }
    /**
     * @brief Return the tuple of i-th values of val, lev, and idx.
     */
    inline std::tuple<T1, T2, int> Get(int i) {
        //return std::make_tuple(val[i], idx[i]);
        return {val[i], lev[i], idx[i]};
    }
    /**
     * @brief Set the value of len to 0.
     */
    inline void Clear() { len = 0; }
};

/*
template <typename T>
class SpVec {
private:
    T *val;
    int *idx;
    int len;
    int mlen;
public:
    SpVec() {
        val = utils::SafeMalloc<T>(SP_LEN);
        idx = utils::SafeMalloc<int>(SP_LEN);
        len = 0;
        mlen = SP_LEN;
    }
    ~SpVec() {
        utils::SafeFree(&val);
        utils::SafeFree(&idx);
    }
    inline int GetLen() { return len; }
    inline void Append(T t_val, int t_idx) {
        if(len == mlen) {
            mlen *= 2;
            val = utils::SafeRealloc(val, mlen);
            idx = utils::SafeRealloc(idx, mlen);
        }
        val[len] = t_val;
        idx[len] = t_idx;
        len++;
    }
    inline void Append(T *t_val, int *t_idx, int t_len) {
        if(len+t_len > mlen) {
            while(len+t_len > mlen) mlen *= 2;
            val = utils::SafeRealloc(val, mlen);
            idx = utils::SafeRealloc(idx, mlen);
        }
        for(int i=0; i<t_len; i++) {
            val[len+i] = t_val[i];
            idx[len+i] = t_idx[i];
        }
        len += t_len;
    }
    inline std::tuple<T, int> Get(int i) {
        //return std::make_tuple(val[i], idx[i]);
        return {val[i], idx[i]};
    }
    inline void Clear() { len = 0; }
};
*/

} // namespace sparse

} // namespace senk

#endif
