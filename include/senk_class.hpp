#ifndef SENK_CLASS_HPP
#define SENK_CLASS_HPP

#include <iostream>
#include <tuple>
#include "senk_utils.hpp"

#define SP_LEN 8

namespace senk {

namespace sparse {
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
template <typename T1, typename T2>
class SpVec {
private:
    T1 *val;
    T2 *lev;
    int *idx;
    int len;
    int mlen;
public:
    SpVec() {
        val = utils::SafeMalloc<T1>(SP_LEN);
        lev = utils::SafeMalloc<T2>(SP_LEN);
        idx = utils::SafeMalloc<int>(SP_LEN);
        len = 0;
        mlen = SP_LEN;
    }
    ~SpVec() {
        utils::SafeFree(&val);
        utils::SafeFree(&lev);
        utils::SafeFree(&idx);
    }
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
    inline int GetLen() { return len; }
    inline T1 GetVal(int i) { return val[i]; }
    inline T2 GetLev(int i) { return lev[i]; }
    inline int GetIdx(int i) { return idx[i]; }
    inline std::tuple<T1, T2, int> Get(int i) {
        //return std::make_tuple(val[i], idx[i]);
        return {val[i], lev[i], idx[i]};
    }
    inline void Clear() { len = 0; }
};

} // namespace sparse

} // namespace senk

#endif
