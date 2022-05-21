#ifndef SENK_TEST_HPP
#define SENK_TEST_HPP

#include "senk_blas1.hpp"
#include "senk_sparse.hpp"

namespace senk {

namespace test {

void RelativeResidualError(
    double *val, int *cind, int *rptr,
    double *b, double *x, int N)
{
    double *t = (double*)malloc(sizeof(double)*N);
    senk::sparse::SpmvCsr<double>(val, cind, rptr, x, t, N);
    senk::blas1::Axpby<double>(1, b, -1, t, N);
    double nrm_t = senk::blas1::Nrm2<double>(t, N);
    double nrm_b = senk::blas1::Nrm2<double>(b, N);
    printf("# test %e\n", nrm_t/nrm_b);
    free(t);
}

}

}

#endif

