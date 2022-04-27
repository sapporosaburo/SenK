#ifndef SENK_TEST_HPP
#define SENK_TEST_HPP

namespace senk {

namespace test {

void RelativeResidualError(
    double *val, int *cind, int *rptr,
    double *b, double *x, int N);

}

}

#endif

