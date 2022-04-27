#ifndef SENK_IO_HPP
#define SENK_IO_HPP

namespace senk {

namespace io {

enum Shape {
    Sym,
    Unsym
};

bool ReadMatrixMarket(
    std::string filename, double **val, int **cind, int **rptr,
    int *N, int *M, Shape *shape, bool removeZeros);

}

}

#endif
