#include <iostream>
#include <string>
#include <chrono>

#include "senk.hpp"
#include "senk_test.hpp"

int main(int argc, char **argv)
{
    std::string filename = "/Users/kengo/HPC/matrix/";
    filename += (argv[1]);

    double *tval;
    int *tcind;
    int *trptr;
    int N, M;
    senk::io::Shape shape;
    bool removeZeros = true;
    senk::io::ReadMatrixMarket(
        filename, &tval, &tcind, &trptr,
        &N, &M, &shape, removeZeros);
    //senk::matrix::RemoveZeros(&tval, &tcind, trptr, N);
    printf("# NNZ : %d\n", trptr[N]);
    double *val;
    int *cind;
    int *rptr;
    if(shape == senk::io::Shape::Sym) {
        senk::matrix::Expand(tval, tcind, trptr, &val, &cind, &rptr, N);
    }else {
        senk::matrix::Duplicate(tval, tcind, trptr, &val, &cind, &rptr, N);
    }
    free(tval);
    free(tcind);
    free(trptr);
    if(!senk::matrix::CheckStructure(val, cind, rptr, N)) exit(1);

    printf("%d\n", rptr[N]);
    int ori_N = N;
    const int bnl = 8;
    const int bnw = 1;
    senk::matrix::Padding(&val, &cind, &rptr, bnl, &N);
    printf("# NNZ : %d -> %d\n", ori_N, M);

    double *b = senk::utils::SafeMalloc<double>(N);
    for(int i=0; i<N; i++) {
        if(i < ori_N) b[i] = 1;
        else b[i] = 0;
    }
    double *x = senk::utils::SafeMalloc<double>(N);
    for(int i=0; i<N; i++) {
        x[i] = 0;
    }

    senk::matrix::Duplicate(val, cind, rptr, &tval, &tcind, &trptr, N);

    double *bval;
    int *bcind;
    int *brptr;
    senk::matrix::Csr2Bcsr<bnl, bnw>(tval, tcind, trptr, &bval, &bcind, &brptr, N);
    free(tval);
    free(tcind);
    free(trptr);
    senk::matrix::Bcsr2Csr<bnl, bnw>(bval, bcind, brptr, &tval, &tcind, &trptr, N);
    
    senk::matrix::Ilu0(tval, tcind, trptr, N);

    //senk::matrix::Ilup(&tval, &tcind, &trptr, N, 5);

    double *lval;
    int *lcind;
    int *lrptr;
    double *uval;
    int *ucind;
    int *urptr;
    senk::matrix::Split(
        tval, tcind, trptr,
        &lval, &lcind, &lrptr,
        &uval, &ucind, &urptr,
        nullptr, N, "L-DU", true);

    double *blval;
    int *blcind;
    int *blrptr;
    double *buval;
    int *bucind;
    int *burptr;

    senk::matrix::Csr2Bcsr<bnl, bnw>(lval, lcind, lrptr, &blval, &blcind, &blrptr, N);
    senk::matrix::Csr2Bcsr<bnl, bnw>(uval, ucind, urptr, &buval, &bucind, &burptr, N);

    int max_iter = 6000;
    double epsilon = 1.0e-8;
    double nrm_b = senk::blas1::Nrm2(b, N);  

    senk::solver::Gmresm<double>(
        val, cind, rptr,
        b, x, nrm_b, max_iter/50, 50, N, epsilon);
    //senk::solver::IluGmresm<double>(
    //    val, cind, rptr,
    //    lval, lcind, lrptr, uval, ucind, urptr,
    //    b, x, nrm_b, max_iter/50, 50, N, epsilon);
    //senk::solver::IlubGmresm<double, bnl, bnw>(
    //    val, cind, rptr,
    //    blval, blcind, blrptr, buval, bucind, burptr,
    //    b, x, nrm_b, max_iter/50, 50, N, epsilon);

    //senk::solver::Bicgstab<double>(
    //    val, cind, rptr,
    //    b, x, nrm_b, max_iter, N, epsilon);
    //senk::solver::IluBicgstab<double>(
    //    val, cind, rptr,
    //    lval, lcind, lrptr, uval, ucind, urptr,
    //    b, x, nrm_b, max_iter, N, epsilon);
    //senk::solver::IlubBicgstab<double, bnl, bnw>(
    //    val, cind, rptr,
    //    blval, blcind, blrptr, buval, bucind, burptr,
    //    b, x, nrm_b, max_iter, N, epsilon);
    
    //senk::solver::Gcrm<double>(
    //    val, cind, rptr,
    //    b, x, nrm_b,
    //    max_iter/50, 50, N, epsilon);
    //senk::solver::IluGcrm<double>(
    //    val, cind, rptr,
    //    lval, lcind, lrptr, uval, ucind, urptr,
    //    b, x, nrm_b,
    //    max_iter/50, 50, N, epsilon);
    //senk::solver::IlubGcrm<double, bnl, bnw>(
    //    val, cind, rptr,
    //    blval, blcind, blrptr, buval, bucind, burptr,
    //    b, x, nrm_b,
    //    max_iter/50, 50, N, epsilon);
    
    senk::test::RelativeResidualError(
        val, cind, rptr, b, x, N);
    
    senk::utils::SafeFree<double>(&x);

    return 0;
}
