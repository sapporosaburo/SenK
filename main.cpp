#include <iostream>
#include <string>
#include <chrono>
#include <random>

#include "senk.hpp"
#include "senk_test.hpp"

int main(int argc, char **argv)
{
    std::string filename = "/path_to/matrix/";
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
    printf("# %d\n", trptr[N]);
    if(!senk::matrix::CheckStructure<double>(tval, tcind, trptr, N)) {
        printf("False: CheckStructure\n"); exit(1);
    }
    double *val;
    int *cind;
    int *rptr;
    if(shape == senk::io::Shape::Sym) {
        senk::matrix::Expand<double>(
            tval, tcind, trptr, &val, &cind, &rptr, N, "L");
    }else {
        senk::matrix::Duplicate<double>(
            tval, tcind, trptr, &val, &cind, &rptr, N);
    }
    free(tval);
    free(tcind);
    free(trptr);

// Padding
    int ori_N = N;
    senk::matrix::Padding(&val, &cind, &rptr, 128, &N);
    printf("%d\n", rptr[N]);
    M = N;
    printf("%d %d\n", N, M);

// Reordering
    int num_color;
    int *size_color;
    int *LP, *RP;
    //senk::graph::GetAMCPermutation(
    //    cind, rptr, &num_color, &size_color, &LP, &RP, N, false);
    senk::graph::GetABMCPermutation(
        cind, rptr, &num_color, &size_color, &LP, &RP, N, 128, false, "connect");
    senk::matrix::Reordering<double>(val, cind, rptr, LP, RP, N);
    printf("# ");
    for(int i=0; i<num_color; i++) {
        printf("%d ", size_color[i+1]-size_color[i]);   
    }
    printf("\n");
    printf("# Reordered!\n");

// Scaling
    //double *tdiag;
    //senk::matrix::GetDiag<double>(val, cind, rptr, &tdiag, N);
    //double *ones = senk::utils::SafeMalloc<double>(N);
    //for(int i=0; i<N; i++) { 
    //    tdiag[i] = 1 / std::sqrt(std::abs(tdiag[i]));
    //}
    //senk::utils::Set<double>(1.0, ones, N);
    //senk::matrix::Scaling<double>(val, cind, rptr, tdiag, tdiag, N);

// ILU factorization
    senk::matrix::Duplicate(val, cind, rptr, &tval, &tcind, &trptr, N);
    const int bnl = 8;
    const int bnw = 1;
    //senk::matrix::AllocBlockZero<double>(&tval, &tcind, &trptr, N, 2, 2);
    //senk::matrix::AllocLevelZero<double>(&tval, &tcind, &trptr, N, 1);    
    senk::matrix::AllocBlockZero<double>(&tval, &tcind, &trptr, N, bnl, bnw);
    printf("# Allocated!\n");

    //senk::matrix::AllocLevelZero(&tval, &tcind, &trptr, N, 2);
    senk::matrix::Ilu0<double>(tval, tcind, trptr, N);
    //senk::matrix::Ilup(&tval, &tcind, &trptr, N, 2);
    printf("# Factrized!\n");
    
    double *lval;
    int *lcind;
    int *lrptr;
    double *uval;
    int *ucind;
    int *urptr;
    bool diagInv = true;
    senk::matrix::Split<double>(
        tval, tcind, trptr,
        &lval, &lcind, &lrptr,
        &uval, &ucind, &urptr,
        nullptr, N, "L-DU", diagInv);
    printf("# Splited!\n");

    double *blval;
    int *blcind;
    int *blrptr;
    double *buval;
    int *bucind;
    int *burptr;

    senk::matrix::Csr2Bcsr<double>(lval, lcind, lrptr, &blval, &blcind, &blrptr, N, bnl, bnw);
    senk::matrix::Csr2Bcsr<double>(uval, ucind, urptr, &buval, &bucind, &burptr, N, bnl, bnw);

    double *x = senk::utils::SafeMalloc<double>(N);
    double *b = senk::utils::SafeMalloc<double>(N);
    senk::utils::Set<double>(0.0, x, N);
    //senk::utils::Set<double>(1.0, b, N);
    std::mt19937 engine(1);
    std::uniform_real_distribution<> dist1(-1.0, 1.0);
    for(int i=0; i<N; i++) { b[i] = dist1(engine); }
    
    for(int i=0; i<N-ori_N; i++) { b[ori_N+i] = 0.0; }

    int max_iter = 6000;
    double epsilon = 1.0e-8;
    double nrm_b = senk::blas1::Nrm2(b, N);  

    auto start = std::chrono::system_clock::now();

    //senk::solver::Gmresm<double>(
    //    val, cind, rptr,
    //    b, x, nrm_b, max_iter/50, 50, N, epsilon);
    senk::solver::IluGmresm<double>(
        val, cind, rptr,
        lval, lcind, lrptr, uval, ucind, urptr,
        b, x, nrm_b, max_iter/50, 50, N, epsilon);
    //senk::solver::AbmcIluGmresm<double>(
    //    val, cind, rptr,
    //    lval, lcind, lrptr, uval, ucind, urptr,
    //    size_color, num_color, 128,
    //    b, x, nrm_b, max_iter/50, 50, N, epsilon);
    //senk::solver::IlubGmresm<double, bnl, bnw>(
    //    val, cind, rptr,
    //    blval, blcind, blrptr, buval, bucind, burptr,
    //    b, x, nrm_b, max_iter/50, 50, N, epsilon);
    //senk::solver::AbmcIlubGmresm<double, bnl, bnw>(
    //    val, cind, rptr,
    //    blval, blcind, blrptr, buval, bucind, burptr,
    //    size_color, num_color, 128,
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
    //    b, x, nrm_b, max_iter/50, 50, N, epsilon);
    //senk::solver::IluGcrm<double>(
    //    val, cind, rptr,
    //    lval, lcind, lrptr, uval, ucind, urptr,
    //    b, x, nrm_b, max_iter/50, 50, N, epsilon);
    //senk::solver::IlubGcrm<double, bnl, bnw>(
    //    val, cind, rptr,
    //    blval, blcind, blrptr, buval, bucind, burptr,
    //    b, x, nrm_b, max_iter/50, 50, N, epsilon);

    auto end = std::chrono::system_clock::now();
    double elapsed 
        = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    printf("%e\n", elapsed);

    senk::test::RelativeResidualError(
        val, cind, rptr, b, x, N);
    
    senk::utils::SafeFree<double>(&x);

    return 0;
}
