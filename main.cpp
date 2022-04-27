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
    printf("# %d\n", trptr[N]);
    //if(!senk::matrix::CheckStructure(tval, tcind, trptr, N)) {
    //    printf("False: CheckStructure\n"); exit(1);
    //}
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
    M = N;
    printf("%d %d\n", N, M);

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

    double *x = senk::utils::SafeMalloc<double>(N);
    double *b = senk::utils::SafeMalloc<double>(N);
    for(int i=0; i<N; i++) {
        x[i] = 0;
        b[i] = 1;
    }

    int max_iter = 6000;
    int gs_iter = 1;
    int cycle_iter = 1;
    double epsilon = 1.0e-8;
    double nrm_b = senk::blas1::Nrm2(b, N);  

    /*
    senk::solver::MgGsCg(
        lu_val, lu_cind, lu_rptr, diag,
        lu_c_val, lu_c_cind, lu_c_rptr, c_diag,
        i_val, i_cind, i_rptr,
        r_val, r_cind, r_rptr,
        b, x, nrm_b,
        max_iter, N, dim, depth,
        gs_iter, cycle_iter, epsilon, true);
    */

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

    /*
    senk::solver::MgGsGcrm(
        lu_val, lu_cind, lu_rptr, diag,
        lu_c_val, lu_c_cind, lu_c_rptr, c_diag,
        i_val, i_cind, i_rptr,
        r_val, r_cind, r_rptr,
        b, x, nrm_b,
        max_iter/10, 10, N, dim, depth,
        gs_iter, cycle_iter, epsilon, true);
    */
/*
    double *csc_val;
    int *csc_rind;
    int *csc_cptr;
    senk::matrix::Csr2Csc(val, cind, rptr, &csc_val, &csc_rind, &csc_cptr, N, N);

    double *zt_val;
    int *zt_cind;
    int *zt_rptr;
    double *wt_val;
    int *wt_cind;
    int *wt_rptr;
    double *diag;

    senk::ainv::left_painv(
        val, cind, rptr,
        csc_val, csc_rind, csc_cptr,
        &zt_val, &zt_cind, &zt_rptr,
        &wt_val, &wt_cind, &wt_rptr,
        &diag, N, 0.1, 1.0);
*/
/*
    senk::ainv::left_ainv(
        val, cind, rptr,
        csc_val, csc_rind, csc_cptr,
        &zt_val, &zt_cind, &zt_rptr,
        &wt_val, &wt_cind, &wt_rptr,
        &diag, N, 0.1);
*/
/*
    double *z_val;
    int *z_cind;
    int *z_rptr;
    senk::matrix::Csr2Csc(zt_val, zt_cind, zt_rptr, &z_val, &z_cind, &z_rptr, N, N);
    for(int i=0; i<N; i++) {
        diag[i] = 1.0 / diag[i];
    }

    free(csc_val);
    free(csc_rind);
    free(csc_cptr);

    double **zt_c_val = senk::utils::SafeMalloc<double*>(depth);
    int **zt_c_cind   = senk::utils::SafeMalloc<int*>(depth);
    int **zt_c_rptr   = senk::utils::SafeMalloc<int*>(depth);
    double **wt_c_val = senk::utils::SafeMalloc<double*>(depth);
    int **wt_c_cind   = senk::utils::SafeMalloc<int*>(depth);
    int **wt_c_rptr   = senk::utils::SafeMalloc<int*>(depth);
    double **c_diag   = senk::utils::SafeMalloc<double*>(depth);
    double **z_c_val  = senk::utils::SafeMalloc<double*>(depth);
    int **z_c_cind    = senk::utils::SafeMalloc<int*>(depth);
    int **z_c_rptr    = senk::utils::SafeMalloc<int*>(depth);
    for(int i=0; i<depth; i++) {
        senk::matrix::Csr2Csc(c_val[i], c_cind[i], c_rptr[i], &csc_val, &csc_rind, &csc_cptr, dim[i], dim[i]);
        senk::ainv::left_painv(
            c_val[i], c_cind[i], c_rptr[i],
            csc_val, csc_rind, csc_cptr,
            &zt_c_val[i], &zt_c_cind[i], &zt_c_rptr[i],
            &wt_c_val[i], &wt_c_cind[i], &wt_c_rptr[i],
            &c_diag[i], dim[i], 0.1, 1.0);
        //senk::ainv::left_ainv(
        //    c_val[i], c_cind[i], c_rptr[i],
        //    csc_val, csc_rind, csc_cptr,
        //    &zt_c_val[i], &zt_c_cind[i], &zt_c_rptr[i],
        //    &wt_c_val[i], &wt_c_cind[i], &wt_c_rptr[i],
        //    &c_diag[i], dim[i], 0.05);
        senk::matrix::Csr2Csc(
            zt_c_val[i], zt_c_cind[i], zt_c_rptr[i],
            &z_c_val[i], &z_c_cind[i], &z_c_rptr[i], dim[i], dim[i]);
        for(int j=0; j<dim[i]; j++) { c_diag[i][j] = 1 / c_diag[i][j]; }
        free(csc_val);
        free(csc_rind);
        free(csc_cptr);
    }
*/
/*
    float *fval = senk::utils::SafeMalloc<float>(rptr[N]);
    float *fb = senk::utils::SafeMalloc<float>(N);
    float *fx = senk::utils::SafeMalloc<float>(N);
    senk::utils::Convert<double, float>(val, fval, rptr[N]);
    senk::utils::Convert<double, float>(b, fb, N);
    senk::utils::Convert<double, float>(x, fx, N);
    float fnrm_b = senk::blas1::Nrm2(fb, N);
    float fepsilon = 1.0e-8;
    //senk::solver::Gmresm(
    //    val, fval, cind, rptr,
    //    b, x, nrm_b, 10, 20, N, epsilon);
*/
/*
    senk::utils::Set<double>(0, x, N);
    senk::solver::InvGmresm(
        val, cind, rptr,
        z_val, z_cind, z_rptr,
        wt_val, wt_cind, wt_rptr, diag,
        b, x, nrm_b, 100, 20, N, epsilon);
*/
/*
    senk::solver::InvGmresm(
        c_val[0], c_cind[0], c_rptr[0],
        z_c_val[0], z_c_cind[0], z_c_rptr[0],
        wt_c_val[0], wt_c_cind[0], wt_c_rptr[0], c_diag[0],
        b, x, nrm_b, 100, 20, dim[0], epsilon);
*/
/*
    senk::solver::MgInvGmresm(
        val, cind, rptr,
        wt_val, wt_cind, wt_rptr,
        z_val, z_cind, z_rptr, diag,
        c_val, c_cind, c_rptr,
        wt_c_val, wt_c_cind, wt_c_rptr,
        z_c_val, z_c_cind, z_c_rptr, c_diag,
        i_val, i_cind, i_rptr,
        r_val, r_cind, r_rptr,
        b, x, nrm_b, 100, 20, N, dim, depth, epsilon);
*/

    //double *t = senk::utils::SafeMalloc<double>(N);
    //double *t2 = senk::utils::SafeMalloc<double>(N);
    //double *t3 = senk::utils::SafeMalloc<double>(N);
/*
    for(int i=0; i<10000; i++) {
        senk::sparse::SpmvCsr(val, cind, rptr, x, t, N);
        senk::blas1::Axpby(1, b, -1, t, N);
        double nrm = senk::blas1::Nrm2(t, N);
        printf("%e\n", nrm/nrm_b);
        senk::solver::InverseSolve(
            wt_val, wt_cind, wt_rptr, z_val, z_cind, z_rptr, diag,
            t, t2, t3, N);
        senk::blas1::Axpy(1, t2, x, N);
    }
*/
    /*
    senk::test::RelativeResidualError(
        val, cind, rptr, b, x, N);
    */
    senk::utils::SafeFree<double>(&x);

    return 0;
}
