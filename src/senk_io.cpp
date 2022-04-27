#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "senk_utils.hpp"
#include "senk_helper.hpp"
#include "senk_io.hpp"

namespace senk {

namespace io {

using senk::utils::SafeMalloc;
using senk::utils::SafeRealloc;

bool ReadMatrixMarket(
    std::string filename, double **val, int **cind, int **rptr,
    int *N, int *M, Shape *shape, bool removeZeros)
{
    std::cout << "# Readfile MatrixMarket" << std::endl;
    std::ifstream ifs(filename);
    if(!ifs) {
        std::cerr << "# Could not open input file." << std::endl;
        return false;
    }
    std::string line;
    std::getline(ifs, line);
    if(line.find("complex") != std::string::npos) {
        std::cerr << "# Header line is not valid." << std::endl;
        return false;
    }
    if(line.find("symmetric") != std::string::npos) { shape[0] = Sym; }
    else if(line.find("general") != std::string::npos) { shape[0] = Unsym; }
    else { std::cerr << "# Header line is not valid." << std::endl; return false; }

    int PE;
    while(std::getline(ifs, line)) {
        if(line[0] == '%') { continue; }
        std::istringstream iss(line);
        iss >> N[0] >> M[0] >> PE;
        break;
    }
    std::cout << "# " << N[0] << " " << M[0] << " " << PE << std::endl;
    *val     = SafeMalloc<double>(PE); //new double[PE];
    *cind    = SafeMalloc<int>(PE); //new int[PE];
    *rptr    = SafeMalloc<int>(N[0]+1); //new int[N[0]+1];
    int *row = SafeMalloc<int>(PE); //new int[PE];

    int nnz = 0;
    double t_val;
    int t_row, t_col;
    for(int i=0; i<PE; i++) {
        ifs >> t_row >> t_col >> t_val;
        if(removeZeros && t_val == 0) continue;
        (*val)[nnz]  = t_val;
        row[nnz]     = t_row-1;
        (*cind)[nnz] = t_col-1;
        nnz++;
    }
    if(nnz != PE) {
        *val  = SafeRealloc<double>(*val, nnz);
        *cind = SafeRealloc<int>(*cind, nnz);
    }
    helper::QuickSort<int, int, double>(row, *cind, *val, 0, nnz-1);
    int cnt = 0;
    (*rptr)[0] = 0;
    for(int i=0; i<N[0]; i++) {
        int num = 0;
        while(row[cnt] == i) {
            num++; cnt++;
            if(cnt == nnz) { break; }
        }
        (*rptr)[i+1] = (*rptr)[i] + num;
        helper::QuickSort<int, double>(*cind, *val, (*rptr)[i], (*rptr)[i+1]-1);
    }
    std::free(row);
    return true;
}
    
}

}
