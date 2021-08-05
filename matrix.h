#pragma once
#include <cassert> 
#include <cmath>
#include <thread>
#include "threadMatrix.h"
class matrix {
private:

    long double** mtx;

    size_t n, m;

    size_t numOfThreads;

    std::vector<cell> splitMtx;

    void allocateMem();

public:


    matrix(size_t N, size_t M);

    matrix(size_t N, size_t M, size_t num);

    matrix(size_t N);

    matrix(const matrix& A);

    friend const matrix operator + (const matrix& A, const matrix& B);

    friend const matrix operator - (const matrix& A, const matrix& B);

    friend const matrix operator * (const matrix& A, const matrix& B);

    friend const matrix operator * (double val, const matrix& B);

    matrix& operator = (const matrix& A);

    long double* operator [] (size_t i);
    
    void pullRandDecimal();

    void pullRandSymDecimal();

    void swapRows(size_t i, size_t j);

    double determinate();

    void print(std::ostream& streamOut);

    void print();

    friend double normInf(const matrix& A);

    friend const matrix toUpperTriangle(const matrix& A);

    friend const matrix toUpperTriangle(const matrix& A, matrix& B);

    friend const matrix toLowerTriangle(const matrix& A); 

    friend const matrix toLowerTriangle(const matrix& A, matrix& B);

    friend const matrix identity(size_t n);

    friend const matrix transpose(const matrix& A);

    friend const matrix inverseIter(const matrix& A);

    friend const matrix inverse(const matrix& A); 

    friend const matrix solve(const matrix& A, const matrix& B); 

    ~matrix(); 
};

const matrix identity(size_t n);
