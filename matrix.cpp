#pragma once
#include <iostream>
#include <mutex>
#include "matrix.h"

void matrix::allocateMem()
{
    this->mtx = new long double* [this->n];
    for (size_t i(0); i < this->n; i++)
        this->mtx[i] = new long double[this->m];
}

matrix::matrix(size_t N, size_t M) : n(N), m(M) {
    numOfThreads = std::thread::hardware_concurrency() - 1;
    splitMtx = createCellForThreads(N, M, numOfThreads);
    allocateMem();
}

matrix::matrix(size_t N, size_t M, size_t num) : n(N), m(M) {
    numOfThreads = num;
    splitMtx = createCellForThreads(N, M, numOfThreads);
    allocateMem();
}

matrix::matrix(size_t N) : n(N), m(N) {
    numOfThreads = std::thread::hardware_concurrency() - 1;
    splitMtx = createCellForThreads(N, N, numOfThreads);
    allocateMem();
}

matrix::matrix(const matrix& A)
{
    numOfThreads = A.numOfThreads;
    splitMtx = A.splitMtx;
    n = A.n;
    m = A.m;
    allocateMem();
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
            mtx[i][j] = A.mtx[i][j];
    }
}

void matrix::swapRows(size_t i, size_t j)
{
    for (int k = 0; k < m; ++k)
        std::swap(mtx[i][k], mtx[j][k]);
}

const matrix operator + (const matrix& A, const matrix& B)
{
    assert(A.n == B.n && A.m == B.m && "Error in OP (+): dimmensions must be equal");
    matrix R(A.n, A.m);
    std::vector<std::thread> threads;
    auto func = [&](int startRow, int stopRow, int startCol, int stopCol)
    {
        for (int i = startRow; i < stopRow; ++i)
        {
            for (int j = startCol; j < stopCol; ++j)
                R[i][j] = A.mtx[i][j] + B.mtx[i][j];
        }
    };
    std::vector<cell> splitMtx(A.splitMtx);
    for (size_t q = 0; q < A.numOfThreads; ++q)
    {
        std::thread tmp(func, splitMtx[q].left.x, splitMtx[q].right.x, splitMtx[q].left.y, splitMtx[q].right.y);
        threads.push_back(std::move(tmp));
    }
    for (std::thread& th : threads)
    {
        if (th.joinable())
            th.join();
    }
    return R;
}

const matrix operator - (const matrix& A, const matrix& B)
{
    assert(A.n == B.n && A.m == B.m && "Error in OP (+): dimmensions must be equal");
    matrix R(A.n, A.m);
    std::vector<std::thread> threads;
    auto func = [&](int startRow, int stopRow, int startCol, int stopCol)
    {
        for (int i = startRow; i < stopRow; ++i)
        {
            for (int j = startCol; j < stopCol; ++j)
                R[i][j] = A.mtx[i][j] - B.mtx[i][j];
        }
    };
    std::vector<cell> splitMtx(A.splitMtx);
    for (size_t q = 0; q < A.numOfThreads; ++q)
    {
        std::thread tmp(func, splitMtx[q].left.x, splitMtx[q].right.x, splitMtx[q].left.y, splitMtx[q].right.y);
        threads.push_back(std::move(tmp));
    }
    for (std::thread& th : threads)
    {
        if (th.joinable())
            th.join();
    }
    return R;
}

const matrix operator * (const matrix& A, const matrix& B)
{
    assert(A.m == B.n && "Error in OP (+): dimmensions must be equal");
    matrix R(A.n, B.m, std::min(A.numOfThreads, B.numOfThreads));
    std::vector<std::thread> threads;
    auto func = [&](int startRow, int stopRow, int startCol, int stopCol)
    {
        for (int i = startRow; i < stopRow; ++i)
        {
            for (int j = startCol; j < stopCol; ++j)
            {
                long double tmp = 0;
                for (int k = 0; k < B.n; ++k)
                    tmp += A.mtx[i][k] * B.mtx[k][j];
                R[i][j] = tmp;
            }
        }
    };
    for (size_t q = 0; q < R.numOfThreads; ++q)
    {
        std::thread tmp(func, R.splitMtx[q].left.x, R.splitMtx[q].right.x, R.splitMtx[q].left.y, R.splitMtx[q].right.y);
        threads.push_back(std::move(tmp));
    }
    for (std::thread& th : threads)
    {
        if (th.joinable())
            th.join();
    }
    return R;
}
const matrix operator * (double val, const matrix& B)
{
    matrix R(B);
    std::vector<std::thread> threads;
    auto func = [&](int startRow, int stopRow, int startCol, int stopCol)
    {
        for (int i = startRow; i < stopRow; ++i)
        {
            for (int j = startCol; j < stopCol; ++j)
            {
                R.mtx[i][j] *= val;
            }
        }
    };
    for (size_t q = 0; q < R.numOfThreads; ++q)
    {
        std::thread tmp(func, R.splitMtx[q].left.x, R.splitMtx[q].right.x, R.splitMtx[q].left.y, R.splitMtx[q].right.y);
        threads.push_back(std::move(tmp));
    }
    for (std::thread& th : threads)
    {
        if (th.joinable())
            th.join();
    }
    return R;
}

matrix& matrix::operator = (const matrix& A)
{
    this->numOfThreads = A.numOfThreads; 
    this->splitMtx = A.splitMtx; 
    assert(A.n == n && A.m == m && "Error in OP(=): dimmension must be equal\n");
    std::vector<std::thread> threads;
    auto func = [&](int startRow, int stopRow, int startCol, int stopCol)
    {
        for (int i = startRow; i < stopRow; ++i)
        {
            for (int j = startCol; j < stopCol; ++j)
            {
                mtx[i][j] = A.mtx[i][j];
            }
        }
    };
    for (size_t q = 0; q < numOfThreads; ++q)
    {
        std::thread tmp(func, A.splitMtx[q].left.x, A.splitMtx[q].right.x, A.splitMtx[q].left.y, A.splitMtx[q].right.y);
        threads.push_back(std::move(tmp));
    }
    for (std::thread& th : threads)
    {
        if (th.joinable())
            th.join();
    }
    return *this;
}

long double* matrix::operator [] (size_t i)
{
    return mtx[i];
}

const matrix toUpperTriangle(const matrix& A)
{
    matrix res(A);
    for (int i = 0; i < std::min(A.n, A.m); ++i)
    {
        double maxKK = std::fabs(res[i][i]);
        int indMax = i;
        for (int j = i + 1; j < A.n; ++j)
        {
            if (maxKK < std::fabs(res[j][i]))
            {
                indMax = j;
                maxKK = std::fabs(res[j][i]);
            }
        }
        if (indMax != i)
        {
            for (int j = 0; j < A.m; ++j)
                res[i][j] *= -1;
            res.swapRows(i, indMax);
        }
        if (res[i][i] == 0)
            continue; 
        for (int j = i + 1; j < A.n; ++j)
        {
            double coef = -res[j][i] / res[i][i];
            res[j][i] = 0;
            for (int k = i + 1; k < A.m; ++k)
            {
                res[j][k] = coef * res[i][k] + res[j][k];
            }
        }
    }
    return res;
}

const matrix toLowerTriangle(const matrix& A)
{
    matrix res(A);
    for (int i = std::min(A.n, A.m)-1; i >=0; --i)
    {
        if (res[i][i] == 0)
            continue; 
        for (int j = i - 1; j>=0; --j)
        {
            double coef = -res[j][i] / res[i][i];
            res[j][i] = 0;
            for (int k = i + 1; k < A.m; ++k)
            {
                res[j][k] = coef * res[i][k] + res[j][k];
            }
        }
    }
    return res;
}

const matrix toLowerTriangle(const matrix& A, matrix& B)
{
    matrix res(A);
    for (int i = std::min(A.n, A.m) - 1; i >= 0; --i)
    {
        if (res[i][i] == 0)
            continue;
        for (int j = i - 1; j >= 0; --j)
        {
            double coef = -res[j][i] / res[i][i];
            res[j][i] = 0;
            B[j][i] = coef * B[i][i] + B[j][i];
            for (int k = i + 1; k < A.m; ++k)
            {
                res[j][k] = coef * res[i][k] + res[j][k];
                B[j][k] = coef * B[i][k] + B[j][k]; 
            }
        }
    }
    return res;
}

const matrix toUpperTriangle(const matrix& A, matrix& B)
{
    matrix res(A);
    for (int i = 0; i < std::min(A.n, A.m); ++i)
    {
        double maxKK = std::fabs(res[i][i]);
        int indMax = i;
        for (int j = i + 1; j < A.n; ++j)
        {
            if (maxKK < std::fabs(res[j][i]))
            {
                indMax = j;
                maxKK = std::fabs(res[j][i]);
            }
        }
        if (indMax != i)
        {
            for (int j = 0; j < A.m; ++j)
            {
                res[i][j] *= -1;
                B[i][j] *= -1; 
            }
            res.swapRows(i, indMax);
            B.swapRows(i, indMax); 
        }
        for (int j = i + 1; j < A.n; ++j)
        {
            double coef = -res[j][i] / res[i][i];
            res[j][i] = 0;
            for (int k=0; k<A.m;++k)
                B[j][k] = coef * B[i][k] + B[j][k];
            for (int k = i + 1; k < A.m; ++k) 
                res[j][k] = coef * res[i][k] + res[j][k];
        }
    }
    return res;
}

const matrix identity(size_t n)
{
    matrix R(n);
    std::vector<std::thread> threads;
    auto func = [&](int startRow, int stopRow, int startCol, int stopCol) {
        for (int i = startRow; i < stopRow; ++i)
        {
            for (int j = startCol; j < stopCol; ++j)
            {
                if (i == j)
                    R[i][j] = 1;
                else
                    R[i][j] = 0;
            }
        }
    };
    for (int q = 0; q < R.numOfThreads; ++q)
    {
        std::thread tmp(func, R.splitMtx[q].left.y, R.splitMtx[q].right.y, R.splitMtx[q].left.x, R.splitMtx[q].right.x);
        threads.push_back(std::move(tmp));
    }
    for (std::thread& th : threads)
    {
        if (th.joinable())
            th.join();
    }
    return R;
}

double normInf(const matrix& A)
{
    double max = std::fabs(A.mtx[0][0]);
    for (int i = 0; i < A.n; ++i)
    {
        double sum = 0;
        for (int j = 0; j < A.m; ++j)
            sum += std::fabs(A.mtx[i][j]);
        if (max < sum)
            max = sum;
    }
    return max;
}

const matrix transpose(const matrix& A)
{
    matrix R(A.n);
    for (int i = 0; i < A.n; ++i)
        for (int j = 0; j < A.m; ++j)
            R[i][j] = A.mtx[j][i];
    return R;
}

void matrix::print(std::ostream& streamOut)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            streamOut << mtx[i][j] << " ";
        streamOut << std::endl;
    }
    streamOut << std::endl;
}

void matrix::print()
{
    return print(std::cout); 
}

const matrix inverseIter(const matrix& A)
{
    double eps = 0.000001;
    matrix iden = identity(A.n);
    matrix U0 = eps * transpose(A);
    matrix discrepancy = iden - A * U0;
    matrix Uk(U0);
    while (normInf(discrepancy) > eps)
    {
        Uk = Uk * (iden + discrepancy);
        discrepancy = iden - A * Uk;
    }
    return Uk;
}

const matrix inverse(const matrix& A)
{
    matrix B = identity(A.n); 
    matrix U = toUpperTriangle(A, B);
    matrix R(A.n, A.n, A.numOfThreads); 

    std::vector<std::thread> threads; 
    std::mutex mute; 
    auto func = [&] (int i0){
        for (int i = A.n - 1; i >= 0; --i)
        {
            double sum = B[i][i0]; 
            for (int j = A.n - 1; j > i; --j)
            {
                sum -= R[j][i0] * U[i][j];
            }
            R[i][i0] = sum / U[i][i];  
        }
    };
    for (int i = 0; i < A.n;)
    {
        //func(i); 
        for (int q = 0; q < A.numOfThreads; ++q)
        {
            std::thread tmp(func, i); 
            threads.push_back(std::move(tmp));
            ++i; 
            if (i >= A.n)
                break; 
        }
        for (std::thread &th : threads)
        {
                th.join(); 
        }
        if (i >= A.n)
            break;
        threads.clear(); 
    }
    return R; 
}

void matrix::pullRandDecimal()
{
    srand(clock());
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            mtx[i][j] = rand() % 10;
}

void matrix::pullRandSymDecimal()
{
    assert(n == m && "Impossible to create Symmetric matrix");
    srand(clock());
    for (int i = 0; i < n; ++i)
        for (int j = i; j < m; ++j)
        {
            mtx[i][j] = rand() % 10;
            mtx[j][i] = mtx[i][j];
        }

}

double matrix::determinate()
{
    assert(n == m && "Determinate defines only for square matrix");
    matrix tmp = toUpperTriangle(*this);
    double res = 1;
    for (int i = 0; i < n; i++)
        res *= tmp[i][i];
    return res;
}

const matrix solve(const matrix& A, const matrix& B)
{
    assert(A.n == B.n && "Cannot solve this: dimmensions must be equal"); 
    matrix B0 = B; 
    matrix U = toUpperTriangle(A, B0);
    return A; 
}

matrix::~matrix()
{
    for (int i = 0; i < n; i++)
        delete[](mtx[i]);
    delete[] mtx;
}
