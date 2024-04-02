
#include "PCAMatrix.h"
#include <iostream>
#include <iomanip>

PCAMatrix PCAMatrix::average()
{
    PCAMatrix pcaMatrix(this->m_rows, 1);
    double sum = 0;
    for (int i = 0; i < this->m_rows; ++i)
    {
        sum = 0;
        for (int j = 0; j < this->m_cols; ++j)
        {
            sum += this->m_data[i * this->m_cols + j];
        }
        pcaMatrix.m_data[i] = sum / this->m_cols;
    }
    return pcaMatrix;
}

PCAMatrix PCAMatrix::operator*(PCAMatrix &matrix_2)
{


    PCAMatrix &matrix_1 = (*this);
    if (matrix_1.m_cols != matrix_2.m_rows)
        throw std::domain_error("");

    int m = matrix_1.m_cols;
    int n = matrix_1.m_rows;
    int p = matrix_2.m_cols;
    PCAMatrix pcaMatrix(n, matrix_2.m_cols);

    double sum = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j)
        {
            sum = 0;
            for (int k = 0; k < m; ++k)
            {
                sum += matrix_1.position(i, k) * matrix_2.position(k, j);
            }
            pcaMatrix.position(i, j) = sum;
        }
    return pcaMatrix;
}

PCAMatrix PCAMatrix::T_Matrix()
{
    PCAMatrix pcaMatrix(this->m_cols, this->m_rows);

    for (int i = 0; i < this->m_cols * this->m_rows; ++i)
        pcaMatrix.m_data[i / this->m_cols + this->m_rows * (i % this->m_cols)] = this->m_data[i];
    return pcaMatrix;
}

PCAMatrix PCAMatrix::operator/(float value)
{
    PCAMatrix pcaMatrix(this->m_rows, this->m_cols);
    for (int i = 0; i < this->m_cols * this->m_rows; ++i)
        pcaMatrix.m_data[i] = this->m_data[i] / value;
    return pcaMatrix;
}

PCAMatrix PCAMatrix::operator*(float value)
{
    PCAMatrix pcaMatrix(this->m_rows, this->m_cols);

    for (int i = 0; i < this->m_cols * this->m_rows; ++i)
        pcaMatrix.m_data[i] = this->m_data[i] * value;
    return pcaMatrix;
}

PCAMatrix PCAMatrix::covariance()
{
    int rows = this->m_rows;
    int cols = this->m_cols;

    PCAMatrix pcaMatrix(rows, cols);

    PCAMatrix avg = this->average();

    for (int band = 0; band < rows; ++band)
    {
        for (int i = 0; i < cols; ++i)
        {
            pcaMatrix.m_data[band * cols + i] = this->m_data[band * cols + i] - avg.m_data[band];
        }
    }

    PCAMatrix T_pcaMatrix = pcaMatrix.T_Matrix();
    return pcaMatrix * T_pcaMatrix / (cols - 1);
}


PCAMatrix &PCAMatrix::operator=(const PCAMatrix &matrix)
{
    const_cast <PCAMatrix &>(matrix).duplicate(*this);
    return *this;
}

PCAMatrix PCAMatrix::operator+(PCAMatrix &matrix)
{
    int rows = this->m_rows;
    int cols = this->m_cols;

    if (rows != matrix.m_rows || cols != matrix.m_cols)
    {
        throw std::domain_error("");
    }

    PCAMatrix pcaMatrix(rows, cols);
    for (int i = 0; i < this->m_cols * this->m_rows; ++i)
        pcaMatrix.m_data[i] = this->m_data[i] + matrix.m_data[i];
    return pcaMatrix;
}

PCAMatrix PCAMatrix::operator+=(PCAMatrix &matrix)
{
    int rows = this->m_rows;
    int cols = this->m_cols;
    if (rows != matrix.m_rows || cols != matrix.m_cols)
    {
        throw std::domain_error("");
    }

    for (int i = 0; i < rows * cols; ++i)
        this->m_data[i] += matrix.m_data[i];
    return (*this);

}