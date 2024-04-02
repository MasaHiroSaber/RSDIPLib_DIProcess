#include "RSDIPLib.h"
#include <stdexcept>

using namespace std;
class PCAMatrix :public CMatrix
{
public:
    PCAMatrix(CMatrix& matrix) { matrix.duplicate(*this); }
    PCAMatrix(PCAMatrix& matrix) { matrix.duplicate(*this); }
    PCAMatrix(PCAMatrix&& matrix) { matrix.duplicate(*this); }
    PCAMatrix(int rows, int cols) { this->create(rows, cols); }
    PCAMatrix() { }

    PCAMatrix operator*(PCAMatrix& matrix_2);

    PCAMatrix operator+(PCAMatrix& matrix);

    PCAMatrix operator+=(PCAMatrix& matrix);

    PCAMatrix& operator=(const PCAMatrix& matrix);

    PCAMatrix operator/(float value);

    PCAMatrix operator*(float value);

    PCAMatrix average();

    PCAMatrix T_Matrix();

    PCAMatrix covariance();

private:
    inline float& position(int row, int col) { return this->m_data[row * this->m_cols + col]; }
};
