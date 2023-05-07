#ifndef MATRIX_H_
#define MATRIX_H_

#include <limits>  // std::numeric_limits
#include <vector>  // std::vector

const double INF = std::numeric_limits<double>::infinity();
const int MAX_DIGITS = std::numeric_limits<double>::digits10;
const double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

// Symmetric matrix with null diagonal values
class Matrix {
public:
    Matrix();
    Matrix(const Matrix& other);
    Matrix(const std::vector<double>& values);
    Matrix(int nRows);
    void setValue(int i, int j, double value);
    double value(int i, int j) const;
    double minValue() const;
    double maxValue() const;
    int numRows() const;
    int precision() const;
private:
    std::vector<double> values;  // Lower triangular values by columns
    int index(int i, int j) const;
};

#endif /* MATRIX_H_ */
