#include "Matrix.h"

#include <algorithm>  // std::max, std::min
#include <cmath>  // std::round, std::sqrt
#include <sstream>  // std::ostringstream
#include <string>  // std::string
#include <vector>  // std::vector

Matrix::Matrix() {}

Matrix::Matrix(const Matrix& other) {
	this->values = other.values;
}

Matrix::Matrix(const std::vector<double>& values) {
	this->values = values;
}

Matrix::Matrix(int nRows) {
	int nValues = (nRows - 1) * nRows / 2;
	this->values = std::vector<double>(nValues, NOT_A_NUMBER);
}

void Matrix::setValue(int i, int j, double value) {
	if (i != j) {
		this->values[index(i, j)] = value;
	}
	return;
}

double Matrix::value(int i, int j) const {
	double vij;
	if (i == j) {
		vij = NOT_A_NUMBER;
	} else {
		vij = this->values[index(i, j)];
	}
	return vij;
}

double Matrix::minValue() const {
	double minv = +INF;
	for (int i = 0; i < (int)this->values.size(); i ++) {
		minv = std::min(minv, this->values[i]);
	}
	return minv;
}

double Matrix::maxValue() const {
	double maxv = -INF;
	for (int i = 0; i < (int)this->values.size(); i ++) {
		maxv = std::max(maxv, this->values[i]);
	}
	return maxv;
}

int Matrix::numRows() const {
	int nValues = (int)this->values.size();
	return (1 + (int)std::round(std::sqrt((double)(1 + 8 * nValues)))) / 2;
}

int Matrix::precision() const {
	std::ostringstream oss;
	oss.precision(MAX_DIGITS);  // Modify the default precision
	int maxDecimals = 0;
	for (int i = 0; i < (int)this->values.size(); i ++) {
		oss.str("");  // Clear string stream
		oss << this->values[i];
		std::string s = oss.str();
		std::size_t found = s.find('.');
		int decimals = (found == std::string::npos)?
				0 : (int)(s.size() - found) - 1;
		maxDecimals = std::max(maxDecimals, decimals);
	}
	return maxDecimals;
}

int Matrix::index(int i, int j) const {
	int k;
	int nRows = numRows();
	if (i == j) {
		k = -1;
	} else if (i > j) {
		k = i + j * nRows - (j + 1) * (j + 2) / 2;
	} else {  // (i < j)
		k = j + i * nRows - (i + 1) * (i + 2) / 2;
	}
	return k;
}
