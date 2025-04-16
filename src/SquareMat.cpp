// baruh.ifraimov@gmail.com
#include <iostream>
#include <stdexcept>
#include "../inc/SquareMat.hpp"

namespace mtrx{

	// Help function to sum any size of matrix
	double matSum( const SquareMat& o){
		int size = o.getSize();
		double result=0;
		for (int i = 0; i < size ; i++)
		{
			for (int j = 0; j < size; j++)
		{
			result+=(o[i][j]);
		}
		}
		return result;
	}

	SquareMat::SquareMat(int size){
		if(size <= 0){
			throw(std::length_error("Negative (or zero) value will not be tolerated!"));
		}
		this->size = size;
		mat_table = new double*[size];
		for (int i = 0; i < size; i++)
		{
			mat_table[i] = new double[size];
		}

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				mat_table[i][j] = 0;
			}
		}
	}

	SquareMat::SquareMat(double** matrix, int size){
		if(size <= 0){
			throw(std::length_error("Negative (or zero) value will not be tolerated!"));
		}

		if(matrix == nullptr || (*matrix) == nullptr){
			throw(std::invalid_argument("Matrix is null"));
		}

		this->size = size;
		mat_table = new double*[size];
		for (int i = 0; i < size; i++)
		{
			mat_table[i] = new double[size];
		}

		for (int i = 0; i < size; i++)
		{
			if (matrix[i] == nullptr) {
				throw std::invalid_argument("Row pointer is null at index " + std::to_string(i));
			}
			for (int j = 0; j < size; j++)
			{
				mat_table[i][j] = matrix[i][j];
			}
		}
		
	}

	SquareMat::SquareMat(const SquareMat & o){

		this->size = o.size;
		mat_table = new double*[size];
		for (int i = 0; i < size; i++)
		{
			mat_table[i] = new double[size];
		}
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				mat_table[i][j] = o.mat_table[i][j];
			}
			
		}
		
	}

	SquareMat::~SquareMat(){
		if(mat_table != nullptr){
		for (int i = 0; i < this->size; i++)
		{
			delete[] mat_table[i];
		}
		delete[] mat_table;
		}
	}

	// ~m1 (Transpose)
	SquareMat SquareMat::operator~(){
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = this->mat_table[j][i];
		}
		}
		return m;
	}


	// m1 = m2
	SquareMat& SquareMat::operator=(const SquareMat& o) {
		if (this == &o) {
			return *this; // are we the same? then dont do nothing.
		}
		
		// Clean up existing fields (to free up unreachable memmory)
		if (mat_table != nullptr) {
			for (int i = 0; i < size; i++) {
				delete[] mat_table[i];
			}
			delete[] mat_table;
		}
		
		// Copy from other matrix
		size = o.size;
		mat_table = new double*[size];
		for (int i = 0; i < size; i++) {
			mat_table[i] = new double[size];
			for (int j = 0; j < size; j++) {
				mat_table[i][j] = o.mat_table[i][j];
			}
		}
		
		return *this;
	}

	// m1+m2
	SquareMat SquareMat::operator+(const SquareMat& o) const{
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = this->mat_table[i][j] + o.mat_table[i][j];
		}
		}
		return m;
	}

	// m1-m2
	SquareMat SquareMat::operator-(const SquareMat& o) const{
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = this->mat_table[i][j] - o.mat_table[i][j];
		}
		}
		return m;
	}

	// m1*m2
	SquareMat SquareMat::operator*(const SquareMat& o) const{
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		SquareMat m(size);
		double curr_result =0;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				for (int t = 0; t < size; t++)
				{
					curr_result += this->mat_table[i][t] * o.mat_table[t][j];
				}
				m.mat_table[i][j] = curr_result;
				curr_result = 0;
			}
		}
		return m;
	}

	// m1%m2
	SquareMat SquareMat::operator%(const SquareMat& o) const{
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = this->mat_table[i][j] * o.mat_table[i][j];
		}
		}
		return m;
	}

	// -m1
	SquareMat SquareMat::operator-(){
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = -(this->mat_table[i][j]);
		}
		}
		return m;
	}

	// m1*<num>
	SquareMat SquareMat::operator*(double num) const{
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = num*(this->mat_table[i][j]);
		}
		}
		return m;
	}

	// <num>*m1
	SquareMat operator*(double num, const SquareMat& o){
		SquareMat m(o.size);
		for (int i = 0; i < o.size; i++)
		{
			for (int j = 0; j < o.size; j++)
		{
			m.mat_table[i][j] = num*(o.mat_table[i][j]);
		}
		}
		return m;
	}

	// m1^<num>
	SquareMat SquareMat::operator^(int num) const{
		SquareMat m(size);
		m = *this;
		if(num == 0){
			for (int i = 0; i < size; i++)
			{
					for (int j = 0; j < size; j++)
				{
					if(i == j)
					{
						m.mat_table[i][j] = 1;
					}
					else
					{
					m.mat_table[i][j] = 0;
					}
				}
			}
		}
		else{
		for (int i = 1; i < num; i++)
		{
			m*=m;
		}
	}
		return m;
	}

	// m1/<num>
	SquareMat SquareMat::operator/(double num) const{
		if(num == 0){
			throw "Can't divide by zero!";
		}
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = (this->mat_table[i][j])/num;
		}
		}
		return m;
	}

	// m1%<num>
	SquareMat SquareMat::operator%(int num) const{
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = (int)(this->mat_table[i][j])%num;
		}
		}
		return m;
	}

	// m1++
	SquareMat SquareMat::operator++(int){
		SquareMat m = *this;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			this->mat_table[i][j]++;
		}
		}
		return m;
	}

	// ++m1
	SquareMat& SquareMat::operator++(){
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			++(this->mat_table[i][j]);
		}
		}
		return *this;
	}

	// m1--
	SquareMat SquareMat::operator--(int){
		SquareMat m = *this;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			(this->mat_table[i][j])--;
		}
		}
		return m;
	}

	// --m1
	SquareMat& SquareMat::operator--(){
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			--(this->mat_table[i][j]);
		}
		}
		return *this;
	}


	// m1[] const
	SquareMat::HiddenArr SquareMat::operator[](int index) const{
		if (index < 0 || index >= size) {
			throw std::out_of_range("Index is out of bounds");
		}
		return HiddenArr(mat_table[index], size);
	}

	// m1[]
	SquareMat::HiddenArr SquareMat::operator[](int index){
		if (index < 0 || index >= size) {
			throw std::out_of_range("Index is out of bounds");
		}
		return HiddenArr(mat_table[index], size);
	}

	// m1 == m2
	bool SquareMat::operator==(const SquareMat& o) const{
		double result1=0,result2=0;
		result1 = mtrx::matSum(*this);
		result2 = mtrx::matSum(o);
		return result1==result2;
	}

	// m1 != m2
	bool SquareMat::operator!=(const SquareMat& o) const{
		if (this->size != o.size){
			return true;
		}
		return !(*this == o);
	}

	// m1 < m2
	bool SquareMat::operator<(const SquareMat& o) const{
		double result1=0,result2=0;
		result1 = mtrx::matSum(*this);
		result2 = mtrx::matSum(o);
		return result1 < result2;
	}

	// m1 > m2
	bool SquareMat::operator>(const SquareMat& o) const{
		double result1=0,result2=0;
		result1 = mtrx::matSum(*this);
		result2 = mtrx::matSum(o);
		return result1 > result2;
	}

	// m1 <= m2
	bool SquareMat::operator<=(const SquareMat& o) const{
		double result1=0,result2=0;
		result1 = mtrx::matSum(*this);
		result2 = mtrx::matSum(o);
		return result1 <= result2;
	}

	// m1 >= m2
	bool SquareMat::operator>=(const SquareMat& o) const{
		double result1=0,result2=0;
		result1 = mtrx::matSum(*this);
		result2 = mtrx::matSum(o);
		return result1 >= result2;
	}

		// This function uses recursive cofactor expansion.
	int determinantHelper(double **mat, int n) {
		// Base case for 1x1 matrix: the determinant is the single element.
		if (n == 1)
			return (int)mat[0][0];

		// Base case for 2x2 matrix: use the formula: det = ad - bc.
		if (n == 2)
			return (int)(mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

		int det = 0;
		int sign = 1;  // Sign multiplier, alternating between 1 and -1.

		// Allocate memory for a temporary submatrix of size (n-1) x (n-1)
		double **temp = new double*[n - 1];
		for (int i = 0; i < n - 1; i++) {
			temp[i] = new double[n - 1];
		}

		// Cofactor expansion along the first row.
		for (int col = 0; col < n; col++) {
			// Build the submatrix by skipping the first row and the current column 'col'.
			int subi = 0; // Submatrix row index.
			for (int i = 1; i < n; i++) {
				int subj = 0; // Submatrix column index.
				for (int j = 0; j < n; j++) {
					if (j == col)
						continue;  // Skip the column used for cofactor expansion.
					temp[subi][subj] = mat[i][j];
					subj++;
				}
				subi++;
			}
			// Recursively calculate the determinant of the submatrix.
			int subDet = determinantHelper(temp, n - 1);
			// Add the term for the current column: element * subDet * sign.
			det += sign * (int)mat[0][col] * subDet;
			// Alternate the sign for the next term.
			sign = -sign;
		}

		// Free the memory allocated for the temporary submatrix.
		for (int i = 0; i < n - 1; i++) {
			delete [] temp[i];
		}
		delete [] temp;

		return det;
	}

	int SquareMat::determinantHelper(double **mat, int n) const {
		// Base case for 1x1 matrix: the determinant is the single element.
		if (n == 1)
			return (int)mat[0][0];
	
		// Base case for 2x2 matrix: use the formula: det = ad - bc.
		if (n == 2)
			return (int)(mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);
	
		int det = 0;
		int sign = 1;  // Sign multiplier, alternating between 1 and -1.
	
		// Allocate memory for a temporary submatrix of size (n-1) x (n-1)
		double **temp = new double*[n - 1];
		for (int i = 0; i < n - 1; i++) {
			temp[i] = new double[n - 1];
		}
	
		// Cofactor expansion along the first row.
		for (int col = 0; col < n; col++) {
			// Build the submatrix by skipping the first row and the current column 'col'.
			int subi = 0; // Submatrix row index.
			for (int i = 1; i < n; i++) {
				int subj = 0; // Submatrix column index.
				for (int j = 0; j < n; j++) {
					if (j == col)
						continue;  // Skip the column used for cofactor expansion.
					temp[subi][subj] = mat[i][j];
					subj++;
				}
				subi++;
			}
			// Recursively calculate the determinant of the submatrix.
			int subDet = determinantHelper(temp, n - 1);
			// Add the term for the current column: element * subDet * sign.
			det += sign * (int)mat[0][col] * subDet;
			// Alternate the sign for the next term.
			sign = -sign;
		}
	
		// Free the memory allocated for the temporary submatrix.
		for (int i = 0; i < n - 1; i++) {
			delete [] temp[i];
		}
		delete [] temp;
	
		return det;
	}

	// !m1 (Determinant)
	int SquareMat::operator!() const {
		return determinantHelper(this->mat_table, this->size);
	}

	// m1 += m2
	SquareMat& SquareMat::operator+=(const SquareMat& o){
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		*this = *this + o;

		return *this;
	}

	// m1 -= m2
	SquareMat& SquareMat::operator-=(const SquareMat& o){
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		*this = *this - o;

		return *this;
	}

	// m1 /= num
	SquareMat& SquareMat::operator/=(double num){
		*this = (*this)/num;

		return *this;
	}

	// m1 *= m2
	SquareMat& SquareMat::operator*=(const SquareMat& o){
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		*this = *this * o;

		return *this;
	}

	// m1 %= m2
	SquareMat& SquareMat::operator%=(const SquareMat& o){
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		*this = *this % o;

		return *this;
	}
	
	// m1 *= num
	SquareMat& SquareMat::operator*=(double num){
		
		*this = *this * num;

		return *this;
	}

	// m1 %= num
	SquareMat& SquareMat::operator%=(double num){
		*this = *this % num;

		return *this;
	}

	// cout << m1 
	std::ostream& operator<<(std::ostream& os, const SquareMat& o){
		int size = o.size;
		for (int i = 0; i < size; i++)
		{
			os << "[";
			for (int j = 0; j < size-1; j++)
			{
				os << o.mat_table[i][j] << ","; 
			}
			os << o.mat_table[i][size-1]; 
			os << "]\n";
		}
		return os;
	}
}