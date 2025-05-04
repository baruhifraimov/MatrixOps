// baruh.ifraimov@gmail.com
#include <iostream>
#include <stdexcept>
#include "../inc/SquareMat.hpp"
#include "../inc/MatrixUtils.hpp"
namespace mtrx{

	// Help function to sum any size of matrix
	double matSum( const SquareMat& o){
		int size = o.getSize();
		double result=0;
		// Loop through the matrix and sum all the elements
		for (int i = 0; i < size ; i++)
		{
			for (int j = 0; j < size; j++)
		{
			result+=(o[i][j]);
		}
		}
		return result;
	}

	// Constructor 1
	SquareMat::SquareMat(int size){
		// Cannot generate a matrix with size below or equal to zero
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

	// Constructor 2
	SquareMat::SquareMat(double** matrix, int size){
		// Cannot generate a matrix with size below or equal to zero
		if(size <= 0){
			throw(std::length_error("Negative (or zero) value will not be tolerated!"));
		}

		// Check if the matrix is null
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
			// Check if the row pointer is null
			if (matrix[i] == nullptr) {
				throw std::invalid_argument("Row pointer is null at index " + std::to_string(i));
			}
			for (int j = 0; j < size; j++)
			{
				mat_table[i][j] = matrix[i][j];
			}
		}
		
	}

	// Copy Constructor
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

	// Destructor
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
	SquareMat SquareMat::operator~() const{
		SquareMat m(size);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
		{
			m.mat_table[i][j] = this->mat_table[j][i]; // Transpose
		}
		}
		return m;
	}


	// m1 = m2
	SquareMat& SquareMat::operator=(const SquareMat& o) {
		// Check for self-assignment
		if (this == &o) {
			return *this; // are we the same? then return myself.
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
		// Check if the sizes of the matrices are equal else throw an exception
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}

		SquareMat m(size);

		// Loop through the matrices add each element by element
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
		// Check if the sizes of the matrices are equal else throw an exception
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		SquareMat m(size);

		// Loop through the matrices subtract each element by element
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
		// Check if the sizes of the matrices are equal else throw an exception
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}

		SquareMat m(size);

		// Multiply the matrices (standard matrix multiplication)
		double curr_result =0;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				for (int t = 0; t < size; t++)
				{
					curr_result += this->mat_table[i][t] * o.mat_table[t][j];
				}
				// Round to 0 if close to 0 (to avoid floating point errors)
				if (curr_result < 0.0001 && curr_result > 0){
					curr_result =0;
				}
				else if(curr_result > -0.0001 && curr_result < 0){
					curr_result =0;
				}
				m.mat_table[i][j] = curr_result;
				curr_result = 0;
			}
		}
		return m;
	}

	// m1%m2
	SquareMat SquareMat::operator%(const SquareMat& o) const{
		// Check if the sizes of the matrices are equal else throw an exception
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		SquareMat m(size);

		// Loop through the matrices multiply each element by element
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
		
		// each element is negated
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
				// Round to 0 if close to 0 (to avoid floating point errors)
				if(m.mat_table[i][j] < 0.00001 && m.mat_table[i][j] > 0){
				m.mat_table[i][j] =0;
			}
			else if(m.mat_table[i][j] > -0.00001 && m.mat_table[i][j] < 0){
				m.mat_table[i][j] =0;
			}
		}
		}
		return m;
	}

	// <num>*m1
	SquareMat operator*(double num, const SquareMat& o){
		SquareMat m(o.size);

		// Loop through the matix and multiply each matrix element by scalar
		for (int i = 0; i < o.size; i++)
		{
			for (int j = 0; j < o.size; j++)
		{
			m.mat_table[i][j] = num*(o.mat_table[i][j]);
			// Round to 0 if close to 0 (to avoid floating point errors)
			if(m.mat_table[i][j] < 0.00001 && m.mat_table[i][j] > 0){
				m.mat_table[i][j] =0;
			}
			else if(m.mat_table[i][j] > -0.00001 && m.mat_table[i][j] < 0){
				m.mat_table[i][j] =0;
			}
		}
		}
		return m;
	}
	
	// Calculates the determinant of a matrix
	double SquareMat::determinantHelper(double** mat, int n) const {
		if (n == 1) // base case
			return mat[0][0];
		if (n == 2) // base case
			return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
	  
		double det = 0, sign = 1;

		// allocate 'temp' as MINOR matrix
		double** temp = new double*[n-1];
		for(int i=0;i<n-1;i++)
		{
			temp[i] = new double[n-1];
		}

		// Determinant by ROW
		for(int col=0; col<n; col++){
		// Build the MINOR, for example ->  for element mat[0][col]
		// This means deleting row 0 and column 'col'
		  int subi = 0;
		  for(int i=1;i<n;i++){
			int subj=0;
			for(int j=0;j<n;j++){
			  if(j==col) continue;
			  temp[subi][subj++] = mat[i][j];
			}
			subi++;
		  }

		// Compute cofactor = sign * element * det(minor)
		  det += sign * mat[0][col] * determinantHelper(temp, n-1);

		// Alternate sign for the next cofactor
		  sign = -sign;
		}
	  
		for(int i=0;i<n-1;i++) delete[] temp[i];
		delete[] temp;
		return det;
	  }
	  

	// !m1 (Determinant)
	double SquareMat::operator!() const {
		return determinantHelper(mat_table, size);
	  }

	// adj (adjugate) of a matrix
	SquareMat SquareMat::adj() const {
		SquareMat adj(size);
		// 1x1 matrix the adj is 1
		if (size == 1) {
			adj.mat_table[0][0] = 1;
			return adj;
		}
		
		// temp matrix for storing cofactors
		double** temp = new double*[size-1];
		for (int i = 0; i < size-1; i++) {
			temp[i] = new double[size-1];
		}
		
		// Calculate the adj matrix
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				// cofactor calculation
				int temp_i = 0, temp_j = 0;
				
				// Looping for each element of the matrix
				for (int row = 0; row < size; row++) {
					for (int col = 0; col < size; col++) {
						if (row != i && col != j) {
							temp[temp_i][temp_j++] = mat_table[row][col];
							// Row is filled, so increase row index and reset col index
							if (temp_j == size - 1) {
								temp_j = 0;
								temp_i++;
							}
						}
					}
				}
				
				// Calculate sign and determinant for adjugate
				int sign = ((i + j) % 2 == 0) ? 1 : -1;
				adj.mat_table[j][i] = sign * determinantHelper(temp, size-1);
			}
		}
		
		// Clean up the temporary matrix
		for (int i = 0; i < size-1; i++) {
			delete[] temp[i];
		}
		delete[] temp;
		
		return adj;
	}

	
	// Calculates the inverse of a matrix
	// based on https://www.mathsisfun.com/algebra/matrix-inverse.html
	SquareMat SquareMat::inverse() const {
		// compute determinant (operator! now returns double)
		double det = !(*this);
		if (det == 0.0) {
			// if det equals to 0 the matrix is not invertible
			throw std::invalid_argument("Matrix is not invertible (determinant is zero)");
		}
	
		// get adjugate (you already have adj())
		SquareMat adj = this->adj();
	
		// scale each element by 1/det
		double invDet = 1.0 / det;
		adj *= invDet;
	
		return adj;
	}

	// m1^<num>
	SquareMat SquareMat::operator^(int num) const{
		SquareMat m(size);
		m = *this;

		if(num <0){
			m = m.inverse();
			SquareMat copyM = m;
			for (int i = 1; i < -num; i++)
			{
				m*=copyM;
			}
		}
		else if(num == 0){
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
			m*=(*this);
		}
	}

		return m;
	}

	// m1/<num>
	SquareMat SquareMat::operator/(double num) const{
		// cannot divide by 0, throw an exception
		if(num == 0){
			throw std::invalid_argument("Can't divide by zero!");
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
		// Check if the index is out of bounds
		if (index < 0 || index >= size) {
			throw std::out_of_range("Index is out of bounds");
		}
		return HiddenArr(mat_table[index], size);
	}

	// m1[]
	SquareMat::HiddenArr SquareMat::operator[](int index){
		// Check if the index is out of bounds
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

	// m1 += m2
	SquareMat& SquareMat::operator+=(const SquareMat& o){
		// Check if the sizes of the matrices are equal else throw an exception
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		*this = *this + o;

		return *this;
	}

	// m1 -= m2
	SquareMat& SquareMat::operator-=(const SquareMat& o){
		// Check if the sizes of the matrices are equal else throw an exception
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
		// Check if the sizes of the matrices are equal else throw an exception
		if (this->size != o.size){
			throw std::invalid_argument("Matrix sizes do not match");
		}
		*this = *this * o;

		return *this;
	}

	// m1 %= m2
	SquareMat& SquareMat::operator%=(const SquareMat& o){
		// Check if the sizes of the matrices are equal else throw an exception
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