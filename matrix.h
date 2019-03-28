#include <vector>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <iostream>

template <typename E, typename T>
class matrixExpression {
	public:
		T operator() (size_t i, size_t j) const {return static_cast<const E&>(*this)(i, j);}
		size_t nrow() const {return static_cast<E const&>(*this).nrow();}
		size_t ncol() const {return static_cast<E const&>(*this).ncol();}
}; 




template <typename T> 
class matrix : public matrixExpression<matrix<T>, T>{

	public:
		matrix(){}

		matrix(size_t m_, size_t n_){
			mat_ = std::vector<T>(m_*n_, 0);
			m = m_;
			n = n_; 
		}

		matrix(const std::vector<std::vector<T>>& v){
			m = v.size();
			n = v.size() == 0? 0 : v[0].size();
			mat_.reserve(m*n);
			mat_.insert(mat_.end(), v[0].begin(), v[0].end());

			for(size_t i = 1; i < m; i++){
			mat_.insert(mat_.end(), v[i].begin(), v[i].end());
			}
		}

		template <typename E>
		matrix(const matrixExpression<E, T> & mat){
			this->m = mat.nrow();
			this->n = mat.ncol();
			this->mat_ = std::vector<T>(this->m*this->n, 0);

			
			for(size_t i = 0; i < m; i++){
				for (size_t j = 0; j < n; j++){
					mat_[i*n + j] = mat(i,j);
				}
			}
		}


		T operator() (size_t i, size_t j) const {return mat_[i*n + j];}
		T & operator() (size_t i, size_t j) {return mat_[i *n + j];}

		size_t nrow() const {return m;}
		size_t ncol() const {return n;}

		//Overload the +, -, = operators
		template <typename E>
		matrix & operator = (const matrixExpression<E, T> & mat) {
			m = mat.nrow();
			n = mat.ncol();
			mat_ = std::vector<T>(m*n, 0);
			for(size_t i = 0; i < m; i++){
				for (size_t j = 0; j < n; j++){
					mat_[i*n + j] = mat(i,j);
				}
			}
			return *this;
		}

		template <typename E>
		matrix & operator += (const matrixExpression<E, T> & mat) { 
			if(mat.ncol() != this->n || mat.nrow() != this->m)
				throw std::length_error("Error in matrix.h: dimensions of matrices do not match");
			for(size_t i = 0; i < m; i++){
				for (size_t j = 0; j < n; j++){
					mat_[i*n + j] += mat(i,j);
				}
			}
			return *this;
		}

		template <typename E>
		matrix & operator -= (const matrixExpression<E, T> & mat) {
			if(mat.ncol() != this->n || mat.nrow() != this->m)
				throw std::length_error("Error in matrix.h: dimensions of matrices do not match");

			for(size_t i = 0; i < m; i++){
				for (size_t j = 0; j < n; j++){
					mat_[i*n + j] -= mat(i,j);
				}
			}
			return *this;
		}

		matrix & operator *= (T scalar) { 
			for(size_t i = 0; i < m*n; i++) mat_[i] *= scalar; 
			return *this;
		}

		matrix & operator /= (T scalar) { 
			for(size_t i = 0; i < m*n; i++) mat_[i] /= scalar;
				return *this;
		}

		void disp() const {
			for(size_t i = 0; i < m; i++){
				for(size_t j = 0; j < n; j++){
					std::cout << mat_[i*n + j] << " ";
				}
				std::cout << std::endl;
			}
		}

	private:
		size_t m;
		size_t n;
		std::vector<T> mat_;
		
};




template <typename E1, typename E2, typename T> 
class matrixBinOp : public matrixExpression<matrixBinOp<E1, E2, T>, T> {
	const matrixExpression<E1, T> & mat1;
	const matrixExpression<E2, T> & mat2;

	std::function<T(T,T)> op;

	public:
		matrixBinOp(const matrixExpression<E1, T> & m1, const matrixExpression<E1, T> & m2, std::function<T(T,T)> op_) : mat1(m1), mat2(m2), op(op_){}
		T operator() (size_t i, size_t j) const {return op(mat1(i,j), mat1(i,j));}
		size_t nrow() const {return this->mat1.nrow();}
		size_t ncol() const {return this->mat1.ncol();}
};




template <typename E1, typename E2, typename T>
class matrixProd : public matrixExpression <matrixProd<E1, E2, T> ,T> {
	const matrixExpression<E1, T> & mat1;
	const matrixExpression<E2, T> & mat2;

	public:
		matrixProd(const matrixExpression<E1, T> & m1, const matrixExpression<E2, T> & m2) : mat1(m1), mat2(m2) {}
		T operator() (size_t i, size_t j) const {
			T res = mat1(i, 0) * mat2(0, j);
			for(size_t k = 1; k < mat1.ncol(); k++){
				res += mat1(i, k) * mat2(k, j);
			}
			return res;
		}

		size_t nrow() const {return this->mat1.nrow();}
		size_t ncol() const {return this->mat2.ncol();}

};




template <typename E1, typename E2, typename T>
matrixBinOp<E1, E2, T> operator + (const matrixExpression<E1, T> & mat1, const matrixExpression<E2, T> & mat2) {
	if(mat2.ncol() != mat1.ncol() || mat2.nrow() != mat1.nrow())
		throw std::length_error("Error in matrix.h: dimensions of matrices do not match");
	return matrixBinOp<E1, E2, T>(mat1, mat2, std::plus<T>());
}



template <typename E1, typename E2, typename T>
matrixBinOp<E1, E2, T> operator - (const matrixExpression<E1, T> & mat1, const matrixExpression<E2, T> & mat2) {
	if(mat2.ncol() != mat1.ncol() || mat2.nrow() != mat1.nrow())
		throw std::length_error("Error in matrix.h: dimensions of matrices do not match");
	return matrixBinOp<E1, E2, T>(mat1, mat2, std::minus<T>());
}



template <typename E1, typename E2, typename T>
matrixProd<E1, E2, T> operator * (const matrixExpression<E1, T> & mat1, const matrixExpression<E2, T> & mat2) {
	if(mat2.nrow() != mat1.ncol())
		throw std::length_error("Error in matrix.h: dimensions of matrices do not match");

	return matrixProd<E1, E2, T>(mat1, mat2); 
}


















