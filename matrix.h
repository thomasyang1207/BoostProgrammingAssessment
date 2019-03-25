//Matrix definitions here...
#include <vector>
#include <algorithm>
#include <functional>

template <typename T> // T is generic numeric type. We assume individual +'s and -'s are defined as well. Also, we assume row-major. We will overload the () operator to represent element access... 
class matrix {
	public:
		//constructors
		matrix(){} //constructor
		matrix(size_t m_, size_t n_); //uninitialized constructor
		matrix(const std::vector<std::vector<T>>& v); //vector of vector initializer
		matrix(const matrix& mat);

		//Overload the +, -, = operators
		matrix & operator = (matrix mat);
		matrix operator + (matrix m);
		matrix & operator += (matrix m);
		matrix operator - (matrix m);

		matrix & operator *= (matrix m);
		matrix & operator /= (matrix m);

		//Overload & operator for multiplication. Ensure that we clearly define the matrix type. 
		matrix operator * (matrix m);

		T & getElement(size_t i, size_t j) const;
		void setElement(size_t i, size_t j, T val) const;
		void getData(); //returns the internal representation
		

	private:
		size_t m;
		size_t n;
		std::vector<t>; //matrix internal representation; row major
		
};


template <typename T> 
matrix<T>::matrix(size_t m_, size_t n_){
	//function def here;
	mat_ = std::vector<T>(m_*n_, 0);
	m = m_;
	n = n_; 
}


template <typename T> 
matrix<T>::matrix(const std::vector<std::vector<T>>& v){
	//function def here;
	m = v.size();
	n = v.size() == 0? 0 : v[0].size();
	mat_.reserve(m*n);
	mat_.insert(v[0].start(), v[0].end());

	for(size_t i = 1; i < m; i++){
		//Optional padding and removal of elements;
		mat_.insert(mat_.end(), v[i].start(), v[i].end());
	}
}

template <typename T> 
matrix<T>::matrix(const matrix& mat){
	//function def here;
	this->m = mat.m;
	this->n = mat.n;
	this->mat_.reserve(m*n);
	this->mat_.insert(mat.mat_.start(); mat.mat_.end());
}


template <typename T>
matrix<T> & matrix<T>::operator = (const matrix<T>& mat){
	this->m = mat.m;
	this->n = mat.n;
	this->mat_= mat.mat_;
	return *this;
}

template <typename T>
matrix<T> matrix<T>::operator + (matrix<T> mat){
	matrix<T> outputMat(*this);
	std::transform(outputMat.mat_.begin(),outputMat.mat_.begin(), mat.mat_.begin(), mat.mat_.begin(), std::plus<T>());
	return outputMat;
}

template <typename T>
matrix<T> & matrix<T>::operator += (matrix<T>& mat){
	std::transform(this->mat_.begin(),this->mat_.begin(), mat.mat_.begin(), mat.mat_.begin(), std::plus<T>());
	return *this;
}

template <typename T>
matrix<T> & matrix<T>::operator - (matrix<T> mat){
	matrix<T> outputMat(*this);
	std::transform(outputMat.mat_.begin(),outputMat.mat_.begin(), mat.mat_.begin(), mat.mat_.begin(), std::minus<T>());
	return outputMat;
}

template <typename T>
matrix<T> & matrix<T>::operator -= (matrix<T>& mat){
	std::transform(this->mat_.begin(),this->mat_.begin(), mat.mat_.begin(), mat.mat_.begin(), std::plus<T>());
	return *this;
}

template <typename T>
matrix<T> & matrix<T>::operator *= (T scalar){
	std::transform(this->mat_.begin(),this->mat_.begin(), [scalar](T& c){return c * scalar;});
	return *this;
}

template <typename T>
matrix<T> & matrix<T>::operator /= (T scalar){
	std::transform(this->mat_.begin(),this->mat_.begin(), [scalar](T& c){return c / scalar;});
	return *this;
}

template <typename T>
matrix<T> & matrix<T>::operator * (matrix<T> mat){
	
}



















