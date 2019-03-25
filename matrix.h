#include <vector>
#include <algorithm>
#include <functional>

template <typename T>
class matrix {
	public:
		matrix(){}
		matrix(size_t m_, size_t n_);
		matrix(const std::vector<std::vector<T>>& v);
		matrix(const matrix& mat);

		//Overload the +, -, = operators
		matrix & operator = (const matrix & mat);
		matrix operator + (const matrix & m) const;
		matrix & operator += (const matrix & m);
		matrix operator - (const matrix & m) const;
		matrix & operator -= (const matrix & m);

		matrix & operator *= (T scalar);
		matrix & operator /= (T scalar);

		matrix operator * (const matrix<T>& mat) const;

		T & getElement(size_t i, size_t j) const;
		void setElement(size_t i, size_t j, T val);

	private:
		size_t m;
		size_t n;
		std::vector<T> mat_;
		
};


template <typename T> 
matrix<T>::matrix(size_t m_, size_t n_){
	mat_ = std::vector<T>(m_*n_, 0);
	m = m_;
	n = n_; 
}


template <typename T> 
matrix<T>::matrix(const std::vector<std::vector<T>>& v){
	m = v.size();
	n = v.size() == 0? 0 : v[0].size();
	mat_.reserve(m*n);
	mat_.insert(mat_.end(), v[0].begin(), v[0].end());

	for(size_t i = 1; i < m; i++){
		mat_.insert(mat_.end(), v[i].begin(), v[i].end());
	}
}

template <typename T> 
matrix<T>::matrix(const matrix& mat){
	//function def here;
	this->m = mat.m;
	this->n = mat.n;
	this->mat_.reserve(m*n);
	this->mat_.insert(this->mat_.end(), mat.mat_.begin(), mat.mat_.end());
}


template <typename T>
matrix<T> & matrix<T>::operator = (const matrix<T>& mat){
	this->m = mat.m;
	this->n = mat.n;
	this->mat_= mat.mat_;
	return *this;
}

template <typename T>
matrix<T> matrix<T>::operator + (const matrix<T> & mat) const{
	matrix<T> outputMat(*this);
	std::transform(outputMat.mat_.begin(),outputMat.mat_.end(), mat.mat_.begin(), outputMat.mat_.begin(), std::plus<T>());
	return outputMat;
}

template <typename T>
matrix<T> & matrix<T>::operator += (const matrix<T>& mat){
	std::transform(this->mat_.begin(),this->mat_.begin(), mat.mat_.begin(), this->mat_.begin(), std::plus<T>());
	return *this;
}

template <typename T>
matrix<T> matrix<T>::operator - (const matrix<T> & mat) const{
	matrix<T> outputMat(*this);
	std::transform(outputMat.mat_.begin(),outputMat.mat_.end(), mat.mat_.begin(), outputMat.mat_.begin(), std::minus<T>());
	return outputMat;
}

template <typename T>
matrix<T> & matrix<T>::operator -= (const matrix<T> & mat){
	std::transform(this->mat_.begin(),this->mat_.begin(), mat.mat_.begin(), this->mat_.begin(), std::plus<T>());
	return *this;
}

template <typename T>
matrix<T> & matrix<T>::operator *= (T scalar){
	std::transform(this->mat_.begin(),this->mat_.end(), this->mat_.begin(), [scalar](T& c){return c * scalar;});
	return *this;
}

template <typename T>
matrix<T> & matrix<T>::operator /= (T scalar){
	std::transform(this->mat_.begin(),this->mat_.end(), this->mat_.begin(), [scalar](T& c){return c / scalar;});
	return *this;
}


template <typename T>
T & matrix<T>::getElement(size_t i, size_t j) const{
	return mat_[i*n + j];
}

template <typename T>
void matrix<T>::setElement(size_t i, size_t j, T val){
	mat_[i*n + j] = val;
}

template <typename T>
matrix<T> matrix<T>::operator * (const matrix<T>& mat) const{
	//check for dimension mismatch... 

	//get correct dimensions
	auto newM = this->m; 
	auto newN = mat.n;

	matrix<T> resMat(newM, newN); // make new matrix; 

	for(size_t i = 0; i < newM; i++){
		auto matStart = resMat.mat_.begin() + i * newN;
		auto matEnd = resMat.mat_.begin() + (i+1) * newN;

		for(size_t j = 0; j < this->n; j++){
			auto targetStart = mat.mat_.begin() + j * newN; 
			std::transform(matStart,matEnd,targetStart,matStart, [this,i,j](T& c1, const T& c2){return c1 + c2 * this->mat_[i*this->n + j];});
		}

	}

	return resMat;
}



















