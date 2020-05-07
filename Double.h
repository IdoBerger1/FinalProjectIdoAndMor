#pragma once
#ifndef  Double_Class
#define Double_Class

#include <immintrin.h>
#include <iostream>

class Double // Type of scalar field
{
	double _num;

public:
	class vec;

	// constructors
	inline Double() : _num(0) {}
	inline Double(double& num) : _num(num) {}
	inline Double(double&& num) : _num(num) {}
	inline Double(Double& F) : _num(F.data()) {}
	inline Double(Double&& F) noexcept : _num(std::move(F.data())) {}

	/* Assignment Operators - START */
	inline Double& operator=(Double& F)
	{
		if (this != &F)
			_num = F.data();
		return *this;
	}

	inline Double& operator=(Double&& F) noexcept
	{
		if (this != &F)
			_num = std::move(F.data());
		return *this;
	}

	inline Double& operator=(double& num)
	{
		if (_num != num)
			_num = num;
		return *this;
	}

	inline Double& operator=(double&& num)
	{
		if (_num != num)
			_num = std::move(num);
		return *this;
	}

	inline Double& operator = (vec& V)
	{
		_mm256_storeu_pd(&_num, V.data());
		return *this;
	}

	inline Double& operator = (vec&& V)
	{
		_mm256_storeu_pd(&_num, V.data());
		return *this;
	}
	/* Assignment Operators - END */

	// Naive Sum
	inline friend Double operator + (Double& A, Double& B) { return Double(A.data() + B.data()); }
	inline friend Double operator + (Double& A, Double&& B) { return Double(A.data() + B.data()); }
	inline friend Double operator + (Double&& A, Double& B) { return Double(A.data() + B.data()); }
	inline friend Double operator + (Double&& A, Double&& B) { return Double(A.data() + B.data()); }

	// Naive Sub
	inline friend Double operator - (Double& A, Double& B) { return Double(A.data() - B.data()); }
	inline friend Double operator - (Double& A, Double&& B) { return Double(A.data() - B.data()); }
	inline friend Double operator - (Double&& A, Double& B) { return Double(A.data() - B.data()); }
	inline friend Double operator - (Double&& A, Double&& B) { return Double(A.data() - B.data()); }

	// Naive Multiplication
	inline friend Double operator * (Double& A, Double& B) { return Double(A.data() * B.data()); }
	inline friend Double operator * (Double& A, Double&& B) { return Double(A.data() * B.data()); }
	inline friend Double operator * (Double&& A, Double& B) { return Double(A.data() * B.data()); }
	inline friend Double operator * (Double&& A, Double&& B) { return Double(A.data() * B.data()); }

	// Naive Division
	inline friend Double operator / (Double& A, Double& B) { return Double(A.data() / B.data()); }
	inline friend Double operator / (Double& A, Double&& B) { return Double(A.data() / B.data()); }
	inline friend Double operator / (Double&& A, Double& B) { return Double(A.data() / B.data()); }
	inline friend Double operator / (Double&& A, Double&& B) { return Double(A.data() / B.data()); }

	// Getters
	inline double data() { return _num; }
	inline double* adress() { return &_num; }

	/* Double class - END */

	/* vec class - START */

	// output operator (most for debug purpose)
	inline friend std::ostream& operator << (std::ostream& out, Double& F)
	{
		out << F.data();
		return out;
	}

	inline friend std::ostream& operator << (std::ostream& out, Double&& F)
	{
		out << F.data();
		return out;
	}

	class vec // Type of AVX vector
	{
		__m256d _v;

	public:
		// constructors
		inline vec() :_v(_mm256_setzero_pd()) {}
		inline vec(__m256d& v) : _v(v) {}
		inline vec(__m256d&& v) : _v(std::move(v)) {}
		inline vec(double* p) : _v(_mm256_loadu_pd(p)) {}
		inline vec(Double* p) : _v(_mm256_loadu_pd(p->adress())) {}
		inline vec(Double& F) : _v(_mm256_broadcast_sd(F.adress())) {}
		inline vec(Double&& F) : _v(_mm256_broadcast_sd(F.adress())) {}
		inline vec(vec& V) : _v(V.data()) {}
		inline vec(vec&& V) noexcept : _v(std::move(V.data())) {}

		// assignment
		inline vec& operator = (vec& V)
		{
			if (this != &V)
				_v = V.data();
			return *this;
		}

		inline vec& operator = (vec&& V) noexcept
		{
			if (this != &V)
				_v = V.data();
			return *this;
		}
		inline vec& operator = (double* p)
		{
			_v = _mm256_loadu_pd(p);
			return *this;
		}

		inline vec& operator = (Double* p)
		{
			_v = _mm256_loadu_pd(p->adress());
			return *this;
		}

		inline vec& operator = (Double& F)
		{
			_v = _mm256_broadcast_sd(F.adress());
			return *this;
		}

		inline vec& operator = (Double&& F)
		{
			_v = _mm256_broadcast_sd(F.adress());
			return *this;
		}

		// geters and seters
		inline __m256d data() { return _v; }
		inline __m256d* adress() { return &_v; }

		// sum operator
		inline friend vec operator+(Double::vec& A, Double::vec& B) { return _mm256_add_pd(A.data(), B.data()); }
		inline friend vec operator+(Double::vec& A, Double::vec&& B) { return _mm256_add_pd(A.data(), B.data()); }
		inline friend vec operator+(Double::vec&& A, Double::vec& B) { return _mm256_add_pd(A.data(), B.data()); }
		inline friend vec operator+(Double::vec&& A, Double::vec&& B) { return _mm256_add_pd(A.data(), B.data()); }

		// sub operator
		inline friend vec operator-(Double::vec& A, Double::vec& B) { return _mm256_sub_pd(A.data(), B.data()); }
		inline friend vec operator-(Double::vec& A, Double::vec&& B) { return _mm256_sub_pd(A.data(), B.data()); }
		inline friend vec operator-(Double::vec&& A, Double::vec& B) { return _mm256_sub_pd(A.data(), B.data()); }
		inline friend vec operator-(Double::vec&& A, Double::vec&& B) { return _mm256_sub_pd(A.data(), B.data()); }

		// multiplication operator
		inline friend vec operator*(Double::vec& A, Double::vec& B) { return _mm256_mul_pd(A.data(), B.data()); }
		inline friend vec operator*(Double::vec& A, Double::vec&& B) { return _mm256_mul_pd(A.data(), B.data()); }
		inline friend vec operator*(Double::vec&& A, Double::vec& B) { return _mm256_mul_pd(A.data(), B.data()); }
		inline friend vec operator*(Double::vec&& A, Double::vec&& B) { return _mm256_mul_pd(A.data(), B.data()); }

		// division operator
		inline friend vec operator/(Double::vec& A, Double::vec& B) { return _mm256_div_pd(A.data(), B.data()); }
		inline friend vec operator/(Double::vec& A, Double::vec&& B) { return _mm256_div_pd(A.data(), B.data()); }
		inline friend vec operator/(Double::vec&& A, Double::vec& B) { return _mm256_div_pd(A.data(), B.data()); }
		inline friend vec operator/(Double::vec&& A, Double::vec&& B) { return _mm256_div_pd(A.data(), B.data()); }


		// A * B + C
		inline vec mul_add(Double::vec& A, Double::vec& B, Double::vec& C) { return _mm256_fmadd_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_add(Double::vec& A, Double::vec& B, Double::vec&& C) { return _mm256_fmadd_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_add(Double::vec& A, Double::vec&& B, Double::vec& C) { return _mm256_fmadd_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_add(Double::vec& A, Double::vec&& B, Double::vec&& C) { return _mm256_fmadd_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_add(Double::vec&& A, Double::vec& B, Double::vec& C) { return _mm256_fmadd_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_add(Double::vec&& A, Double::vec& B, Double::vec&& C) { return _mm256_fmadd_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_add(Double::vec&& A, Double::vec&& B, Double::vec& C) { return _mm256_fmadd_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_add(Double::vec&& A, Double::vec&& B, Double::vec&& C) { return _mm256_fmadd_pd((A.data()), (B.data()), (C.data())); }


		// A * B - C
		inline vec mul_sub(Double::vec& A, Double::vec& B, Double::vec& C) { return _mm256_fmsub_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_sub(Double::vec& A, Double::vec& B, Double::vec&& C) { return _mm256_fmsub_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_sub(Double::vec& A, Double::vec&& B, Double::vec& C) { return _mm256_fmsub_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_sub(Double::vec& A, Double::vec&& B, Double::vec&& C) { return _mm256_fmsub_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_sub(Double::vec&& A, Double::vec& B, Double::vec& C) { return _mm256_fmsub_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_sub(Double::vec&& A, Double::vec& B, Double::vec&& C) { return _mm256_fmsub_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_sub(Double::vec&& A, Double::vec&& B, Double::vec& C) { return _mm256_fmsub_pd((A.data()), (B.data()), (C.data())); }
		inline vec mul_sub(Double::vec&& A, Double::vec&& B, Double::vec&& C) { return _mm256_fmsub_pd((A.data()), (B.data()), (C.data())); }
	};
};
#endif // ! Double_Class