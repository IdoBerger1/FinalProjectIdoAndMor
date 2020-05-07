#pragma once
#ifndef Int_Class
#define Int_Class

#include <immintrin.h>
#include <iostream>

class Int // Type of scalar field
{
	int _num;

public:
	class vec;

	// constructors
	inline Int() : _num(0) {}
	inline Int(int& num) : _num(num) {}
	inline Int(int&& num) : _num(num) {}
	inline Int(Int& I) : _num(I._num) {}
	inline Int(Int&& I) noexcept : _num(move(I._num)) {}

	/* Assignment Operators - START */
	inline Int& operator=(Int& I)
	{
		if (this != &I)
			_num = I._num;
		return *this;
	}

	inline Int& operator=(Int&& I) noexcept
	{
		if (this != &I)
			_num = move(I._num);
		return *this;
	}

	inline Int& operator=(int& num)
	{
		if (_num != num)
			_num = num;
		return *this;
	}

	inline Int& operator=(int&& num)
	{
		if (_num != num)
			_num = move(num);
		return *this;
	}

	inline Int& operator = (vec& V)
	{
		_mm256_storeu_si256((__m256i*) & _num, V.data());
		return *this;
	}

	inline Int& operator = (vec&& V)
	{
		_mm256_storeu_si256((__m256i*) & _num, V.data());
		return *this;
	}
	/* Assignment Operators - END */

	// Naive Sum
	inline friend Int operator + (Int& A, Int& B) { return Int(A.data() + B.data()); }
	inline friend Int operator + (Int& A, Int&& B) { return Int(A.data() + B.data()); }
	inline friend Int operator + (Int&& A, Int& B) { return Int(A.data() + B.data()); }
	inline friend Int operator + (Int&& A, Int&& B) { return Int(A.data() + B.data()); }

	// Naive Sub
	inline friend Int operator - (Int& A, Int& B) { return Int(A.data() - B.data()); }
	inline friend Int operator - (Int& A, Int&& B) { return Int(A.data() - B.data()); }
	inline friend Int operator - (Int&& A, Int& B) { return Int(A.data() - B.data()); }
	inline friend Int operator - (Int&& A, Int&& B) { return Int(A.data() - B.data()); }

	// Naive Multiplication
	inline friend Int operator * (Int& A, Int& B) { return Int(A.data() * B.data()); }
	inline friend Int operator * (Int& A, Int&& B) { return Int(A.data() * B.data()); }
	inline friend Int operator * (Int&& A, Int& B) { return Int(A.data() * B.data()); }
	inline friend Int operator * (Int&& A, Int&& B) { return Int(A.data() * B.data()); }

	// Naive Division
	inline friend Int operator / (Int& A, Int& B) { return Int(A.data() / B.data()); }
	inline friend Int operator / (Int& A, Int&& B) { return Int(A.data() / B.data()); }
	inline friend Int operator / (Int&& A, Int& B) { return Int(A.data() / B.data()); }
	inline friend Int operator / (Int&& A, Int&& B) { return Int(A.data() / B.data()); }

	// Getters
	inline int data() { return _num; }
	inline int* adress() { return &_num; }

	/* Double class - END */

	/* vec class - START */

	// output operator (most for debug purpose)
	inline friend std::ostream& operator << (std::ostream& out, Int& I)
	{
		out << I.data();
		return out;
	}

	inline friend std::ostream& operator << (std::ostream& out, Int&& I)
	{
		out << I.data();
		return out;
	}

	class vec // Type of AVX vector
	{
		__m256i _v;

	public:
		// constructors
		inline vec() :_v(_mm256_setzero_si256()) {}
		inline vec(__m256i& v) : _v(v) {}
		inline vec(__m256i&& v) : _v(move(v)) {}
		inline vec(int* p) : _v(_mm256_loadu_si256((__m256i*)p)) {}
		inline vec(Int* p) : _v(_mm256_loadu_si256((__m256i*)p->adress())) {}
		inline vec(Int& I) : _v(_mm256_set1_epi32(I.data())) {} // broadcast
		inline vec(Int&& I) : _v(_mm256_set1_epi32(I.data())) {} // broadcast
		inline vec(vec& V) : _v(V.data()) {}
		inline vec(vec&& V) noexcept : _v(move(V.data())) {}

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
		inline vec& operator = (int* p)
		{
			_v = _mm256_loadu_si256((__m256i*)p);
			return *this;
		}

		inline vec& operator = (Int* p)
		{
			_v = _mm256_loadu_si256((__m256i*)p->adress());
			return *this;
		}

		inline vec& operator = (Int& I)
		{
			_v = _mm256_set1_epi32(I.data()); // broadcast
			return *this;
		}

		inline vec& operator = (Int&& I)
		{
			_v = _mm256_set1_epi32(I.data()); // broadcast
			return *this;
		}

		// geters and seters
		inline __m256i data() { return _v; }
		inline __m256i* adress() { return &_v; }

		// sum operator
		inline friend vec operator + (Int::vec& A, Int::vec& B) { return vec(_mm256_add_epi32(A.data(), B.data())); }
		inline friend vec operator + (Int::vec& A, Int::vec&& B) { return vec(_mm256_add_epi32(A.data(), B.data())); }
		inline friend vec operator + (Int::vec&& A, Int::vec& B) { return vec(_mm256_add_epi32(A.data(), B.data())); }
		inline friend vec operator + (Int::vec&& A, Int::vec&& B) { return vec(_mm256_add_epi32(A.data(), B.data())); }

		// sub operator
		inline friend vec operator - (Int::vec& A, Int::vec& B) { return vec(_mm256_sub_epi32(A.data(), B.data())); }
		inline friend vec operator - (Int::vec& A, Int::vec&& B) { return vec(_mm256_sub_epi32(A.data(), B.data())); }
		inline friend vec operator - (Int::vec&& A, Int::vec& B) { return vec(_mm256_sub_epi32(A.data(), B.data())); }
		inline friend vec operator - (Int::vec&& A, Int::vec&& B) { return vec(_mm256_sub_epi32(A.data(), B.data())); }

		// multiplication operator
		inline friend vec operator * (Int::vec& A, Int::vec& B) { return vec(_mm256_mullo_epi32(A.data(), B.data())); }
		inline friend vec operator * (Int::vec& A, Int::vec&& B) { return vec(_mm256_mullo_epi32(A.data(), B.data())); }
		inline friend vec operator * (Int::vec&& A, Int::vec& B) { return vec(_mm256_mullo_epi32(A.data(), B.data())); }
		inline friend vec operator * (Int::vec&& A, Int::vec&& B) { return vec(_mm256_mullo_epi32(A.data(), B.data())); }

		// division operator
#ifdef _WIN32 || _WIN64 
		inline friend vec operator / (Int::vec& A, Int::vec& B) { return vec(_mm256_div_epi32(A.data(), B.data())); }
		inline friend vec operator / (Int::vec& A, Int::vec&& B) { return vec(_mm256_div_epi32(A.data(), B.data())); }
		inline friend vec operator / (Int::vec&& A, Int::vec& B) { return vec(_mm256_div_epi32(A.data(), B.data())); }
		inline friend vec operator / (Int::vec&& A, Int::vec&& B) { return vec(_mm256_div_epi32(A.data(), B.data())); }

#elif __linux__ || __unix || __unix__
		inline friend vec operator / (vec& A, vec& B) { return vec(_mm256_castps_si256(_mm256_div_ps(_mm256_castsi256_ps(A.data()), _mm256_castsi256_ps(B.data())))); }
		inline friend vec operator / (vec& A, vec&& B) { return vec(_mm256_castps_si256(_mm256_div_ps(_mm256_castsi256_ps(A.data()), _mm256_castsi256_ps(B.data())))); }
		inline friend vec operator / (vec&& A, vec& B) { return vec(_mm256_castps_si256(_mm256_div_ps(_mm256_castsi256_ps(A.data()), _mm256_castsi256_ps(B.data())))); }
		inline friend vec operator / (vec&& A, vec&& B) { return vec(_mm256_castps_si256(_mm256_div_ps(_mm256_castsi256_ps(A.data()), _mm256_castsi256_ps(B.data())))); }
#endif

		// A * B + C
		inline vec mul_add(Int::vec& A, Int::vec& B, Int::vec& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_add_epi32(product, C.data()));
		}
		inline vec mul_add(Int::vec& A, Int::vec& B, Int::vec&& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_add_epi32(product, C.data()));
		}
		inline vec mul_add(Int::vec& A, Int::vec&& B, Int::vec& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_add_epi32(product, C.data()));
		}
		inline vec mul_add(Int::vec& A, Int::vec&& B, Int::vec&& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_add_epi32(product, C.data()));
		}
		inline vec mul_add(Int::vec&& A, Int::vec& B, Int::vec& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_add_epi32(product, C.data()));
		}
		inline vec mul_add(Int::vec&& A, Int::vec& B, Int::vec&& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_add_epi32(product, C.data()));
		}
		inline vec mul_add(Int::vec&& A, Int::vec&& B, Int::vec& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_add_epi32(product, C.data()));
		}
		inline vec mul_add(Int::vec&& A, Int::vec&& B, Int::vec&& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_add_epi32(product, C.data()));
		}


		// A * B - C
		inline vec mul_sub(Int::vec& A, Int::vec& B, Int::vec& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_sub_epi32(product, C.data()));
		}
		inline vec mul_sub(Int::vec& A, Int::vec& B, Int::vec&& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_sub_epi32(product, C.data()));
		}
		inline vec mul_sub(Int::vec& A, Int::vec&& B, Int::vec& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_sub_epi32(product, C.data()));
		}
		inline vec mul_sub(Int::vec& A, Int::vec&& B, Int::vec&& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_sub_epi32(product, C.data()));
		}
		inline vec mul_sub(Int::vec&& A, Int::vec& B, Int::vec& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_sub_epi32(product, C.data()));
		}
		inline vec mul_sub(Int::vec&& A, Int::vec& B, Int::vec&& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_sub_epi32(product, C.data()));
		}
		inline vec mul_sub(Int::vec&& A, Int::vec&& B, Int::vec& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_sub_epi32(product, C.data()));
		}
		inline vec mul_sub(Int::vec&& A, Int::vec&& B, Int::vec&& C)
		{
			__m256i product = _mm256_mullo_epi32(A.data(), B.data());
			return  vec(_mm256_sub_epi32(product, C.data()));
		}
	};
};
#endif // ! INT_Class