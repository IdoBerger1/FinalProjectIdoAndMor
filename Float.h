#ifndef  Float_Class
#define Float_Class

#include <immintrin.h>
#include <iostream>
using namespace std;

class Float // Type of scalar field
{
	float _num;

public:
	class vec;

	// constructors
	inline Float() : _num(0) {}
	inline Float(float& num) : _num(num) {}
	inline Float(float&& num) : _num(move(num)) {}
	inline Float(Float& F) : _num(F.data()) {}
	inline Float(Float&& F) noexcept : _num(move(F.data())) {}

	// assignment
	inline Float& operator = (Float& F)
	{
		if (this != &F)
			_num = F.data();
		return *this;
	}

	inline Float& operator = (Float&& F) noexcept
	{
		if (this != &F)
			_num = move(F.data());
		return *this;
	}

	inline Float& operator = (float& num)
	{
		if (&_num != &num)
			_num = num;
		return *this;
	}

	inline Float& operator = (float&& num)
	{
		if (&_num != &num)
			_num = num;
		return *this;
	}

	inline Float& operator = (vec& V)
	{
		_mm256_storeu_ps(&_num, V.data());
		return *this;
	}

	inline Float& operator = (vec&& V)
	{
		_mm256_storeu_ps(&_num, V.data());
		return *this;
	}

	// naive sum operator
	inline friend Float operator + (Float& A, Float& B) { return Float(A.data() + B.data()); }
	inline friend Float operator + (Float& A, Float&& B) { return Float(A.data() + B.data()); }
	inline friend Float operator + (Float&& A, Float& B) { return Float(A.data() + B.data()); }
	inline friend Float operator + (Float&& A, Float&& B) { return Float(A.data() + B.data()); }

	// naive sub operator
	inline friend Float operator - (Float& A, Float& B) { return Float(A.data() - B.data()); }
	inline friend Float operator - (Float& A, Float&& B) { return Float(A.data() - B.data()); }
	inline friend Float operator - (Float&& A, Float& B) { return Float(A.data() - B.data()); }
	inline friend Float operator - (Float&& A, Float&& B) { return Float(A.data() - B.data()); }

	// naive multiplication operator
	inline friend Float operator * (Float& A, Float& B) { return Float(A.data() * B.data()); }
	inline friend Float operator * (Float& A, Float&& B) { return Float(A.data() * B.data()); }
	inline friend Float operator * (Float&& A, Float& B) { return Float(A.data() * B.data()); }
	inline friend Float operator * (Float&& A, Float&& B) { return Float(A.data() * B.data()); }

	// naive division operator
	inline friend Float operator / (Float& A, Float& B) { return Float(A.data() / B.data()); }
	inline friend Float operator / (Float& A, Float&& B) { return Float(A.data() / B.data()); }
	inline friend Float operator / (Float&& A, Float& B) { return Float(A.data() / B.data()); }
	inline friend Float operator / (Float&& A, Float&& B) { return Float(A.data() / B.data()); }

	// accessors
	inline float data() { return _num; }
	inline float* adress() { return &_num; }

	// output operator (most for debug purpose)
	inline friend ostream& operator << (ostream& out, Float& F) 
	{
		out << F.data();
		return out;
	}

	inline friend ostream& operator << (ostream& out, Float&& F)
	{
		out << F.data();
		return out;
	}

	class vec // Type of AVX vector
	{
		__m256 _v;

	public:
		// constructors
		inline vec() : _v(_mm256_setzero_ps()) {}
		inline vec(__m256& v) : _v(v) {}
		inline vec(__m256&& v) : _v(move(v)) {}
		inline vec(float* p) : _v(_mm256_loadu_ps(p)) {}
		inline vec(Float* p) : _v(_mm256_loadu_ps(p->adress())) {}
		inline vec(Float& F) : _v(_mm256_broadcast_ss(F.adress())) {}
		inline vec(Float&& F) : _v(move(_mm256_broadcast_ss(F.adress()))) {}
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
		inline vec& operator = (float* p)
		{
			_v = _mm256_loadu_ps(p);
			return *this;
		}

		inline vec& operator = (Float* p)
		{
			_v = _mm256_loadu_ps(p->adress());
			return *this;
		}

		inline vec& operator = (Float& F)
		{
			_v = _mm256_broadcast_ss(F.adress());
			return *this;
		}

		inline vec& operator = (Float&& F)
		{
			_v = _mm256_broadcast_ss(F.adress());
			return *this;
		}

		// geters and seters
		inline __m256 data() { return _v; }
		inline __m256* adress() { return &_v; }

		// sum operator
		inline friend vec operator + (vec& A, vec& B) { return vec(_mm256_add_ps(A.data(), B.data())); }
		inline friend vec operator + (vec& A, vec&& B) { return vec(_mm256_add_ps(A.data(), B.data())); }
		inline friend vec operator + (vec&& A, vec& B) { return vec(_mm256_add_ps(A.data(), B.data())); }
		inline friend vec operator + (vec&& A, vec&& B) { return vec(_mm256_add_ps(A.data(), B.data())); }

		// sub operator
		inline friend vec operator - (vec& A, vec& B) { return vec(_mm256_sub_ps(A.data(), B.data())); }
		inline friend vec operator - (vec& A, vec&& B) { return vec(_mm256_sub_ps(A.data(), B.data())); }
		inline friend vec operator - (vec&& A, vec& B) { return vec(_mm256_sub_ps(A.data(), B.data())); }
		inline friend vec operator - (vec&& A, vec&& B) { return vec(_mm256_sub_ps(A.data(), B.data())); }

		// multiplication operator
		inline friend vec operator * (vec& A, vec& B) { return vec(_mm256_mul_ps(A.data(), B.data())); }
		inline friend vec operator * (vec& A, vec&& B) { return vec(_mm256_mul_ps(A.data(), B.data())); }
		inline friend vec operator * (vec&& A, vec& B) { return vec(_mm256_mul_ps(A.data(), B.data())); }
		inline friend vec operator * (vec&& A, vec&& B) { return vec(_mm256_mul_ps(A.data(), B.data())); }

		// division operator
		inline friend vec operator / (vec& A, vec& B) { return vec(_mm256_div_ps(A.data(), B.data())); }
		inline friend vec operator / (vec& A, vec&& B) { return vec(_mm256_div_ps(A.data(), B.data())); }
		inline friend vec operator / (vec&& A, vec& B) { return vec(_mm256_div_ps(A.data(), B.data())); }
		inline friend vec operator / (vec&& A, vec&& B) { return vec(_mm256_div_ps(A.data(), B.data())); }

		// A * B + C
		inline vec mul_add(vec& A, vec& B, vec& C) { return vec(_mm256_fmadd_ps(A.data(), B.data(), C.data())); }
		inline vec mul_add(vec&& A, vec& B, vec& C) { return vec(_mm256_fmadd_ps(A.data(), B.data(), C.data())); }
		inline vec mul_add(vec& A, vec&& B, vec& C) { return vec(_mm256_fmadd_ps(A.data(), B.data(), C.data())); }
		inline vec mul_add(vec& A, vec& B, vec&& C) { return vec(_mm256_fmadd_ps(A.data(), B.data(), C.data())); }
		inline vec mul_add(vec&& A, vec&& B, vec& C) { return vec(_mm256_fmadd_ps(A.data(), B.data(), C.data())); }
		inline vec mul_add(vec&& A, vec& B, vec&& C) { return vec(_mm256_fmadd_ps(A.data(), B.data(), C.data())); }
		inline vec mul_add(vec& A, vec&& B, vec&& C) { return vec(_mm256_fmadd_ps(A.data(), B.data(), C.data())); }
		inline vec mul_add(vec&& A, vec&& B, vec&& C) { return vec(_mm256_fmadd_ps(A.data(), B.data(), C.data())); }

		// A * B - C
		inline vec mul_sub(vec& A, vec& B, vec& C) { return vec(_mm256_fmsub_ps(A.data(), B.data(), C.data())); }
		inline vec mul_sub(vec&& A, vec& B, vec& C) { return vec(_mm256_fmsub_ps(A.data(), B.data(), C.data())); }
		inline vec mul_sub(vec& A, vec&& B, vec& C) { return vec(_mm256_fmsub_ps(A.data(), B.data(), C.data())); }
		inline vec mul_sub(vec& A, vec& B, vec&& C) { return vec(_mm256_fmsub_ps(A.data(), B.data(), C.data())); }
		inline vec mul_sub(vec&& A, vec&& B, vec& C) { return vec(_mm256_fmsub_ps(A.data(), B.data(), C.data())); }
		inline vec mul_sub(vec&& A, vec& B, vec&& C) { return vec(_mm256_fmsub_ps(A.data(), B.data(), C.data())); }
		inline vec mul_sub(vec& A, vec&& B, vec&& C) { return vec(_mm256_fmsub_ps(A.data(), B.data(), C.data())); }
		inline vec mul_sub(vec&& A, vec&& B, vec&& C) { return vec(_mm256_fmsub_ps(A.data(), B.data(), C.data())); }
	};
};
#endif // ! Float_Class