#ifndef Vector_Class
#define Vector_Class

#include "Matrix.h"

template <typename scalar>
class Vector : public Matrix<scalar>
{
public:
	using Matrix<scalar>::Matrix;

	// constructors
	Vector() : Matrix<scalar>::Matrix() {} // empty constructor
	Vector(size_t len) : Matrix<scalar>::Matrix(len, 1) {} // Vector of size len
	Vector(size_t len, scalar val) : Matrix<scalar>::Matrix(len, 1, val) {} // Vector filled by val
	Vector(size_t len, string type) : Matrix<scalar>::Matrix(len, 1, type) {} //rand Vector ...
	Vector(initializer_list<scalar> list) 
	{
		initializer_list<initializer_list<scalar>> listOfLists = { list };
		Matrix<scalar>::Matrix(listOfLists);
	}
	Vector(vector<scalar>& vec)
	{
		vector<vector<scalar>> vectorOfVectors(vec);
		Matrix<scalar>::Matrix(vectorOfVectors);
	}
	Vector(vector<scalar>&& vec)
	{
		vector<vector<scalar>> vectorOfVectors(vec);
		Matrix<scalar>::Matrix(vectorOfVectors);
	}
	Vector(const Vector& V) : Matrix<scalar>::Matrix(V) {} // lvalue copy constructor
	Vector(Vector& V) : Matrix<scalar>::Matrix(V) {} // lvalue copy constructor
	Vector(Vector&& V) : Matrix<scalar>::Matrix(V) {} // rvalue copy constructor
};
#endif // !Vector_Class