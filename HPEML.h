#ifndef HPEML
#define HPEML // stands for: High Performance Extendable Math Library

// Scalars
#include "Float.h"
#include "Double.h"
#include "Int.h"

// Math Operations
#include "Memory_Block.h"
#include "Matrix.h"
#include "Vector.h"

// System Interface
#include "Initialization.h"

using namespace std;
using namespace System;

#define	IN
#define OUT
#define INOUT
#define EPS 10e-08
#define PHI 3.14

typedef struct complex
{
	double a;
	double b;
}complex;

void copyVec(IN Matrix<Double>& mat1, OUT Matrix<Double>& mat2, IN Int startLocation);
void normalize(INOUT Matrix<Double>& mat1, IN Int startLocation);
Double dotProduct(IN Matrix<Double> vecA, IN Matrix<Double> vecB, IN Int vecAIndex, IN Int vecBIndex);
void findQ(IN Matrix<Double>& mat, OUT Matrix<Double>& Q, IN Int n);
void subVec(double* vecA, double* vecB, int n);
double* scalarMult(IN Double scalar, IN double* Q, IN int n);
void findR(IN Matrix<Double>& mat, IN Matrix<Double>& Q, OUT Matrix<Double>& R, IN Int n);
complex* qrAlgorithm(IN Matrix<Double>& mat, IN Int n, IN Int k);
void multiplyMatrix(IN Matrix<Double>& matrix1, IN Matrix<Double>& matrix2, OUT Matrix<Double>& calculatedMatrix1);
void JacobiAlgorithm(IN Matrix<Double>& matrixA, IN Int n, OUT Matrix<Double>& calculatedMatrix, OUT double* eigenValues);
void svdAlgorithm(IN Matrix<Double>& matrixA, IN Int n, OUT Matrix<Double>& U, OUT Matrix<Double>& sig, OUT Matrix<Double>& V);

#endif