#include "HPEML.h"

complex* qrAlgorithm(IN Matrix<Double>& mat, IN Int n, IN Int k)
{
	Int i(0);
	Matrix<Double> Q(n.data(), n.data());
	Matrix<Double> R(n.data(), n.data());
	Matrix<Double> A(n.data(), n.data());
	Matrix<Double> matrixEigenVectors(n.data(), n.data(), "Identity");
	Matrix<Double> tempMatrix(n.data(), n.data());
	complex* eigenValues = (complex*)malloc(n.data() * sizeof(complex));

	A = mat;

	for (i = 1; i.data() < k.data(); i = i + 1)
	{
		findQ(A, Q, n);
		tempMatrix = matrixEigenVectors * Q;
		matrixEigenVectors = tempMatrix;
		findR(A, Q, R, n);
		A = R * Q;
	}

	for (i = 0; i.data() < n.data(); i = i + 1)
	{
		// complex numbers (a+bi) and (a-bi)
		// a = A[i][i], b = A[i][n-i-1]
		// a = A[n-i-1][n-i-1], -b = A[n-i-1][i]
		if (abs(A(i.data(), n.data() - i.data() - 1).data()) > EPS)
		{
			eigenValues[i.data()].a = (double)A(i.data(), i.data()).data();
			eigenValues[i.data()].b = (double)A(i.data(), n.data() - i.data() - 1).data();
		}
		// Real numbers
		else
		{
			eigenValues[i.data()].a = (double)A(i.data(), i.data()).data();
			eigenValues[i.data()].b = 0;
		}
	}
	return eigenValues;
}

// subVec - right vector - left vector
// scalarMult - mult vector in scalar
// dotProduct - mult two vectors
void findQ(IN Matrix<Double>& mat, OUT Matrix<Double>& Q, IN Int n)
{
	Int i(0), j(0);
	//Matrix<Int> tempVector(1, n.data(), 0);
	double* tempVector = NULL;
	mat.trans(true);
	copyVec(Q, mat, 0);
	normalize(Q, 0);
	for (i = 1; i.data() < n.data(); i = i + 1)
	{
		copyVec(Q, mat, i);

		for (j = 0; j.data() < i.data(); j = j + 1)
		{
			tempVector = scalarMult(dotProduct(mat, Q, i, j), (double*)Q[j.data()], Q.rows());
			subVec((double*)Q[i.data()], tempVector, Q.rows());
			free(tempVector);
		}

		normalize(Q, i);
	}

	Q.trans(true);
}

// multiplyMat - allocate nxn memory
// return from the function: pointer to mem allocation.
void findR(IN Matrix<Double>& mat, IN Matrix<Double>& Q, OUT Matrix<Double>& R, IN Int n)
{
	Q.trans(true);
	R = Q * mat;
	Q.trans(true);
}

void copyVec(OUT Matrix<Double>& mat1, IN Matrix<Double>& mat2, IN Int startLocation)
{
	Int i(0);
	for (i = 0; i.data() < mat2.rows(); i = i + 1)
	{
		mat1(startLocation.data(), i.data()) = mat2(startLocation.data(), i.data());
	}
}

void normalize(INOUT Matrix<Double>& mat1, IN Int startLocation)
{
	Int i(0);
	Double sum1(0);

	for (i = 0; i.data() < mat1.rows(); i = i + 1)
	{
		sum1 = sum1 + mat1(startLocation.data(), i.data()) * mat1(startLocation.data(), i.data());
	}

	sum1 = sqrt(sum1.data());

	for (i = 0; i.data() < mat1.rows(); i = i + 1)
	{
		mat1(startLocation.data(), i.data()) = mat1(startLocation.data(), i.data()) / sum1.data();
	}
}

// dotProduct - mult two vectors
Double dotProduct(IN Matrix<Double> vecA, IN Matrix<Double> vecB, IN Int vecAIndex, IN Int vecBIndex)
{
	Int i(0);
	Double sum(0);

	for (i = 0; i.data() < vecA.rows(); i = i + 1)
	{
		sum = sum + vecA(vecAIndex.data(), i.data()) * vecB(vecBIndex.data(), i.data());
	}
	return sum;
}

double* scalarMult(IN Double scalar, IN double* Q, IN int n)
{
	Int i(0);

	double* vector = (double*)malloc(n * sizeof(double));

	for (i = 0; i.data() < n; i = i + 1)
	{
		vector[i.data()] = scalar.data() * Q[i.data()];
	}
	return vector;
}

// subVec - right vector - left vector
void subVec(double* vecA, double* vecB, int n)
{
	int i;

	for (i = 0; i < n; i++)
	{
		vecA[i] = vecA[i] - vecB[i];
	}
}

void multiplyMatrix(IN Matrix<Double>& matrix1, IN Matrix<Double>& matrix2, OUT Matrix<Double>& calculatedMatrix1)
{
	Int i(0), j(0), k(0);

	for (i = 0; i.data() < matrix1.rows(); i = i + 1)
		for (j = 0; j.data() < matrix1.cols(); j = j + 1)
			calculatedMatrix1(i.data(), j.data()) = 0;

	for (i = 0; i.data() < matrix1.rows(); i = i + 1)
	{
		for (j = 0; j.data() < matrix1.rows(); j = j + 1)
		{
			for (k = 0; k.data() < matrix1.rows(); k = k + 1)
			{
				calculatedMatrix1(i.data(), j.data()) = calculatedMatrix1(i.data(), j.data()) + matrix1(i.data(), k.data()) * matrix2(k.data(), j.data());
			}
		}
	}
}

void JacobiAlgorithm(IN Matrix<Double>& matrixA, IN Int n, OUT Matrix<Double>& calculatedMatrix, OUT double* eigenValues)
{
	Matrix<Double> D(n.data(), n.data());
	Matrix<Double> S(n.data(), n.data(), "Identity");
	Matrix<Double> S1(n.data(), n.data(), "Identity");
	Int i(0), j(0), indexI(0), indexJ(0);
	double tetha = 0, maxValue = 0;

	D = matrixA;

	// Step2
	do
	{
		// Step3: Find the largest off diagonal value
		for (i = 0; i.data() < n.data(); i = i + 1)
		{
			for (j = i + 1;j.data() < n.data(); j = j + 1)
			{
				if (D(i.data(), j.data()).data() > maxValue)
				{
					maxValue = D(i.data(), j.data()).data();
					indexI = i;
					indexJ = j;
				}
			}
		}

		// Step4: Find the rotation angle
		if (D(indexI.data(), indexI.data()).data() == D(indexJ.data(), indexJ.data()).data())
		{
			if (D(indexI.data(), indexJ.data()).data() > 0)
			{
				tetha = PHI / 4;
			}
			else
			{
				tetha = -PHI / 4;
			}
		}
		else
		{
			tetha = 0.5 * atan((2 * D(indexI.data(), indexJ.data()).data()) / (D(indexI.data(), indexI.data()).data() - D(indexJ.data(), indexJ.data()).data()));
		}

		// Step5: Compute the matrix S1
		S1(indexI.data(), indexI.data()) = S1(indexJ.data(), indexJ.data()) = cos(tetha);
		S1(indexI.data(), indexJ.data()) = -sin(tetha);
		S1(indexJ.data(), indexI.data()) = sin(tetha);

		// Step6: Find D and S
		S = S * S1;
		D = S1.trans(false) * D * S1;

	} while (maxValue < EPS /**NORMM_FROBENIUS(S)*/);
	// Step7 ended: D is diagonal.

	// Step8: Diagonal element of D are the eigenValues 
	// The columns of S are the eigenVectors.

	// Find eigenValues
	for (i = 0; i.data() < n.data(); i = i + 1)
	{
		eigenValues[i.data()] = D(i.data(), i.data()).data();
	}

	calculatedMatrix = S;
}

void svdAlgorithm(IN Matrix<Double>& matrixA, IN Int n, OUT Matrix<Double>& U, OUT Matrix<Double>& sig, OUT Matrix<Double>& V)
{
	Matrix<Double> AtA(n.data(), n.data());
	Matrix<Double> AAt(n.data(), n.data());
	Matrix<Double> matrixAt = matrixA.trans(false);
	Matrix<Double> A(n.data(), n.data());

	double* eigenValues = NULL;
	Int i(0), j(0);

	// Zero sigma matrix.
	for (i = 0; i.data() < sig.rows(); i = i + 1)
		for (j = 0; j.data() < sig.cols(); j = j + 1)
			sig(i.data(), j.data()) = 0;

	eigenValues = (double*)malloc(sizeof(double) * n.data());

	AtA = matrixAt * matrixA;
	AAt = matrixA * matrixAt;

	JacobiAlgorithm(AtA, n, V, eigenValues);
	JacobiAlgorithm(AAt, n, U, eigenValues);

	for (i = 0; i.data() < n.data(); i = i + 1)
	{
		sig(i.data(), i.data()) = eigenValues[i.data()];
	}

	V.trans(true);
	A = U * sig * V;
}


