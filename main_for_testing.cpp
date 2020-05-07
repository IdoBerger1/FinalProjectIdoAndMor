#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <x86intrin.h>
#include <mkl.h>
#include <thread>
#include <mutex>
#include <immintrin.h>

#include "t_timer.h"
#include <unistd.h>
using namespace std;

/////////////////////////////////////////Definition Section/////////////////////////////////////////
float *sum(float *A, float *B, size_t row, size_t col);
float *sub(float *A, float *B, size_t row, size_t col);
void sub_matrix_row(float *A, size_t colA, size_t r_1, size_t r_2, size_t c_1, size_t c_2, float *B);
void collect_mat_row(float *A0, float *A1, float *A2, float *A3, float *B, size_t row, size_t col);
void print_matrix(float *A, size_t row, size_t col);
bool matrix_equal(float *A, float *B, size_t row, size_t col);
inline float *Mult_AVX(float *A, float *B, size_t rowA, size_t colA_rowB, size_t colB);
// float *Mult_AVX_multithread_S(float *A, float *B, size_t n, int iter);
float *Mult_AVX_multithread_SW(float *A, float *B, size_t rowA, size_t colA_rowB, size_t colB, int iter);
inline tuple<float *, size_t, size_t> paddingBlockSize(float *matrix, size_t row, size_t col);
inline tuple<float *, size_t, size_t> paddSingleRowAndColumn(float *matrix, size_t row, size_t col);
inline float *removePadding(float *matrix, size_t row, size_t col, size_t colAfterPadding);
void test(size_t rowA, size_t colA_rowB, size_t colB);
///////////////////////////////////////End Definition Section///////////////////////////////////////

/////////////////////////////////////////////Constance//////////////////////////////////////////////
#define VECSIZE 8
#define I_BLOCKSIZE 256
#define J_BLOCKSIZE 256
#define K_BLOCKSIZE 256
#define UNROLL_1 4 // UNROLL_1 * VECSIZE <= min{I_BLOCKSIZE, J_BLOCKSIZE, K_BLOCKSIZE}
#define UNROLL_2 2 //UNROLL_2 * VECSIZE <= min{I_BLOCKSIZE, J_BLOCKSIZE, K_BLOCKSIZE}
#define ITER_NUM 1

mutex mutex1;
//////////////////////////////////////////End Constance/////////////////////////////////////////////

/////////////////////////////////////////Sum of Matrices//////////////////////////////////////////
float *sum(float *A, float *B, size_t row, size_t col)
{
	size_t vecsize = VECSIZE, SizeC = row * col, i;
	float *C = new float[SizeC];

	if (SizeC >= vecsize)
	{
		for (i = 0; i < SizeC - vecsize; i += vecsize)
		{
			_mm256_storeu_ps(C + i, _mm256_add_ps(_mm256_loadu_ps(A + i), _mm256_loadu_ps(B + i)));
		}

		for (; i < SizeC; i++)
		{
			C[i] = A[i] + B[i];
		}
	}
	else
	{
		for (i = 0; i < SizeC; i++)
		{
			C[i] = A[i] + B[i];
		}
	}

	return C;
}
///////////////////////////////////////End Sum of Matrices////////////////////////////////////////

/////////////////////////////////////Substration of Matrices//////////////////////////////////////
float *sub(float *A, float *B, size_t row, size_t col)
{
	size_t vecsize = VECSIZE, SizeC = row * col, i;
	float *C = new float[SizeC];

	if (SizeC >= vecsize)
	{
		for (i = 0; i < SizeC - vecsize; i += vecsize)
		{
			_mm256_storeu_ps(C + i, _mm256_sub_ps(_mm256_loadu_ps(A + i), _mm256_loadu_ps(B + i)));
		}

		for (; i < SizeC; i++)
		{
			C[i] = A[i] - B[i];
		}
	}
	else
	{
		for (i = 0; i < SizeC; i++)
			C[i] = A[i] - B[i];
	}

	return C;
}
////////////////////////////////////End Substration of Matrices////////////////////////////////////

////////////////////////////////////Extraction of sub-matrices/////////////////////////////////////
void sub_matrix_row(float *A, size_t colA, size_t r_1, size_t r_2, size_t c_1, size_t c_2, float *B)
{
	size_t vecsize = VECSIZE, i, j, k;
	if (c_2 >= vecsize)
	{
		for (i = r_1, k = 0; i < r_2; i++)
		{
			for (j = c_1; j < c_2 - vecsize; j += vecsize, k += vecsize)
			{
				_mm256_storeu_ps(B + k, _mm256_loadu_ps(A + i * colA + j));
			}

			for (; j < c_2; j++, k++)
			{
				B[k] = A[i * colA + j];
			}
		}
	}
	else
	{
		for (i = r_1, k = 0; i < r_2; i++)
		{
			for (j = c_1; j < c_2; j++, k++)
			{
				B[k] = A[i * colA + j];
			}
		}
	}
}
//////////////////////////////////End Extraction of sub-matrices///////////////////////////////////

//////////////////////////////Assembeling of matrix from sub-matrices//////////////////////////////
void collect_mat_row(float *A0, float *A1, float *A2, float *A3, float *A, size_t row, size_t col)
{
	size_t vecsize = VECSIZE, i, j;

	if (col / 2 >= vecsize)
	{
		for (i = 0; i < row / 2; i++)
		{
			for (j = 0; j < col / 2 - vecsize; j += vecsize)
			{
				_mm256_storeu_ps(A + i * col + j, _mm256_loadu_ps(A0 + i * (col / 2) + j));
				_mm256_storeu_ps(A + i * col + col / 2 + j, _mm256_loadu_ps(A1 + i * (col / 2) + j));
				_mm256_storeu_ps(A + (i + row / 2) * col + j, _mm256_loadu_ps(A2 + i * (col / 2) + j));
				_mm256_storeu_ps(A + (i + row / 2) * col + col / 2 + j, _mm256_loadu_ps(A3 + i * (col / 2) + j));
			}

			for (; j < col / 2; j++)
			{
				A[i * col + j] = A0[i * (col / 2) + j];
				A[i * col + col / 2 + j] = A1[i * (col / 2) + j];
				A[(i + row / 2) * col + j] = A2[i * (col / 2) + j];
				A[(i + row / 2) * col + col / 2 + j] = A3[i * (col / 2) + j];
			}
		}
	}
	else
	{
		for (i = 0; i < row / 2; i++)
		{
			for (j = 0; j < col / 2; j++)
			{
				A[i * col + j] = A0[i * (col / 2) + j];
				A[i * col + col / 2 + j] = A1[i * (col / 2) + j];
				A[(i + row / 2) * col + j] = A2[i * (col / 2) + j];
				A[(i + row / 2) * col + col / 2 + j] = A3[i * (col / 2) + j];
			}
		}
	}
}
/////////////////////////////End Assembeling of matrix from sub-matrices////////////////////////////

//////////////////////////////////////Printing of matrices//////////////////////////////////////
void print_matrix(float *A, size_t row, size_t col)
{
	for (size_t i = 0; i < row; i++)
	{
		for (size_t j = 0; j < col; j++)
		{
			cout << A[i * col + j] << " ";
		}
		cout << endl;
	}
}
//////////////////////////////////////End Printing of matrices//////////////////////////////////////

/////////////////////////////////////////Mult_AVX Row Major/////////////////////////////////////////
inline float *Mult_AVX(float *A, float *B, size_t rowA, size_t colA_rowB, size_t colB)
{
	//mutex1.lock();
	//cout << "rowA = " << rowA << ", colA_rowB = " << colA_rowB << ", colB = " << colB << endl;
	//mutex1.unlock();

	float *C = new float[rowA * colB];				   // return matrix
	size_t i_block = std::min(I_BLOCKSIZE, (int)rowA); // sizes of blocks.
	size_t j_block = std::min(J_BLOCKSIZE, (int)colB); // Different values are possible.
	size_t k_block = std::min(K_BLOCKSIZE, (int)colA_rowB);

	float *tmp_B = new float[colA_rowB * colB]; //Temporary matrix for storing matrix B.
	size_t index = 0;
	size_t vecsize = VECSIZE;		  // number of floats in __m256
	const size_t unroll_1 = UNROLL_1; // unroll_1 loops parametr
	const size_t unroll_2 = UNROLL_2; // unroll_2 loops parametr

	__m256 sum[unroll_2][unroll_1], vecA[unroll_1], vecB[unroll_2];

	//	mutex1.lock();
		//cout << "before" << endl;
		//print_matrix(B, colA_rowB, colB);

		//Storing matrix B in tmp_B (in other order)
	for (size_t jj = 0; jj < colB; jj += j_block)
	{
		for (size_t kk = 0; kk < colA_rowB; kk += k_block)
		{
			for (size_t j = jj; j < jj + j_block; j += unroll_2 * vecsize)
			{
				for (size_t k = kk; k < kk + k_block; k++)
				{
					for (size_t x = 0; x < unroll_2; x++)
					{
						vecB[x] = _mm256_loadu_ps(B + k * colB + (j + x * vecsize));
						_mm256_storeu_ps(tmp_B + index + x * vecsize, vecB[x]);
					}
					index += unroll_2 * vecsize;
				}
			}
		}
	}
	//cout << "after B" << endl;
	//print_matrix(B, colA_rowB, colB);
	//cout << "after" << endl;
	//print_matrix(tmp_B, colA_rowB, colB);
	//mutex1.unlock();


	//Matrix multiplication
	for (size_t ii = 0; ii < rowA; ii += i_block)
	{
		for (size_t jj = 0; jj < colB; jj += j_block)
		{
			for (size_t kk = 0; kk < colA_rowB; kk += k_block)
			{
				for (size_t i = ii; i < ii + i_block; i += unroll_1)
				{
					for (size_t j = jj; j < jj + j_block; j += unroll_2 * vecsize)
					{
						index = (j - jj) * k_block + kk * j_block + jj * colA_rowB;

						if (kk == 0)
						{
							for (size_t x = 0; x < unroll_1; x++)
							{
								for (size_t y = 0; y < unroll_2; y++)
								{
									sum[y][x] = _mm256_setzero_ps();
								}
							}
						}
						else
						{
							for (size_t x = 0; x < unroll_1; x++)
							{
								for (size_t y = 0; y < unroll_2; y++)
								{
									sum[y][x] = _mm256_loadu_ps(C + (i + x) * colB + (j + y * vecsize));
								}
							}
						}

						for (size_t k = kk; k < kk + k_block; k++)
						{
							for (size_t x = 0; x < unroll_2; x++)
							{
								vecB[x] = _mm256_loadu_ps(tmp_B + index + x * vecsize);
							}

							for (size_t x = 0; x < unroll_1; x++)
							{
								vecA[x] = _mm256_broadcast_ss(A + (i + x) * colA_rowB + k);
								for (size_t y = 0; y < unroll_2; y++)
								{
									sum[y][x] = _mm256_fmadd_ps(vecA[x], vecB[y], sum[y][x]);
								}
							}
							index += unroll_2 * vecsize;
						}

						for (size_t x = 0; x < unroll_1; x++)
						{
							for (size_t y = 0; y < unroll_2; y++)
							{
								_mm256_storeu_ps(C + (i + x) * colB + (j + y * vecsize), sum[y][x]);
							}
						}
					}
				}
			}
		}
	}

	delete[] tmp_B;
	return C;
}

/////////////////////////////Checking that the matrices are equal/////////////////////////////
bool matrix_equal(float *A, float *B, size_t row, size_t col)
{
	float eps = 0.00001;

	for (size_t i = 0; i < row; i++)
	{
		for (size_t j = 0; j < col; j++)
		{
			if (fabs(A[i * col + j] - B[i * col + j]) > eps)
			{
				return false;
			}
		}
	}

	return true;
}
///////////////////////////End Checking that the matrices are equal///////////////////////////

/*
////////////////////////////////Mult_AVX multithreaded_S - Strassen/////////////////////////////////
float *Mult_AVX_multithread_S(float *A, float *B, size_t n, int iter)
{
	static size_t N = n;
	if (n <= N / pow(2, iter))
		return Mult_AVX(A, B, n, n, n);

	thread t[7];

	float *AA[4], *BB[4], *CC[4], *P[7], *C = new float[n * n];

	for (int i = 0; i < 4; i++)
	{
		AA[i] = new float[(n / 2) * (n / 2)];
		BB[i] = new float[(n / 2) * (n / 2)];
	}

	sub_matrix_row(A, n, 0, n / 2 - 1, 0, n / 2 - 1, AA[0]);
	sub_matrix_row(A, n, n / 2, n - 1, 0, n / 2 - 1, AA[1]);
	sub_matrix_row(A, n, 0, n / 2 - 1, n / 2, n - 1, AA[2]);
	sub_matrix_row(A, n, n / 2, n - 1, n / 2, n - 1, AA[3]);

	sub_matrix_row(B, n, 0, n / 2 - 1, 0, n / 2 - 1, BB[0]);
	sub_matrix_row(B, n, n / 2, n - 1, 0, n / 2 - 1, BB[1]);
	sub_matrix_row(B, n, 0, n / 2 - 1, n / 2, n - 1, BB[2]);
	sub_matrix_row(B, n, n / 2, n - 1, n / 2, n - 1, BB[3]);

	t[0] = thread([&P, AA, BB, n, iter]() { P[0] = Mult_AVX_multithread_S(sum(AA[0], AA[3], n / 2), sum(BB[0], BB[3], n / 2), n / 2, iter); });
	t[1] = thread([&P, AA, BB, n, iter]() { P[1] = Mult_AVX_multithread_S(sum(AA[1], AA[3], n / 2), BB[0], n / 2, iter); });
	t[2] = thread([&P, AA, BB, n, iter]() { P[2] = Mult_AVX_multithread_S(AA[0], sub(BB[2], BB[3], n / 2), n / 2, iter); });
	t[3] = thread([&P, AA, BB, n, iter]() { P[3] = Mult_AVX_multithread_S(AA[3], sub(BB[1], BB[0], n / 2), n / 2, iter); });
	t[4] = thread([&P, AA, BB, n, iter]() { P[4] = Mult_AVX_multithread_S(sum(AA[0], AA[2], n / 2), BB[3], n / 2, iter); });
	t[5] = thread([&P, AA, BB, n, iter]() { P[5] = Mult_AVX_multithread_S(sub(AA[1], AA[0], n / 2), sum(BB[0], BB[2], n / 2), n / 2, iter); });
	t[6] = thread([&P, AA, BB, n, iter]() { P[6] = Mult_AVX_multithread_S(sub(AA[2], AA[3], n / 2), sum(BB[1], BB[3], n / 2), n / 2, iter); });

	for (int i = 0; i < 7; i++)
		t[i].join();

	for (int i = 0; i < 4; i++)
	{
		delete[] AA[i];
		delete[] BB[i];
	}

	CC[0] = sum(sub(sum(P[0], P[3], n / 2), P[4], n / 2), P[6], n / 2);
	CC[2] = sum(P[2], P[4], n / 2);
	CC[1] = sum(P[1], P[3], n / 2);
	CC[3] = sum(sub(sum(P[0], P[2], n / 2), P[1], n / 2), P[5], n / 2);

	for (int i = 0; i < 7; i++)
		delete[] P[i];

	collect_mat_row(CC[0], CC[1], CC[2], CC[3], C, n);
	for (int i = 0; i < 4; i++)
		delete[] CC[i];

	return C;
}
////////////////////////////////////End Mult_AVX multithreaded_S////////////////////////////////////
*/

////////////////////////////Mult_AVX multithreaded_SW - Strassen-Winograd///////////////////////////
float *Mult_AVX_multithread_SW(float *A, float *B, size_t rowA, size_t colA_rowB, size_t colB, int iter)
{
	static size_t N = rowA % 2 == 0 ? rowA : rowA + 1;
	static size_t M = colA_rowB % 2 == 0 ? colA_rowB : colA_rowB + 1;
	static size_t K = colB % 2 == 0 ? colB : colB + 1;
	size_t resultRowWithoutPadding = rowA, resultColWithoutPadding = colB;
	size_t blockSize = I_BLOCKSIZE;
	float *paddedA = A, *paddedB = B;
	bool isPadA = false, isPadB = false;

	tuple<float *, size_t, size_t> newA, newB; // padded matrices and their new diemnsions

	if (rowA <= N / pow(2, iter) || colA_rowB <= M / pow(2, iter) || colB <= K / pow(2, iter))
	{
		if ((rowA % blockSize != 0) || (colA_rowB % blockSize != 0))
		{
			newA = paddingBlockSize(A, rowA, colA_rowB);
			paddedA = get<0>(newA);
			rowA = get<1>(newA);
			isPadA = true;
		}

		if ((colA_rowB % blockSize != 0) || (colB % blockSize != 0))
		{
			newB = paddingBlockSize(B, colA_rowB, colB);
			paddedB = get<0>(newB);
			colA_rowB = get<1>(newB);
			colB = get<2>(newB);
			isPadB = true;
		}
		float *resultMatrix = Mult_AVX(paddedA, paddedB, rowA, colA_rowB, colB);

		if (isPadA || isPadB)
		{
			float *resultMatrixWithoutPadding = removePadding(resultMatrix, resultRowWithoutPadding, resultColWithoutPadding, colB);
			delete[] resultMatrix;

			if (isPadA)
			{
				delete[] paddedA;
			}

			if (isPadB)
			{
				delete[] paddedB;
			}

			return resultMatrixWithoutPadding;
		}
		else
		{
			return resultMatrix;
		}
	}

	// pad A matrix into even dimensions
	if (rowA % 2 != 0 || colA_rowB % 2 != 0)
	{
		newA = paddSingleRowAndColumn(A, rowA, colA_rowB);
		paddedA = get<0>(newA);
		rowA = get<1>(newA);
		isPadA = true;
	}

	// pad B matrix into even dimensions
	if (colA_rowB % 2 != 0 || colB % 2 != 0)
	{
		newB = paddSingleRowAndColumn(B, colA_rowB, colB);
		paddedB = get<0>(newB);
		colA_rowB = get<1>(newB);
		colB = get<2>(newB);
		isPadB = true;
	}

	// divide the A and B into sub matrices
	float *AA[4], *BB[4];

	for (size_t i = 0; i < 4; i++)
	{
		AA[i] = new float[(rowA / 2) * (colA_rowB / 2)];
		BB[i] = new float[(colA_rowB / 2) * (colB / 2)];
	}

	// split A matrix
	sub_matrix_row(paddedA, colA_rowB, 0, rowA / 2, 0, colA_rowB / 2, AA[0]);
	sub_matrix_row(paddedA, colA_rowB, 0, rowA / 2, colA_rowB / 2, colA_rowB, AA[1]);
	sub_matrix_row(paddedA, colA_rowB, rowA / 2, rowA, 0, colA_rowB / 2, AA[2]);
	sub_matrix_row(paddedA, colA_rowB, rowA / 2, rowA, colA_rowB / 2, colA_rowB, AA[3]);

	if (isPadA)
	{
		delete[] paddedA;
	}

	// split B matrix
	sub_matrix_row(paddedB, colB, 0, colA_rowB / 2, 0, colB / 2, BB[0]);
	sub_matrix_row(paddedB, colB, 0, colA_rowB / 2, colB / 2, colB, BB[1]);
	sub_matrix_row(paddedB, colB, colA_rowB / 2, colA_rowB, 0, colB / 2, BB[2]);
	sub_matrix_row(paddedB, colB, colA_rowB / 2, colA_rowB, colB / 2, colB, BB[3]);

	if (isPadB)
	{
		delete[] paddedB;
	}

	float *S[8];

	// sums from A sub matrices
	S[0] = sum(AA[2], AA[3], rowA / 2, colA_rowB / 2);
	S[1] = sub(S[0], AA[0], rowA / 2, colA_rowB / 2);
	S[2] = sub(AA[0], AA[2], rowA / 2, colA_rowB / 2);
	S[3] = sub(AA[1], S[1], rowA / 2, colA_rowB / 2);

	// sums from B sub matrices
	S[4] = sub(BB[1], BB[0], colA_rowB / 2, colB / 2);
	S[5] = sub(BB[3], S[4], colA_rowB / 2, colB / 2);
	S[6] = sub(BB[3], BB[1], colA_rowB / 2, colB / 2);
	S[7] = sub(S[5], BB[2], colA_rowB / 2, colB / 2);

	float *P[7];
	thread t[7];
	t[0] = thread([&P, AA, BB, S, rowA, colA_rowB, colB, iter]() { P[0] = Mult_AVX_multithread_SW(S[1], S[5], rowA / 2, colA_rowB / 2, colB / 2, iter); });
	t[1] = thread([&P, AA, BB, S, rowA, colA_rowB, colB, iter]() { P[1] = Mult_AVX_multithread_SW(AA[0], BB[0], rowA / 2, colA_rowB / 2, colB / 2, iter); });
	t[2] = thread([&P, AA, BB, S, rowA, colA_rowB, colB, iter]() { P[2] = Mult_AVX_multithread_SW(AA[1], BB[2], rowA / 2, colA_rowB / 2, colB / 2, iter); });
	t[3] = thread([&P, AA, BB, S, rowA, colA_rowB, colB, iter]() { P[3] = Mult_AVX_multithread_SW(S[2], S[6], rowA / 2, colA_rowB / 2, colB / 2, iter); });
	t[4] = thread([&P, AA, BB, S, rowA, colA_rowB, colB, iter]() { P[4] = Mult_AVX_multithread_SW(S[0], S[4], rowA / 2, colA_rowB / 2, colB / 2, iter); });
	t[5] = thread([&P, AA, BB, S, rowA, colA_rowB, colB, iter]() { P[5] = Mult_AVX_multithread_SW(S[3], BB[3], rowA / 2, colA_rowB / 2, colB / 2, iter); });
	t[6] = thread([&P, AA, BB, S, rowA, colA_rowB, colB, iter]() { P[6] = Mult_AVX_multithread_SW(AA[3], S[7], rowA / 2, colA_rowB / 2, colB / 2, iter); });

	for (size_t i = 0; i < 7; i++)
	{
		t[i].join();
	}

	for (size_t i = 0; i < 8; i++)
	{
		delete[] S[i];
	}

	float *T[2];
	T[0] = sum(P[0], P[1], rowA / 2, colB / 2);
	T[1] = sum(T[0], P[3], rowA / 2, colB / 2);

	for (size_t i = 0; i < 4; i++)
	{
		delete[] AA[i];
		delete[] BB[i];
	}

	float *CC[4];
	CC[0] = sum(P[1], P[2], rowA / 2, colB / 2);
	CC[1] = sum(sum(T[0], P[4], rowA / 2, colB / 2), P[5], rowA / 2, colB / 2);
	CC[2] = sub(T[1], P[6], rowA / 2, colB / 2);
	CC[3] = sum(T[1], P[4], rowA / 2, colB / 2);

	for (size_t i = 0; i < 7; i++)
	{
		delete[] P[i];
	}

	delete[] T[0];
	delete[] T[1];

	// collect into result matrix
	float *C = new float[rowA * colB];
	collect_mat_row(CC[0], CC[1], CC[2], CC[3], C, rowA, colB);
	for (size_t i = 0; i < 4; i++)
	{
		delete[] CC[i];
	}
	//cout << "11" << endl;

	if ((rowA == resultRowWithoutPadding + 1) || (colB == resultColWithoutPadding + 1))
	{
		float *C_withoutPadding = removePadding(C, resultRowWithoutPadding, resultColWithoutPadding, colB);
		delete[] C;
		return C_withoutPadding;
	}
	else
	{
		return C;
	}
}
////////////////////////////////////End Mult_AVX multithreaded_SW///////////////////////////////////

///////////////////Add single row or column of zeros to complete for even matrix dimensions//////////////////////////
inline tuple<float *, size_t, size_t> paddSingleRowAndColumn(float *matrix, size_t row, size_t col)
{
	size_t newRow = (row % 2 == 0) ? row : (row + 1);
	size_t newCol = (col % 2 == 0) ? col : (col + 1);
	size_t i, j, vecsize = VECSIZE, newSize = newRow * newCol;
	float *paddedMatrix = new float[newSize];

	if (col >= vecsize)
	{
		// copy the original matrix
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col - vecsize; j += vecsize)
			{
				_mm256_storeu_ps(paddedMatrix + i * newCol + j, _mm256_loadu_ps(matrix + i * col + j));
			}

			for (; j < col; j++)
			{
				paddedMatrix[i * newCol + j] = matrix[i * col + j];
			}
		}
	}
	else
	{
		// copy the original matrix
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				paddedMatrix[i * newCol + j] = matrix[i * col + j];
			}
		}
	}

	// padd new row
	if (newRow != row)
	{
		if (newSize >= vecsize)
		{
			for (i = row * newCol; i < newSize - vecsize; i += vecsize)
			{
				_mm256_storeu_ps(paddedMatrix + i, _mm256_setzero_ps());
			}

			for (; i < newSize; i++)
			{
				paddedMatrix[i] = 0;
			}
		}
		else
		{
			for (i = row * newCol; i < newSize; i++)
			{
				paddedMatrix[i] = 0;
			}
		}
	}

	// padd new column
	if (newCol != col)
	{
		for (j = col; j < newSize; j += newCol)
		{
			paddedMatrix[j] = 0;
		}
	}

	return { paddedMatrix, newRow, newCol };
}
//////////////////END Add single row or column of zeros to complete for even matrix dimensions/////////////////////////

////////////////////////////Expand the matrix to multiple of block size/////////////////////////////////////
inline tuple<float *, size_t, size_t> paddingBlockSize(float *matrix, size_t row, size_t col)
{
	size_t blockSize = I_BLOCKSIZE, vecSize = VECSIZE, paddedRow = row, paddedCol = col, i, j;

	if ((row % blockSize) != 0)
	{
		paddedRow = ((row / blockSize) + 1) * blockSize;
	}

	if ((col % blockSize) != 0)
	{
		paddedCol = ((col / blockSize) + 1) * blockSize;
	}

	size_t paddedSize = paddedRow * paddedCol;
	float *newA = new float[paddedSize];

	// copy the original matrix
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col - vecSize; j += vecSize)
		{
			_mm256_storeu_ps(newA + i * paddedCol + j, _mm256_loadu_ps(matrix + i * col + j));
		}

		for (; j < col; j++)
		{
			newA[i * paddedCol + j] = matrix[i * col + j];
		}
	}

	// pad the rows with zeros
	if (row != paddedRow)
	{
		for (i = row * paddedCol; i < paddedSize; i += blockSize)
		{
			_mm256_storeu_ps(newA + i, _mm256_setzero_ps());
		}
	}

	// pad the columns with zeros
	if (col != paddedCol)
	{
		for (i = 0; i < paddedRow; i++)
		{
			for (j = col; j < paddedCol; j++)
			{
				newA[i * paddedCol + j] = 0;
			}
		}
	}
	return { newA, paddedRow, paddedCol };
}
/////////////////////////////End Expand the matrix to multiple of block size/////////////////////////////////////

////////////////////////////Remove additional zeros from matrix//////////////////////////////////////////////////
inline float *removePadding(float *matrix, size_t row, size_t col, size_t colAfterPadding)
{
	size_t i, j, vecsize = VECSIZE, sizeWithoutPadding = row * col;
	float *newMatrix = new float[sizeWithoutPadding];

	if (col >= vecsize)
	{
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col - vecsize; j += vecsize)
			{
				_mm256_storeu_ps(newMatrix + i * col + j, _mm256_loadu_ps(matrix + i * colAfterPadding + j));
			}

			for (; j < col; j++)
			{
				newMatrix[i * col + j] = matrix[i * colAfterPadding + j];
			}
		}
	}
	else
	{
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				newMatrix[i * col + j] = matrix[i * colAfterPadding + j];
			}
		}
	}

	return newMatrix;
}
////////////////////////////END Remove additional zeros from matrix//////////////////////////////////////////////////

void test(size_t rowA, size_t colA_rowB, size_t colB)
{
	ofstream result;
	result.open("result.csv", ios_base::app);
	const char comma = ',';
	t_timer tt; // timer
	//cout << "I_BLOCKSIZE = " << I_BLOCKSIZE << comma << " ITER_NUM = " << ITER_NUM << comma << " rowA = " << rowA << comma << " colA_rowB = " << colA_rowB << comma << " colB = " << colB << endl;

	float *A = new float[rowA * colA_rowB];
	float *B = new float[colA_rowB * colB];
	float *C_SW = nullptr, *C_MKL = new float[rowA * colB];
	//int cnt = 0;

	// initialization of matrix A
	for (size_t i = 0; i < rowA; i++)
	{
		for (size_t j = 0; j < colA_rowB; j++)
		{
			//A[i * colA_rowB + j] = 1;
			//A[i * colA_rowB + j] = (i + j) % 5 + 1;
			A[i * colA_rowB + j] = ((float)rand()) / (float)RAND_MAX;
		}
	}

	int cnt = 0;

	// initialization of matrix B
	for (size_t i = 0; i < colA_rowB; i++)
	{
		for (size_t j = 0; j < colB; j++)
		{
			//B[i * colB + j] = 1;
			//B[i * colB + j] = ++cnt;
			B[i * colB + j] = ((float)rand()) / (float)RAND_MAX;
		}
	}

	//Mult_AVX(A, B, rowA, colA_rowB, colB);

	/*************** Check with Dan if we still need this shit *********************/

	//// Mult_AVX multithreaded - Strassen
	//for (int it = 0; it < count_it; it++)
	//{
	//	tt.start();														 //start timer
	//	C1 = Mult_AVX_multithread_S(B2, A2, originalSizes[4], ITER_NUM); //multipy matrix
	//	tt.stop();														 //stop timer
	//	t += tt.get_time();
	//}
	//cout << "Mult_AVX Strassen time = " << t / count_it << endl; //print time
	//double t_avx = t / count_it;

	/*************** Check with Dan if we still need this shit *********************/

	// Mult_AVX multithreaded - Strassen Winograd
	double t = 0;
	size_t count_it = 10;

	for (size_t it = 0; it < count_it; it++)
	{
		tt.start();															   // start timer
		C_SW = Mult_AVX_multithread_SW(A, B, rowA, colA_rowB, colB, ITER_NUM); // multipy matrices
		tt.stop();															   // stop timer
		t += tt.get_time();
	}
	double t_avx_SW = t / count_it;
	cout << "Mult_AVX Strassen-Winograd time = " << t_avx_SW << endl; //print time

	// Intel MKL Realization Row Major
	t = 0;
	for (size_t it = 0; it < count_it; it++)
	{
		tt.start();																																   // start timer
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, colB, colA_rowB, (float)1.0, A, colA_rowB, B, colB, (float)0.0, C_MKL, colB); //https://developer.apple.com/documentation/accelerate/1513264-cblas_sgemm?language=objc
		tt.stop();																																   // stop timer
		t += tt.get_time();
	}
	double t_mkl = t / count_it;
	cout << "MKL time = " << t_mkl << endl; // print time
	/*
	// cout << "Mult_AVX Strassen time / MKL time = " << t_avx / t_mkl << endl;
	cout << "Mult_AVX Strassen-Winograd time / MKL time = " << t_avx_SW / t_mkl << endl;


	//// Checking that the matrices are equal.
	//if (matrix_equal(C_SW, C_MKL, rowA, colB))
	//{
	//	cout << "Matrices are equal." << endl;
	//}
	//else
	//{
	//	cout << "Error of calculations. Matrices are not equal." << endl;
	//}

	result << "I_BLOCKSIZE = " << I_BLOCKSIZE << comma
		<< "J_BLOCKSIZE = " << J_BLOCKSIZE << comma
		<< "K_BLOCKSIZE = " << K_BLOCKSIZE << comma
		<< " ITER_NUM = " << ITER_NUM << comma
		<< " rowA = " << rowA << comma
		<< " colA_rowB = " << colA_rowB << comma
		<< " colB = " << colB << comma
		<< " Mult_AVX Strassen-Winograd time = " << t_avx_SW << comma
		<< " MKL time = " << t_mkl << comma
		<< " Mult_AVX Strassen-Winograd time / MKL time = " << t_avx_SW / t_mkl << endl;
		*/
		//result << "A = " << n << "x" << m << comma << "B = " << k << "x" << l << comma << "Mult_AVX_multithread_SW = " << t_avx_1 << comma << "Intel MKL = " << t_mkl << comma << "Mult_AVX Strassen-Winograd time / MKL time = " << t_avx_1 / t_mkl << endl;
	//size_t rowA, size_t colA_rowB, size_t colB
	result << I_BLOCKSIZE << comma << t_avx_SW << comma << t_mkl << endl;

	//cout << "L1 Instructions Cache Size = " << sysconf(_SC_LEVEL1_ICACHE_SIZE) / pow(2, 10) << "KB" << endl;
	//cout << "L1 Data Cache Size = " << sysconf(_SC_LEVEL1_DCACHE_SIZE) / pow(2, 10) << "KB" << endl;
	//cout << "L2 Cache Size = " << sysconf(_SC_LEVEL2_CACHE_SIZE) / pow(2, 10) << "KB" << endl;
	//cout << "L3 Cache Size = " << sysconf(_SC_LEVEL3_CACHE_SIZE) / pow(2, 10) << "KB" << endl;

// Memory free

	delete[] A;
	delete[] B;
	delete[] C_SW;
	delete[] C_MKL;

	result.close();
}
////////////////////////////////////////////////Main////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	mkl_set_num_threads(1);

	try
	{
		//std::stoi( str )
		size_t rowA = std::stoi(argv[1]);
		size_t colA_rowB = std::stoi(argv[2]);
		size_t colB = std::stoi(argv[3]);
		test(rowA, colA_rowB, colB);

	}
	catch (const std::exception &exc)
	{
		cout << exc.what();
	}

	return 0;
}
//////////////////////////////////////////////End Main///////////////////////////////