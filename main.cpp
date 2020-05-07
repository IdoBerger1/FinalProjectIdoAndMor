#include "HPEML.h"

class test
{
public:
	float num;

	test() : num(0) {}
	test(float _num) : num(_num) {}

	test& operator = (__m256 _register) { _mm256_storeu_ps(&num, _register); return *this; }
	//void method(__m256 _register)
	//{
	//	_mm256_storeu_ps(&num, _register);
	//}

	float* address() { return &num; }
};


std::pair<int, double>  Initialization::times[5];

int main()
{
	try
	{
		size_t row = 2, col = 2;
		//double values[4] = { 3,4,4,9 };
		double values[4] = { 2,-1,1,2 };
		Matrix<Double> mat(row, col, "Customize", values);
		Matrix<Double> calculatedMatrix(row, col);
		Matrix<Double> U(row, col);
		Matrix<Double> sig(row, col);
		Matrix<Double> V(row, col);

		double* eigenValuesJacobi = (double*)malloc(sizeof(double) * calculatedMatrix.rows());
		Int k(10);
		//complex* eigenValues = (complex*)malloc(2*sizeof(complex));

		//eigenValues =  qrAlgorithm(mat, mat.rows(), k);
		
		//JacobiAlgorithm(mat, mat.rows(), calculatedMatrix, eigenValuesJacobi);

		//svdAlgorithm(mat, mat.rows(), U, sig, V);
		size_t row1 = 37, row2 = 39;
		Int F1(8), F2(2);
		//out << "start" << endl;
		Matrix<Int> A(row1, col, F1), B(col, row2, F1);
		cout << A / F2;


		////Matrix<Float> C = A.trans(true);
		//double t = 0;
		//t_timer tt;
		//for (size_t k = 0; k < 10; k++)
		//{
		//	tt.start();		//start timer
		//	A *= B;
		//	tt.stop();	//stop timer
		//	t += tt.get_time();
		//}
		//double time = t / 10;
		//cout << time << endl;
		//cout << A << endl << B;
		//clock_t t = clock();
		//for (size_t i = 0; i < 10; ++i)
		//	A *= B;
		//t = clock() - t;
		//t /= 10;
		//cout << (float)t / CLOCKS_PER_SEC;
		//cout << "Init Start" << endl;
		////Initialization::Init();
		////cout << "Init Done" << endl;
		//cout << "readData Start" << endl;
		//Initialization::readData();
		//cout << "readData Done" << endl;
		//Initialization::setup();
		//initializer_list< initializer_list<int>> l = { {1, 2, 3}, {4, 5, 6} };
//cout << *((l.begin() + 1)->begin() + 1);
		Initialization::getOS();
		//Initialization::readCache();

		return 0;
	}
	catch (const char* caught)
	{
		std::cout << "Got " << caught << std::endl;
	}
}
