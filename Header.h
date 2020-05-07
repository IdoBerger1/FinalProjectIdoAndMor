//
//  Header.h
//  Classes
//
//  Created by Dan Lemberg on 18/05/2019.
//  Copyright © 2019 Dan Lemberg. All rights reserved.
//

#ifndef Header_h
#define Header_h

#include <x86intrin.h>
using namespace std;

#define architecture_256
////////////////////////////////////////////////////////////////////////////////////
#ifdef architecture_128
    #define vectypefloat __m128
#else
#ifdef architecture_256
    #define vectypefloat __m256
#else
#ifdef architecture_512
    #define vectypefloat __m512
#else
    //??????????????????
#endif
#endif
#endif
////////////////////////////////////////////////////////////////////////////////////
class System{
    //...
};
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
class Float{ //Type of scalar field
    float _num;
public:
    class vec;
    
    //constructors
    inline Float();
    inline Float(float& num);
    inline Float(float&& num);
    inline Float(Float& F);
    inline Float(Float&& F);
    
    //assignment
    inline Float& operator = (Float& F);
    inline Float& operator = (Float&& F);
    inline Float& operator = (float& num);
    inline Float& operator = (float&& num);
    inline Float& operator = (vec& V);  //_mm256_storeu_ps
    inline Float& operator = (vec&& V); //_mm256_storeu_ps
    
    //accessors
    inline float data();
    inline float* adress();
    
    //...
    
    class vec{ //Type of AVX vector
        vectypefloat _v;
        
    public:
        
        //constructors
        inline vec();
        inline vec(vectypefloat& v);
        inline vec(vectypefloat&& v);
        inline vec(float *p); // or inline void load(float *p); //_mm256_loadu_ps
        inline vec(Float *p); // or inline void load(Float *p); //_mm256_loadu_ps
        inline vec(vec& V);
        inline vec(vec&& V);
        
        //assignment
        inline vec& operator = (vec& V);
        inline vec& operator = (vec&& V);
        inline vec& operator = (float *p); //_mm256_loadu_ps
        inline vec& operator = (Float *p); //_mm256_loadu_ps
        inline vec& operator = (Float& F); // or inline void set(Float& x); //_mm256_broadcast_ss
        inline vec& operator = (Float&& F);//or inline void set(Float&& x); //_mm256_broadcast_ss
        
        //geters and seters
        inline vectypefloat data();
        inline vectypefloat* adress();
        
        //operators arithmetic
        inline friend vec operator + (vec& A, vec& B); //_mm256_add_ps        //4 times
        inline friend vec operator - (vec& A, vec& B); //_mm256_sub_ps        //4 times
        inline friend vec operator * (vec& A, vec& B); //_mm256_mul_ps        //4 times
        inline friend vec operator / (vec& A, vec& B); //_mm256_div_ps        //4 times
        inline friend vec mul_add(vec& A, vec& B, vec& C); //_mm256_fmadd_ps // return A * B + C
        inline friend vec mul_sub(vec& A, vec& B, vec& C); //_mm256_fmadd_ps // return A * B - C
        
        //...
    };
};

//class Int
//class Long
//template <size_t base> class GF; GF<base>
//class Bool
//class Double
//class Float
//class Complex_Double
//class Complex_Float
////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////
//CRTP (Curiously Recurring Template Pattern)
template <typename scalar, typename T>
class Memory_Block{ //Base type of Matrix
protected:
    scalar* _mat;
    size_t _row;
    size_t _col;
public:
    //constructors
    Memory_Block(); //empty constructor
    Memory_Block( size_t row, size_t col ); //matrix of size (rows_, cols_)
    Memory_Block( size_t row, size_t col, scalar val ); //matrix filled by val
    Memory_Block( size_t row, size_t col, std::string type ); //rand or one matrix ...
    Memory_Block( std::initializer_list<std::initializer_list<scalar>> list_lists );
    Memory_Block( std::vector<std::vector<scalar>>& vec_vecs );
    Memory_Block( std::vector<std::vector<scalar>>&& vec_vecs );
    Memory_Block( const Memory_Block& M ); //lvalue copy constructor
    Memory_Block( Memory_Block& M ); //lvalue copy constructor
    Memory_Block( Memory_Block&& M ); //rvalue copy constructor
    
    //accessors
    inline scalar& operator () ( size_t i, size_t j ); //access to element (i, j)
    inline scalar* operator [] ( size_t i ); //access to row i
    
    //extractors
    T operator () ( size_t row_1, size_t row_2, size_t col_1, size_t col_2 );  //sub-matrix: (row_1 : row_2, col_1 : col_2)
    T sub( std::vector<size_t> row_list,  std::vector<size_t> col_list );
    //sub-matrix: row_list - is a list of row nambers, col_list - is a list of column nambers
    //if (row_list.size() == 0) then - all rows
    //if (col_list.size() == 0) then - all columns
    
    //geters and seters
    inline scalar* data();
    inline size_t rows();
    inline size_t cols();
    
    //assignment
    inline Memory_Block& operator = ( const Memory_Block& M );
    inline Memory_Block& operator = ( Memory_Block& M );
    inline Memory_Block& operator = ( std::vector< std::vector<scalar> >& vec_vecs );
    inline Memory_Block& operator = ( std::vector< std::vector<scalar> >&& vec_vecs );
    inline Memory_Block& operator = ( Memory_Block&& M );
    
    //input, output operators
    friend std::ostream& operator << ( std::ostream& out, Memory_Block& m );
    friend std::ostream& operator << ( std::ostream& out, Memory_Block&& m );
    friend Memory_Block& operator << ( Memory_Block& M, T x);
    friend Memory_Block& operator , ( Memory_Block& M, T x);
    
    //destructor
    ~Memory_Block();
    
    //Slices ??? "Probably class Slice"
    //...
};
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
template <typename scalar>
class Matrix : public Memory_Block<scalar, Matrix<scalar>>{ //Type of Matrix
    //constructors
    using Memory_Block<scalar, Matrix<scalar>>::Memory_Block;
    
    Matrix( const Matrix& M );
    Matrix( Matrix& M );
    Matrix( Matrix&& M );
    
    //operators assignment
    inline Matrix& operator = ( const Matrix& M );
    inline Matrix& operator = ( Matrix& M );
    inline Matrix& operator = ( std::vector<std::vector<scalar> >& vec_vecs );
    inline Matrix& operator = ( std::vector<std::vector<scalar> >&& vec_vecs );
    inline Matrix& operator = ( Matrix&& M );
    
    //operators arithmetic
    friend Matrix& operator += ( Matrix& A, Matrix& B ); //4 times
    friend Matrix& operator -= ( Matrix& A, Matrix& B ); //4 times
    friend Matrix& operator *= ( Matrix& M, scalar c );  //2 times
    friend Matrix& operator += ( Matrix& M, scalar c );  //2 times
    friend Matrix& operator -= ( Matrix& M, scalar c );  //2 times
    friend Matrix& operator /= ( Matrix& M, scalar c );  //2 times
    
    friend Matrix operator + ( Matrix& A, Matrix& B ); //4 times
    friend Matrix operator - ( Matrix& A, Matrix& B ); //4 times
    friend Matrix operator * ( Matrix& A, Matrix& B ); //4 times
    friend Matrix operator + ( Matrix& A, scalar c ); //4 times
    friend Matrix operator - ( Matrix& A, scalar c ); //4 times
    friend Matrix operator * ( Matrix& A, scalar c ); //4 times
    friend Matrix operator / ( Matrix& A, scalar c ); //2 times

    //product with transpose
    friend Matrix product( Matrix&& A, char mode_a, Matrix&& B, char mode_b ); //4 times
    
    //Element-wise product
    friend Matrix dot_product( Matrix&& a, Matrix&& b ); //4 times
    
    Matrix trans(bool inplace);
    Matrix conj(bool inplace);
    Matrix diag();
    //...
};
////////////////////////////////////////////////////////////////////////////////////

template <typename scalar>
class Vector : public Matrix<scalar>{
    //constructors
    Vector(); //empty constructor
    Vector( size_t len ); //Vector of size len
    Vector( size_t len, scalar val ); //Vector filled by val
    Vector( size_t row, size_t col, std::string type ); //rand Vector ...
    Vector( std::initializer_list<scalar> list );
    Vector( std::vector<scalar>& vec );
    Vector( std::vector<scalar>&& vec_vecs );
    Vector( const Vector& V ); //lvalue copy constructor
    Vector( Vector& V ); //lvalue copy constructor
    Vector( Vector&& V ); //rvalue copy constructor
    
    //...
};
#endif /* Header_h */
