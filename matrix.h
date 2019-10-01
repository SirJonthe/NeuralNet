#ifndef MATRIX_H
#define MATRIX_H

#include <stdarg.h>
#include <assert.h>
#include <cstdint>
#include "MiniLib/MTL/mtlArray.h"

class Matrix
{
private:
	mtlArray<double> e;
	union
	{
		uint32_t m_width;
		uint32_t m_columns;
	};
	union
	{
		uint32_t m_height;
		uint32_t m_rows;
	};

public:
	Matrix( void );
	Matrix(uint32_t width, uint32_t height);
	Matrix(uint32_t width, uint32_t height, const double e00, ...);
	Matrix(uint32_t width, uint32_t height, const double * const mat);
	Matrix(const Matrix &m);
	Matrix &operator=(const Matrix &m);

	void Create(uint32_t width, uint32_t height);
	void Destroy( void );

	void PrintToConsole( void ) const;

	double       &operator()(uint32_t x, uint32_t y);
	const double &operator()(uint32_t x, uint32_t y) const;

	uint32_t GetWidth( void ) const;
	uint32_t GetHeight( void ) const;

	void ToVector(mtlArray<double> &out) const;
};

Matrix Transpose(const Matrix &m);
Matrix Identity(uint32_t width, uint32_t height);
Matrix Random(uint32_t width, uint32_t height);
Matrix Collapse(const Matrix &m);
double Trace(const Matrix &m);
double Determinant(const Matrix &m);
Matrix Invert(const Matrix &m);

Matrix &operator+=(Matrix &l, const Matrix &r);
Matrix &operator-=(Matrix &l, const Matrix &r);
Matrix &operator*=(Matrix &l, Matrix r);
Matrix &operator*=(Matrix &l, const double r);
Matrix  operator+(Matrix l, const Matrix &r);
Matrix  operator-(Matrix l, const Matrix &r);
Matrix  operator*(const Matrix &l, const Matrix &r);
Matrix  operator*(Matrix l, const double &r);
Matrix  operator*(const double &l, Matrix r);

#endif // MATRIX_H
