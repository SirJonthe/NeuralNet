#include <ctime>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include "MiniLib/MML/mmlRandom.h"
#include "MiniLib/MML/mmlMath.h"
#include "matrix.h"

Matrix::Matrix( void ) : e(), m_width(0), m_height(0)
{}

Matrix::Matrix(uint32_t width, uint32_t height) : e(), m_width(0), m_height(0)
{
	Create(width, height);
}

Matrix::Matrix(uint32_t width, uint32_t height, const double e00, ...) {
	Create(width, height);
	va_list vl;
	va_start(vl, e00);
	(*this)(0,0) = double(e00);
	for (uint32_t x = 1; x < m_width; ++x) { (*this)(x,0) = va_arg(vl, double); }
	for (uint32_t y = 1; y < m_height; ++y) {
		for (uint32_t x = 0; x < m_width; ++x) {
			(*this)(x,y) = va_arg(vl, double);
		}
	}
	va_end(vl);
}

Matrix::Matrix(uint32_t width, uint32_t height, const double * const mat) : e(), m_width(0), m_height(0)
{
	Create(width, height);
	for (uint32_t y = 0; y < m_height; ++y) {
		for (uint32_t x = 0; x < m_width; ++x) {
			(*this)(x, y) = mat[y*m_width + x];
		}
	}
}

Matrix::Matrix(const Matrix &m) : e(m.e), m_width(m.m_width), m_height(m.m_height)
{}

Matrix &Matrix::operator=(const Matrix &m)
{
	if (this != &m) {
		e = m.e;
		m_width = m.m_width;
		m_height = m.m_height;
	}
	return *this;
}

void Matrix::Create(uint32_t width, uint32_t height)
{
	m_width = width;
	m_height = height;
	e.Create(int(m_width * m_height));
	for (uint32_t x = 0; x < m_width; ++x) {
		for (uint32_t y = 0; y < m_height; ++y) {
			(*this)(x, y) = 0.0;
		}
	}
}

void Matrix::Destroy( void )
{
	e.Free();
	m_width = 0;
	m_height = 0;
}

void Matrix::PrintToConsole( void ) const
{
	std::cout << std::setprecision(3) << std::fixed;
	for (uint32_t y = 0; y < m_height; ++y) {
		for (uint32_t x = 0; x < m_width; ++x) {
			std::cout << "  " << (*this)(x, y);
		}
		if (y < m_height - 1) {
			std::cout << std::endl;
		}
	}
}

double &Matrix::operator()(uint32_t x, uint32_t y)
{
	assert(x < m_width && y < m_height);
	return e[y * m_width + x];
}

const double &Matrix::operator()(uint32_t x, uint32_t y) const
{
	assert(x < m_width && y < m_height);
	return e[y * m_width + x];
}

uint32_t Matrix::GetWidth( void ) const
{
	return m_width;
}

uint32_t Matrix::GetHeight( void ) const
{
	return m_height;
}

void Matrix::ToVector(mtlArray<double> &out) const
{
	out.Copy(e);
}

Matrix Transpose(const Matrix &m)
{
	Matrix a(m.GetHeight(), m.GetWidth());
	for (uint32_t y = 0; y < m.GetHeight(); ++y) {
		for (uint32_t x = 0; x < m.GetWidth(); ++x) {
			a(y, x) = m(x, y);
		}
	}
	return a;
}

Matrix Identity(uint32_t width, uint32_t height)
{
	Matrix m(width, height);
	for (uint32_t y = 0; y < m.GetHeight(); ++y) {
		for (uint32_t x = 0; x < m.GetWidth(); ++x) {
			m(x,y) = (x == y) ? 1.0 : 0.0;
		}
	}
	return m;
}

Matrix Random(uint32_t width, uint32_t height)
{
	Matrix m(width, height);
	mmlRandom rand;
	rand.SetSeed(uint64_t(time(nullptr)));
	for (uint32_t y = 0; y < m.GetHeight(); ++y) {
		for (uint32_t x = 0; x < m.GetWidth(); ++x) {
			m(x,y) = double(rand.GetFloat());
		}
	}
	return m;
}

Matrix Collapse(const Matrix &m, uint32_t x, uint32_t y)
{
	Matrix a(m.GetWidth() - 1, m.GetHeight() - 1);
	uint32_t i2 = 0;
	for (; i2 < y; ++i2) {
		uint32_t j2 = 0;
		for (; j2 < x; ++j2)              { a(j2,i2) = m(j2,i2); }
		for (; j2 < m.GetWidth()-1; ++j2) { a(j2,i2) = m(j2+1,i2); }
	}
	for (; i2 < m.GetHeight()-1; ++i2) {
		uint32_t j2 = 0;
		for (; j2 < x; ++j2)              { a(j2,i2) = m(j2,i2+1); }
		for (; j2 < m.GetWidth()-1; ++j2) { a(j2,i2) = m(j2+1,i2+1); }
	}
	return a;
}

double Trace(const Matrix &m)
{
	assert(m.GetWidth() == m.GetHeight());

	double trace = 0.0;
	for (uint32_t i = 0; i < m.GetWidth(); ++i) {
		trace += m(i,i);
	}
	return trace;
}

double Determinant(const Matrix &m)
{
	assert(m.GetWidth() == m.GetHeight());

	// Stopping conditions
	if (m.GetHeight() < 1)  { return 0.0; }
	if (m.GetHeight() == 1) { return m(0,0); }

	// Transform to upper triangular matrix
	Matrix upperTri = m;
	uint32_t counti = 0;
	for(uint32_t y = 0; y < m.GetHeight() - 1; ++y) {
		// Elementary Row Operation I
		if(upperTri(y,y) == 0.0) {
			for(uint32_t k = y; k < m.GetHeight(); ++k) {
				if(upperTri(y,k) != 0.0) {
					for(uint32_t x = 0; x < m.GetWidth(); ++x) {
						mmlSwap(upperTri(x, y), upperTri(x, k));
					}
					k = m.GetHeight();
				}
			}
			++counti;
		}
		// Elementary Row Operation III
		if(upperTri(y,y) != 0.0) {
			for(uint32_t k = y + 1; k < m.GetHeight(); ++k) {
				double factor = -1.0 * upperTri(y,k) /  upperTri(y,y);
				for(uint32_t x = y; x < m.GetWidth(); ++x) {
					upperTri(x,k) = upperTri(x,k) + (factor * upperTri(x, y));
				}
			}
		}
	}

	// Calculate determinant
	double det = 1.0;
	for (uint32_t y = 0; y < m.GetHeight(); ++y) {
		det *= upperTri(y,y);
	}
	if (counti % 2 != 0) {
		det = -det;
	}

	return det;
}

Matrix Invert(const Matrix &m)
{
	assert(m.GetWidth() == m.GetHeight());

	// Hasn't checked the determinant. If det=0 then there is no inverse.
	Matrix inv = m;
	if (inv(0,0) == 0.0) {
		for (uint32_t y = 1; y < m.GetHeight(); ++y) {
			if (inv(0,y) != 0.0) {
				for (uint32_t x = 0; x < m.GetWidth(); ++x) {
					mmlSwap(inv(x,0), inv(x,y));
				}
				break;
			}
		}
	}
	for (uint32_t y = 1; y < m.GetHeight(); ++y) { inv(y,0) /= inv(0,0); } // normalize row 0
	for (uint32_t y = 1; y < m.GetHeight(); ++y) {
		for (uint32_t x = y; x < m.GetWidth(); ++x) { // do a column of L
			double sum = 0.0;
			for (uint32_t k = 0; k < y; ++k) {
				sum += inv(k,x) * inv(y,k);
			}
			inv(y,x) -= sum;
		}
		if (y == m.GetHeight()-1) { continue; }
		for (uint32_t x = y+1; x < m.GetHeight(); ++x) { // do a row of U
			double sum = 0.0;
			for (uint32_t k = 0; k < y; ++k) {
				sum += inv(k,y) * inv(x,k);
			}
			inv(x, y) = (inv(x, y)-sum) / inv(y,y);
		}
	}
	for (uint32_t y = 0; y < m.GetHeight(); ++y) { // invert L
		for (uint32_t x = y; x < m.GetHeight(); ++x) {
			double X = 1.0;
			if (y != x) {
				X = 0.0;
				for (uint32_t k = y; k < x; ++k) {
					X -= inv(k,x) * inv(y,k);
				}
			}
			inv(y,x) = X / inv(x,x);
		}
	}
	for (uint32_t y = 0; y < m.GetHeight(); ++y) { // invert U
		for (uint32_t x = y; x < m.GetHeight(); ++x) {
			if (y == x) { continue; }
			double sum = 0.0;
			for (uint32_t k = y; k < x; ++k) {
				sum += inv(x,k) * ((y==k) ? 1.0 : inv(k,y));
			}
			inv(x, y) = -sum;
		}
	}
	for (uint32_t y = 0; y < m.GetHeight(); ++y) { // final inversion
		for (uint32_t x = 0; x < m.GetHeight(); ++x) {
			double sum = 0.0;
			for (uint32_t k = ((y > x) ? y : x); k < m.GetHeight(); ++k) {
				sum += ((x==k) ? 1.0 : inv(k,x)) * inv(y,k);
			}
			inv(y, x) = sum;
		}
	}
	return inv;
}

double Dot(const double *a, const double *b, uint32_t n)
{
	double d = 0.0;
	for (uint32_t j = 0; j < n; ++j) {
		d += a[j] * b[j];
	}
	return d;
}

Matrix &operator+=(Matrix &l, const Matrix &r)
{
	assert(l.GetWidth() == r.GetWidth() && l.GetHeight() == r.GetHeight());

	for (uint32_t y = 0; y < l.GetHeight(); ++y) {
		for (uint32_t x = 0; x < l.GetWidth(); ++x) {
			l(x,y) += r(x,y);
		}
	}

	return l;
}

Matrix &operator-=(Matrix &l, const Matrix &r)
{
	assert(l.GetWidth() == r.GetWidth() && l.GetHeight() == r.GetHeight());

	for (uint32_t y = 0; y < l.GetHeight(); ++y) {
		for (uint32_t x = 0; x < l.GetWidth(); ++x) {
			l(x,y) -= r(x,y);
		}
	}
	return l;
}

Matrix &operator*=(Matrix &l, Matrix r)
{
	assert(l.GetWidth() == l.GetHeight() && l.GetWidth() == r.GetWidth() && l.GetHeight() == r.GetHeight());

	r = Transpose(r); // transpose right
	for (uint32_t y = 0; y < l.GetHeight(); ++y) {
		const double *a = &l(0, y);
		for (uint32_t x = 0; x < l.GetWidth(); ++x) {
			const double *b = &r(0, x);
			l(x, y) = Dot(a, b, l.GetWidth());
		}
	}
	return l;
}
Matrix &operator*=(Matrix &l, const double r)
{
	for (uint32_t y = 0; y < l.GetHeight(); ++y) {
		for (uint32_t x = 0; x < l.GetWidth(); ++x) {
			l(x,y) *= r;
		}
	}
	return l;
}

Matrix operator+(Matrix l, const Matrix &r) { return (l+=r); }
Matrix operator-(Matrix l, const Matrix &r) { return (l-=r); }
Matrix operator*(const Matrix &l, const Matrix &r)
{
	assert(l.GetWidth() == r.GetHeight());

	const Matrix tr = Transpose(r); // transpose right
	Matrix m(r.GetWidth(), l.GetHeight());
	for (uint32_t y = 0; y < m.GetHeight(); ++y) {
		for (uint32_t x = 0; x < m.GetWidth(); ++x) {
			m(x,y) = Dot(&l(0,y), &tr(0,x), l.GetWidth());
		}
	}
	return m;
}
Matrix operator*(Matrix l, const double &r) { return (l*=r); }
Matrix operator*(const double &l, Matrix r) { return (r*=l); }
