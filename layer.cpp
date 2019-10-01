#include <stdarg.h>
#include "layer.h"

Layer::Layer( void ) : m_neurons()
{}

Layer::Layer(uint32_t size)
{
	Create(size);
}

void Layer::Create(uint32_t size)
{
	m_neurons.Create(int(size));
}

void Layer::Destroy( void )
{
	m_neurons.Free();
}

uint32_t Layer::GetSize( void ) const
{
	return uint32_t(m_neurons.GetSize());
}

void Layer::SetInput(uint32_t i, double x)
{
	m_neurons[i].SetInput(x);
}

void Layer::SetInputs(const double *x)
{
	for (uint32_t i = 0; i < uint32_t(m_neurons.GetSize()); ++i) {
		m_neurons[i].SetInput(x[i]);
	}
}

void Layer::SetInputs(double x0, ...)
{
	va_list vl;
	va_start(vl, x0);
	m_neurons[0] = x0;
	for (int i = 1; i < m_neurons.GetSize(); ++i) {
		m_neurons[i] = va_arg(vl, double);
	}
	va_end(vl);
}

double Layer::GetInput(uint32_t i) const
{
	return m_neurons[i].GetInput();
}

double Layer::GetOutput(uint32_t i) const
{
	return m_neurons[i].GetOutput();
}

double Layer::GetDerivedOutput(uint32_t i) const
{
	return m_neurons[i].GetDerivedOutput();
}

Matrix Layer::GetInputMatrix( void ) const
{
	Matrix m(uint32_t(m_neurons.GetSize()), 1, false);
	for (uint32_t i = 0; i < m.GetWidth(); ++i) {
		m(i, 0) = m_neurons[i].GetInput();
	}
	return m;
}

Matrix Layer::GetOutputMatrix( void ) const
{
	Matrix m(uint32_t(m_neurons.GetSize()), 1, false);
	for (uint32_t i = 0; i < m.GetWidth(); ++i) {
		m(i, 0) = m_neurons[i].GetOutput();
	}
	return m;
}

Matrix Layer::GetDerivedOutputMatrix( void ) const
{
	Matrix m(uint32_t(m_neurons.GetSize()), 1, false);
	for (uint32_t i = 0; i < m.GetWidth(); ++i) {
		m(i, 0) = m_neurons[i].GetDerivedOutput();
	}
	return m;
}
