#include "neuron.h"

// Fast sigmoid function
// f(x) = x / (1 + abs(x))
void Neuron::CalculateOutputs( void )
{
	y = x / (1.0 + std::abs(x)); // Fast sigmoid function, f(x) = x / (1 + abs(x))
	yp = y * (1.0 - y); // Derived fast sigmoid function, f'(x) = f(x) * (1 - f(x))
}

Neuron::Neuron( void )
{
	SetInput(0.0);
}

Neuron::Neuron(double x)
{
	SetInput(x);
}

double Neuron::GetInput( void ) const
{
	return x;
}

void Neuron::SetInput(double x)
{
	this->x = x;
	CalculateOutputs();
}

double Neuron::GetOutput( void ) const
{
	return y;
}

double Neuron::GetDerivedOutput( void ) const
{
	return yp;
}
