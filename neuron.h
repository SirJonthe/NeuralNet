#ifndef NEURON_H
#define NEURON_H

#include <cmath>

class Neuron
{
private:
	double m_value; // linear x
	double m_activated_value; // f(x) = [0, 1]
	double m_derived_value; // approx deriv of activated value

	double x;  // function input to f(x)        - m_value
	double y;  // functoun output f(x) = [0, 1] - m_activated_value
	double yp; // derivation of y               - m_derived_value

private:
	void CalculateOutputs( void ); // Activate & Derive

public:
	Neuron( void );
	Neuron(double x);

	double GetInput( void ) const; // GetValue
	void   SetInput(double x);     // SetValue

	double GetOutput( void ) const;        // GetActivatedValue
	double GetDerivedOutput( void ) const; // GetDerivedValue

};

#endif // NEURON_H
