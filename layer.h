#ifndef LAYER_H
#define LAYER_H

#include <cstdint>
#include "MiniLib/MTL/mtlArray.h"
#include "neuron.h"
#include "matrix.h"

class Layer
{
private:
	mtlArray<Neuron> m_neurons;

public:
	Layer( void );
	explicit Layer(uint32_t size);

	void Create(uint32_t size);
	void Destroy( void );

	uint32_t GetSize( void ) const;

	void   SetInput(uint32_t i, double x); // SetValue
	void   SetInputs(const double *xs);
	void   SetInputs(double x0, ...);
	double GetInput(uint32_t i) const; // GetValue
	double GetOutput(uint32_t i) const; // GetActivatedValue
	double GetDerivedOutput(uint32_t i) const; // GetDerivedValue

	Matrix GetInputMatrix( void ) const; // MatrixifyValues
	Matrix GetOutputMatrix( void ) const; // MatrixifyActivatedValues
	Matrix GetDerivedOutputMatrix( void ) const; // MatrixifyDerivedValues
};

#endif // LAYER_H
