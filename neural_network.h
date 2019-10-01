#ifndef NEURAL_NETWORK_H
#define NEURAL_NETWORK_H

#include "layer.h"
#include "matrix.h"
#include "MiniLib/MTL/mtlList.h"

class NeuralNetwork
{
private:
	mtlArray<uint32_t> m_topology; // number of neurons in each layer (already stored in layers(n).size ???
	mtlArray<Layer>    m_layers;
	mtlArray<Matrix>   m_weights;
	double             m_bias;

	mtlArray<double> m_target_output;
	mtlArray<double> m_errors; // error of each output neuron
	double           m_total_error;
	mtlList<double>  m_historical_errors;

	mtlArray<Matrix> m_gradients;

private:
	Matrix       GetInputMatrix(uint32_t i) const;
	Matrix       GetOutputMatrix(uint32_t i) const;
	Matrix       GetDerivedOutputMatrix(uint32_t i) const;
	Layer       &GetInputLayer( void );
	const Layer &GetInputLayer( void ) const;
	Layer       &GetOutputLayer( void );
	const Layer &GetOutputLayer( void ) const;
	void         FeedForward( void );
	void         UpdateErrors( void );
	void         PropagateBackward( void );

public:
	NeuralNetwork( void );
	explicit NeuralNetwork(const mtlArray<uint32_t> &topology);
	NeuralNetwork(const mtlArray<uint32_t> &topology, double *x);

	void Create(const mtlArray<uint32_t> &topology);
	void Destroy( void );

	uint32_t GetInputLayerSize( void ) const;
	uint32_t GetOutputLayerSize( void ) const;

	void SetInputs(const double *x);
	void SetExpectedOutputs(const double *y);

	double GetInput(uint32_t i) const;
	double GetOutput(uint32_t i) const;

	void Train( void );

	void PrintToConsole( void );

	double                  GetTotalError( void ) const;
	const mtlArray<double> &GetOutputErrors( void ) const;
};

#endif // NEURAL_NETWORK_H
