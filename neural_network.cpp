#include <iostream>
#include "neural_network.h"

Matrix NeuralNetwork::GetInputMatrix(uint32_t i) const
{
	return m_layers[i].GetInputMatrix();
}

Matrix NeuralNetwork::GetOutputMatrix(uint32_t i) const
{
	return m_layers[i].GetOutputMatrix();
}

Matrix NeuralNetwork::GetDerivedOutputMatrix(uint32_t i) const
{
	return m_layers[i].GetDerivedOutputMatrix();
}

Layer &NeuralNetwork::GetInputLayer( void )
{
	return m_layers[0];
}

const Layer &NeuralNetwork::GetInputLayer( void ) const
{
	return m_layers[0];
}

Layer &NeuralNetwork::GetOutputLayer( void )
{
	return m_layers[m_layers.GetSize() - 1];
}

const Layer &NeuralNetwork::GetOutputLayer( void ) const
{
	return m_layers[m_layers.GetSize() - 1];
}

void NeuralNetwork::UpdateErrors( void )
{
	m_total_error = 0.0;
	for (uint32_t i = 0; i < GetOutputLayerSize(); ++i) {
		double error_delta = GetOutput(i) - m_target_output[i]; // NOTE: This is the COST FUNCTION (you may want to employ more advanced ones later)
		m_errors[i] = error_delta;
		m_total_error += error_delta;
	}
	m_historical_errors.AddLast(m_total_error);
}


void NeuralNetwork::PropagateBackward( void )
{
	mtlArray<Matrix> new_weights;
	new_weights.SetCapacity(m_weights.GetSize());
	Matrix gradient;

	// OUTPUT TO LAST HIDDEN LAYER
	Matrix derived_values_y_to_z = GetOutputLayer().GetDerivedOutputMatrix();
	Matrix gradients_y_to_z = Matrix(derived_values_y_to_z.GetWidth(), derived_values_y_to_z.GetHeight());
	for (uint32_t i = 0; i < uint32_t(m_errors.GetSize()); ++i) {
		gradients_y_to_z(i, 0) = derived_values_y_to_z(i, 0) * m_errors[i];
	}

	int32_t  last_hidden_layer_index     = m_layers.GetSize() - 2;
	Layer   *last_hidden_layer           = &m_layers[last_hidden_layer_index];
	Matrix  *weights_output_to_hidden    = &m_weights[last_hidden_layer_index];
	Matrix   delta_output_to_hidden      = Transpose(Transpose(gradients_y_to_z) * last_hidden_layer->GetOutputMatrix());
	Matrix   new_weights_output_to_hidden = (*weights_output_to_hidden) - delta_output_to_hidden;

	new_weights.Add(new_weights_output_to_hidden);

	gradient = gradients_y_to_z;

	// Moving from last hidden layer down to input layer
	for (uint32_t i = uint32_t(last_hidden_layer_index); i > 0; --i) {
		Layer  *l                 = &m_layers[i];
		Matrix  derived_hidden    = l->GetDerivedOutputMatrix();
		Matrix  derived_gradients = Matrix(l->GetSize(), 1);
		Matrix *weight_matrix     = &m_weights[i];
		Matrix *original_weight   = &m_weights[i - 1];
		for (uint32_t y = 0; y < weight_matrix->GetHeight(); ++y) {
			double sum = 0.0;
			for (uint32_t x = 0; x < weight_matrix->GetWidth(); ++x) {
				sum += gradient(x, 0) * (*weight_matrix)(x, y);
			}
			derived_gradients(y, 0) = sum * l->GetOutput(y);
		}

		Matrix left_neurons = (i - 1 == 0) ? m_layers[i - 1].GetInputMatrix() : m_layers[i - 1].GetOutputMatrix();
		Matrix delta_weights = Transpose(Transpose(derived_gradients) * left_neurons);
		Matrix new_weights_hidden = *original_weight - delta_weights;
		new_weights.Add(new_weights_hidden);
		gradient = derived_gradients;
	}

	assert(new_weights.GetSize() == m_weights.GetSize());

	for (uint32_t i = 0; i < uint32_t(m_weights.GetSize()); ++i) {
		uint32_t j = uint32_t(m_weights.GetSize()) - i - 1;

		assert(m_weights[i].GetWidth() == new_weights[j].GetWidth() && m_weights[i].GetHeight() == new_weights[j].GetHeight());
		m_weights[i] = new_weights[j];
	}
}

NeuralNetwork::NeuralNetwork( void ) : m_topology(), m_layers(), m_weights(), m_bias(0.0), m_total_error(0.0)
{}

NeuralNetwork::NeuralNetwork(const mtlArray<uint32_t> &topology) : NeuralNetwork()
{
	Create(topology);
}

NeuralNetwork::NeuralNetwork(const mtlArray<uint32_t> &topology, double *x) : NeuralNetwork()
{
	Create(topology);
	SetInputs(x);
}

void NeuralNetwork::Create(const mtlArray<uint32_t> &topology)
{
	Destroy();

	assert(topology.GetSize() >= 2);

	m_topology = topology;
	m_weights.Create(topology.GetSize() - 1);
	m_layers.Create(topology.GetSize());
	for (int i = 0; i < m_layers.GetSize(); ++i) {
		m_layers[i].Create(topology[i]);
	}
	for (int i = 0; i < m_weights.GetSize(); ++i) {
		m_weights[i] = Random(m_topology[i+1], m_topology[i]);
	}
	m_target_output.Create(int(GetOutputLayerSize()));
	m_errors.Create(int(GetOutputLayerSize()));
	m_gradients.Create(int(GetOutputLayerSize()));
	for (uint32_t i = 0; i < GetOutputLayerSize(); ++i) {
		m_target_output[i] = 0.0;
		m_errors[i] = 0.0;
	}
	m_bias = 0.0;
}

void NeuralNetwork::Destroy( void )
{
	m_topology.Free();
	m_weights.Free();
	m_layers.Free();
	m_errors.Free();
	m_target_output.Free();
	m_historical_errors.RemoveAll();
	m_bias = 0.0;
	m_total_error = 0.0;
	m_gradients.Free();
}

uint32_t NeuralNetwork::GetInputLayerSize( void ) const
{
	return GetInputLayer().GetSize();
}

uint32_t NeuralNetwork::GetOutputLayerSize( void ) const
{
	return GetOutputLayer().GetSize();
}

void NeuralNetwork::SetInputs(const double *x)
{
	for (uint32_t i = 0; i < GetInputLayerSize(); ++i) {
		GetInputLayer().SetInput(i, x[i]);
	}
}

void NeuralNetwork::SetExpectedOutputs(const double *y)
{
	for (uint32_t i = 0; i < GetOutputLayerSize(); ++i) {
		m_target_output[i] = y[i];
	}
}

double NeuralNetwork::GetInput(uint32_t i) const
{
	return GetInputLayer().GetInput(i);
}

double NeuralNetwork::GetOutput(uint32_t i) const
{
	return GetOutputLayer().GetOutput(i);
}

void NeuralNetwork::FeedForward( void )
{
	for (uint32_t i = 0; i < uint32_t(m_weights.GetSize()); ++i) {

		Matrix neuron_matrix = (i != 0) ? GetOutputMatrix(i) : GetInputMatrix(i);
		Matrix c = neuron_matrix * m_weights[i];

		for (uint32_t x = 0; x < c.GetWidth(); ++x) {
			m_layers[i + 1].SetInput(uint32_t(x), c(x, 0) + m_bias);
		}
	}
}

void NeuralNetwork::Train( void )
{
	FeedForward();
	UpdateErrors();
	PropagateBackward();
}

void NeuralNetwork::PrintToConsole( void )
{
	for (int i = 0; i < m_layers.GetSize(); ++i) {
		Matrix m = i == 0 ? m_layers[i].GetInputMatrix() : m_layers[i].GetOutputMatrix();
		std::cout << "layer " << i << ": {" << std::endl;
		std::cout << "->";
		m.PrintToConsole();
		std::cout << std::endl;
		if (i < m_weights.GetSize()) {
			std::cout << "weight" << std::endl;
			m_weights[i].PrintToConsole();
			std::cout << std::endl;
		}
		std::cout << "}" << std::endl;
	}
	std::cout << "total error: " << m_total_error << std::endl;
}

double NeuralNetwork::GetTotalError( void ) const
{
	return m_total_error;
}

const mtlArray<double> &NeuralNetwork::GetOutputErrors( void ) const
{
	return m_errors;
}
