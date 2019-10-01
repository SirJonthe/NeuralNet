#include <ctime>
#include <cstdlib>
#include <iostream>
#include "neuron.h"
#include "matrix.h"
#include "layer.h"
#include "neural_network.h"
#include "MiniLib/MTL/mtlArray.h"
#include <SDL/SDL.h>

// https://www.youtube.com/channel/UCg9rw36CJztvSEjmei0hSPQ/videos

int main(int, char**)
{
	mtlArray<uint32_t> topology;
	topology.Create(3);
	topology[0] = 3;
	topology[1] = 2;
	topology[2] = 3;

	double inputs[3] = { 1.0, 0.0, 1.0 };

	NeuralNetwork neural_network(topology);
	neural_network.SetInputs(inputs);
	neural_network.SetExpectedOutputs(inputs); // just to test the learning process (aka. auto-encoder neural network)

	for (uint64_t i = 0; i < 100; ++i) {
		neural_network.Train();
		std::cout << "Total error: " << neural_network.GetTotalError() << std::endl;
	}

	return 0;
}
