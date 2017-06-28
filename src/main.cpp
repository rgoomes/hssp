#include "test.h" // TODO: later remove this include
#include "hssp.h"

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char **argv){
	std::vector<std::string> args(argv+1, argv+argc);

	// TODO: create main funcion in test.cpp and remove this. needs new make instruction
	if(contains(args, std::string("--test")))
		run_validation_tests();
	else {
		long int nodes = 0;
		std::vector<int> solution;
		double volume = hssp(args, solution, nodes);

		if(volume >= 0)
			std::cout << std::scientific << std::setprecision(15) << volume << std::endl;
		for(int p : solution)
			std::cout << p << " ";

		std::cout << std::endl;
	}

	return 0;
}
