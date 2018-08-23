/*
	jDE algorithm implementation with different randomisation methods.
	Author: Luiza Engler Stadelhofer
	
	To compile: make
	To run: ./test
	
*/

#include "jDE.hpp"

int main() {
	
	jDE ex;
	/* First parameter is number of dimentions, second is type of
	   distribution and third is number of generations */
	ex.initPopulation(100, 1, 100000);

	return 0;
}

