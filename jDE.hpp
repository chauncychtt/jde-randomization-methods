/*
	jDE algorithm implementation with different randomization methods.
	Author: Luiza Engler Stadelhofer
	
	To compile: make
	To run: ./test
	
*/

#define Fl 0.1
#define Fu 0.9
#define T1 0.1
#define T2 0.1
#define UNIFORM 1
#define GAUSSIAN 2
#define CAUCHY 3
#define LOGISTIC 4
#define KENT 5
#define RAND_MAX 1.0
#define FUNC_MIN -5.12 // For the sphere function
#define FUNC_MAX 5.12 // For the sphere function

#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <numeric>
#include <cmath>

using namespace std;

struct Individual {
	vector<double> x;
	vector<double> trialVector;
	double F;
	double CR;
};

class jDE {

private:
	int NP;
	vector<Individual> population;
public:
	jDE();
	double calculateNewF(double Fi, int dist);
	double calculateNewCR(double CRi, int dist);
	void initPopulation(int numberDim, int dist, int generations);
	void mutationOperation();
	void crossoverOperation();
	void selectionOperation();
	double sphereFunction(vector<double> x);
	double ackleyFunction(vector<double> x);
	vector<Individual> getPopulation();
	double kentMap(double randomNum, double m);
	double logisticMap(double randomNum, double r);
	double rastriginFunction(vector<double> x);
	double rosenbrockFunction(vector<double> x);
	double eggholderFunction(vector<double> x);
};