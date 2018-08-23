/*
	jDE algorithm implementation with different randomisation methods.
	Author: Luiza Engler Stadelhofer

	To compile: make
	To run: ./test
	
*/

#include "jDE.hpp"

jDE::jDE() {
	NP = 50;
}

// Calculating new F parameter
double jDE::calculateNewF(double Fi, int dist) {

	double rand, rand2;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	if(dist == UNIFORM) {
		uniform_real_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == GAUSSIAN) {
		normal_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == CAUCHY) {
		cauchy_distribution<double> distribution(5.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == LOGISTIC) {
		uniform_real_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
		rand = logisticMap(rand, 4.0);
		rand2 = logisticMap(rand2, 4.0);
	} else if(dist == KENT) {
		uniform_real_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
		rand = kentMap(rand, 0.7);
		rand2 = kentMap(rand2, 0.7);
	}

	rand = (fmod(fabs(rand),RAND_MAX));
	rand2 = (fmod(fabs(rand2),RAND_MAX));

	//cout << rand  << " " << rand2 << endl;

	if(rand2 < T1) {
		return Fl + rand * Fu;
	} else {
		return Fi;
	}
}

// Calculating new CR parameter
double jDE::calculateNewCR(double CRi, int dist) {

	double rand, rand2;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);

	if(dist == UNIFORM) {
		uniform_real_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == GAUSSIAN) {
		normal_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == CAUCHY) {
		cauchy_distribution<double> distribution(5.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == LOGISTIC) {
		uniform_real_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
		rand = logisticMap(rand, 4.0);
		rand2 = logisticMap(rand2, 4.0);
	} else if(dist == KENT) {
		uniform_real_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
		rand = kentMap(rand, 0.7);
		rand2 = kentMap(rand2, 0.7);
	}

	rand = (fmod(fabs(rand),RAND_MAX));
	rand2 = (fmod(fabs(rand2),RAND_MAX));

	//cout << rand  << " " << rand2 << endl;

	if(rand2 < T2) {
		return rand;
	} else {
		return CRi;
	}

}

void jDE::initPopulation(int numberDim, int dist, int generations) {

	int i, j;
	vector<double> x;
	Individual ind;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	double rand;
	
	/* Initializing parameters F and CR */
	ind.F = 0.5;
	ind.CR = 0.1;

	/* Randomising initial population */
	for(i = 0; i < NP; i++) {
		for(j = 0; j < numberDim; j++) {
			if(dist == UNIFORM) {
				uniform_real_distribution<double> distribution(0.0,1.0);
				x.push_back(distribution(generator));
			} else if(dist == GAUSSIAN) {
				normal_distribution<double> distribution(0.0,1.0);
				rand = distribution(generator);
				while(rand < FUNC_MIN || rand > FUNC_MAX) {
					rand = distribution(generator);
				}
				x.push_back(rand);
			} else if(dist == CAUCHY) {
				cauchy_distribution<double> distribution(5.0,1.0);
				rand = distribution(generator);
				while(rand < FUNC_MIN || rand > FUNC_MAX) {
					rand = distribution(generator);
				}
				x.push_back(rand);
			} else if(dist == LOGISTIC) {
				uniform_real_distribution<double> distribution(0.0,1.0);
				rand = distribution(generator);
				rand = logisticMap(rand, 4.0);
				while(rand < FUNC_MIN || rand > FUNC_MAX) {
					rand = distribution(generator);
				}
				x.push_back(rand);
			} else if(dist == KENT) {
				uniform_real_distribution<double> distribution(0.0,1.0);
				rand = distribution(generator);
				rand = kentMap(rand, 0.7);
				while(rand < FUNC_MIN || rand > FUNC_MAX) {
					rand = distribution(generator);
				}
				x.push_back(rand);
			}
		}
		ind.x = x;
		population.push_back(ind);
		x.clear();
		ind.x.clear();
	}

	for(i = 0; i < generations; i++) {
		cout << "Generation: " << i << endl;
		for(j = 0; j < NP; j++) {
			population[j].F = calculateNewF(population[j].F, dist);
			population[j].CR = calculateNewCR(population[j].CR, dist);
		}

		mutationOperation();
		crossoverOperation();
		selectionOperation();

		for (int i = 0; i < NP; ++i){
			for (int j = 0; j < numberDim; ++j){
				cout << population[i].x[j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void jDE::mutationOperation() {

	int i, j;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	uniform_real_distribution<double> distribution(0.0,(double)NP);
	uniform_real_distribution<double> secondDistribution(-5.12,5.12);

	/* Using rand/1 strategy */
	for(i = 0; i < NP; i++) {
		population[i].trialVector.clear();
		for(j = 0; j < population[i].x.size(); j++) {
			/* Vi = Xrand1 + F * (Xrand2 - Xrand3) */
			double vi = population[i].x[distribution(generator)] + population[i].F * (population[i].x[distribution(generator)] - population[i].x[distribution(generator)]);
			while(vi < FUNC_MIN || vi > FUNC_MAX) {
				vi = secondDistribution(generator);
			}
			population[i].trialVector.push_back(vi);
		}
	}

}
	
void jDE::crossoverOperation() {
	int i, j;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	uniform_real_distribution<double> distribution(0.0,1.0);
	uniform_real_distribution<double> secondDistribution(0.0,(double)NP);

	/* Binary crossover */
	for(i = 0; i < NP; i++) {
		for(j = 0; j < population[i].x.size(); j++) {
			if(distribution(generator) <= population[i].CR || j == secondDistribution(generator)) {
				population[i].trialVector[j] = population[i].trialVector[j];
			} else {
				population[i].trialVector[j] = population[i].x[j];
			}
		}
	}
}

void jDE::selectionOperation() {

	int i, j;

	for(i = 0; i < NP; i++) {
		if(sphereFunction(population[i].trialVector) < sphereFunction(population[i].x)) {
			population[i].x.clear();
			population[i].x = population[i].trialVector;
		} else {
			population[i].x = population[i].x;
		}
	}
}

double jDE::sphereFunction(vector<double> x) {

	/* Defining fitness function - Sphere's function */
	int numDim = x.size();
	double sum1 = 0, result;

	for(int i = 0; i < numDim; i++) {
		sum1 += pow(x[i], 2);
	}
	result = sum1;
   	return result;
}

double jDE::ackleyFunction(vector<double> x) {

	/* Defining fitness function - Ackley's function */
	int numDim = x.size();
	double b = 0.2, c = 2.0*M_PI, sum1 = 0, sum2 = 0, result = 0, term1 = 0, term2 = 0, term3 = 0, a = 20;

	for(int i = 0; i < numDim; i++) {
		sum1 += pow(x[i], 2);
	}

	for(int i = 0; i < numDim; i++) {
		sum2 += cos(c * x[i]);
	}

	term1 = -1 * a * exp(-1*b*sqrt(1/(double)numDim * sum1));
	term2 = exp(1/(double)numDim * sum2);
	term3 = a + exp(1);
	result = term1 - term2 + term3;

   	return result;

}

vector<Individual> jDE::getPopulation(){
	return population;
}

double jDE::kentMap(double randomNum, double m) {
	if(randomNum <= m) {
		randomNum = randomNum / m;
	} else {
		randomNum = (1 - randomNum) / (1 - m);
	}
	return randomNum;
}

double jDE::logisticMap(double randomNum, double r) {
	randomNum = r * randomNum * (1 - randomNum);
	return randomNum;
}
