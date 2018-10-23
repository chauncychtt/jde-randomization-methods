/*
	jDE algorithm implementation with different randomization methods.
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
	mt19937 generator(seed);

	if(dist == UNIFORM) {
		uniform_real_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == GAUSSIAN) {
		uniform_real_distribution<double> secondDistribution(FUNC_MIN,FUNC_MAX);
		normal_distribution<double> distribution(secondDistribution(generator),1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == CAUCHY) {
		uniform_real_distribution<double> secondDistribution(FUNC_MIN,FUNC_MAX);
		cauchy_distribution<double> distribution(secondDistribution(generator),1.0);
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
	mt19937 generator(seed);

	if(dist == UNIFORM) {
		uniform_real_distribution<double> distribution(0.0,1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == GAUSSIAN) {
		uniform_real_distribution<double> secondDistribution(FUNC_MIN,FUNC_MAX);
		normal_distribution<double> distribution(secondDistribution(generator),1.0);
		rand = distribution(generator);
		rand2 = distribution(generator);
	} else if(dist == CAUCHY) {
		uniform_real_distribution<double> secondDistribution(FUNC_MIN,FUNC_MAX);
		cauchy_distribution<double> distribution(secondDistribution(generator),1.0);
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
	mt19937 generator(seed);
	double rand;
	
	/* Initializing parameters F and CR */
	ind.F = 1.0;
	ind.CR = 0.9;

	/* Randomizing initial population */
	for(i = 0; i < NP; i++) {
		for(j = 0; j < numberDim; j++) {
			if(dist == UNIFORM) {
				uniform_real_distribution<double> distribution(FUNC_MIN,FUNC_MAX);
				x.push_back(distribution(generator));
			} else if(dist == GAUSSIAN) {
				uniform_real_distribution<double> secondDistribution(FUNC_MIN,FUNC_MAX);
				normal_distribution<double> distribution(secondDistribution(generator),1.0);
				rand = distribution(generator);
				while(rand < FUNC_MIN || rand > FUNC_MAX) {
					rand = distribution(generator);
				}
				x.push_back(rand);
			} else if(dist == CAUCHY) {
				uniform_real_distribution<double> secondDistribution(FUNC_MIN,FUNC_MAX);
				cauchy_distribution<double> distribution(secondDistribution(generator),1.0);
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
					rand = logisticMap(rand, 4.0);
				}
				x.push_back(rand);
			} else if(dist == KENT) {
				uniform_real_distribution<double> distribution(0.0,1.0);
				rand = distribution(generator);
				rand = kentMap(rand, 0.7);
				while(rand < FUNC_MIN || rand > FUNC_MAX) {
					rand = distribution(generator);
					rand = kentMap(rand, 0.7);
				}
				x.push_back(rand);
			}
		}
		ind.x = x;
		population.push_back(ind);
		x.clear();
		ind.x.clear();
	}

	/* Main algorithm */
	ofstream outdata;
	outdata.open("Results.ods");
	for(i = 0; i < generations; i++) {
		//cout << "Generation: " << i << endl;
		for(j = 0; j < NP; j++) {
			population[j].F = calculateNewF(population[j].F, dist);
			population[j].CR = calculateNewCR(population[j].CR, dist);
		}

		mutationOperation(dist);
		crossoverOperation(dist);
		selectionOperation();

		if(i == generations-1) {
			for (int i = 0; i < NP; ++i){
				for (int j = 0; j < numberDim; ++j){
					outdata << population[i].x[j] << "\t";
				}
				outdata << endl;			
			}
		}
	}
}

void jDE::mutationOperation(int dist) {

	int i, j, ind1, ind2, ind3;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	uniform_real_distribution<double> secondDistribution(FUNC_MIN,FUNC_MAX);
	uniform_real_distribution<double> uniformDistribution(0.0,(double)NP);
	normal_distribution<double> normalDistribution(secondDistribution(generator),1.0);
	cauchy_distribution<double> cauchyDistribution(secondDistribution(generator),1.0);

	if(dist == UNIFORM) {
		/* Using rand/1 strategy */
		for(i = 0; i < NP; i++) {
			population[i].trialVector.clear();
			for(j = 0; j < population[i].x.size(); j++) {
				/* Vi = Xrand1 + F * (Xrand2 - Xrand3) */
				ind1 = uniformDistribution(generator);
				while(ind1 == j) {
					ind1 = uniformDistribution(generator);
				}
				ind2 = uniformDistribution(generator);
				while(ind2 == j || ind2 == ind1) {
					ind2 = uniformDistribution(generator);
				}
				ind3 = uniformDistribution(generator);
				while(ind3 == j || ind3 == ind1 || ind3 == ind2) {
					ind3 = uniformDistribution(generator);
				}
				double vi = population[i].x[ind1] + population[i].F * (population[i].x[ind2] - population[i].x[ind3]);
				if(vi < FUNC_MIN || vi > FUNC_MAX) {
					vi = secondDistribution(generator);
				}
				population[i].trialVector.push_back(vi);
			}
		}
	} else if(dist == GAUSSIAN) {
		/* Using rand/1 strategy */
		for(i = 0; i < NP; i++) {
			population[i].trialVector.clear();
			for(j = 0; j < population[i].x.size(); j++) {
				/* Vi = Xrand1 + F * (Xrand2 - Xrand3) */
				ind1 = ((int)fmod(fabs(normalDistribution(generator)),NP));
				while(ind1 == j) {
					ind1 = ((int)fmod(fabs(normalDistribution(generator)),NP));
				}
				ind2 = ((int)fmod(fabs(normalDistribution(generator)),NP));
				while(ind2 == j || ind2 == ind1) {
					ind2 = ((int)fmod(fabs(normalDistribution(generator)),NP));
				}
				ind3 = ((int)fmod(fabs(normalDistribution(generator)),NP));
				while(ind3 == j || ind3 == ind1 || ind3 == ind2) {
					ind3 = ((int)fmod(fabs(normalDistribution(generator)),NP));
				}
				double vi = population[i].x[ind1] + population[i].F * (population[i].x[ind2] - population[i].x[ind3]);
				while(vi < FUNC_MIN || vi > FUNC_MAX) {
					vi = normalDistribution(generator);
				}
				//cout << "--------------- " << vi << " -------------" <<endl;
				population[i].trialVector.push_back(vi);
			}
		}
	} else if(dist == CAUCHY) {
		/* Using rand/1 strategy */
		for(i = 0; i < NP; i++) {
			population[i].trialVector.clear();
			for(j = 0; j < population[i].x.size(); j++) {
				/* Vi = Xrand1 + F * (Xrand2 - Xrand3) */
				ind1 = ((int)fmod(fabs(cauchyDistribution(generator)),NP));
				while(ind1 == j) {
					ind1 = ((int)fmod(fabs(cauchyDistribution(generator)),NP));
				}
				ind2 = ((int)fmod(fabs(cauchyDistribution(generator)),NP));
				while(ind2 == j || ind2 == ind1) {
					ind2 = ((int)fmod(fabs(cauchyDistribution(generator)),NP));
				}
				ind3 = ((int)fmod(fabs(cauchyDistribution(generator)),NP));
				while(ind3 == j || ind3 == ind1 || ind3 == ind2) {
					ind3 = ((int)fmod(fabs(cauchyDistribution(generator)),NP));
				}
				double vi = population[i].x[ind1] + population[i].F * (population[i].x[ind2] - population[i].x[ind3]);
				while(vi < FUNC_MIN || vi > FUNC_MAX) {
					vi = cauchyDistribution(generator);
				}
				population[i].trialVector.push_back(vi);
			}
		}
	} else if(dist == LOGISTIC) {
		/* Using rand/1 strategy */
		for(i = 0; i < NP; i++) {
			population[i].trialVector.clear();
			for(j = 0; j < population[i].x.size(); j++) {
				/* Vi = Xrand1 + F * (Xrand2 - Xrand3) */
				ind1 = ((int)fmod(fabs(logisticMap(uniformDistribution(generator), 4.0)),NP));
				cout << ind1 << endl;
				while(ind1 == j) {
					ind1 = ((int)fmod(fabs(logisticMap(uniformDistribution(generator), 4.0)),NP));
				}
				ind2 = ((int)fmod(fabs(logisticMap(uniformDistribution(generator), 4.0)),NP));
				while(ind2 == j || ind2 == ind1) {
					ind2 = ((int)fmod(fabs(logisticMap(uniformDistribution(generator), 4.0)),NP));
				}
				ind3 = ((int)fmod(fabs(logisticMap(uniformDistribution(generator), 4.0)),NP));
				while(ind3 == j || ind3 == ind1 || ind3 == ind2) {
					ind3 = ((int)fmod(fabs(logisticMap(uniformDistribution(generator), 4.0)),NP));
				}
				double vi = population[i].x[ind1] + population[i].F * (population[i].x[ind2] - population[i].x[ind3]);
				while(vi < FUNC_MIN || vi > FUNC_MAX) {
					vi = logisticMap(uniformDistribution(generator), 4.0);
				}
				population[i].trialVector.push_back(vi);
			}
		}
	} else if(dist == KENT) {
		/* Using rand/1 strategy */
		for(i = 0; i < NP; i++) {
			population[i].trialVector.clear();
			for(j = 0; j < population[i].x.size(); j++) {
				/* Vi = Xrand1 + F * (Xrand2 - Xrand3) */
				ind1 = ((int)fmod(fabs(kentMap(uniformDistribution(generator), 0.7)),NP));
				while(ind1 == j) {
					ind1 = ((int)fmod(fabs(kentMap(uniformDistribution(generator), 0.7)),NP));
				}
				ind2 = ((int)fmod(fabs(kentMap(uniformDistribution(generator), 0.7)),NP));
				while(ind2 == j || ind2 == ind1) {
					ind2 = ((int)fmod(fabs(kentMap(uniformDistribution(generator), 0.7)),NP));
				}
				ind3 = ((int)fmod(fabs(kentMap(uniformDistribution(generator), 0.7)),NP));
				while(ind3 == j || ind3 == ind1 || ind3 == ind2) {
					ind3 = ((int)fmod(fabs(kentMap(uniformDistribution(generator), 0.7)),NP));
				}
				double vi = population[i].x[ind1] + population[i].F * (population[i].x[ind2] - population[i].x[ind3]);
				while(vi < FUNC_MIN || vi > FUNC_MAX) {
					vi = kentMap(uniformDistribution(generator), 0.7);
				}
				population[i].trialVector.push_back(vi);
			}
		}
	}

}
	
void jDE::crossoverOperation(int dist) {
	int i, j, jrand;
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	//uniform_real_distribution<double> distribution(0.0,1.0);
	//uniform_int_distribution<int> secondDistribution(0,NP);
	uniform_real_distribution<double> uniformDistribution(0.0,1.0);
	uniform_real_distribution<double> secondDistribution(0.0,(double)NP);
	uniform_real_distribution<double> thirdDistribution(FUNC_MIN,FUNC_MAX);
	normal_distribution<double> normalDistribution(thirdDistribution(generator),1.0);
	cauchy_distribution<double> cauchyDistribution(thirdDistribution(generator),1.0);
	double rand;

	/* Binary crossover */
	for(i = 0; i < NP; i++) {
		if(dist == UNIFORM) {
			jrand = (int)secondDistribution(generator);
		} else if(dist == GAUSSIAN) {
			jrand = (int)normalDistribution(generator);
			jrand = abs(jrand) % NP;
		} else if(dist == CAUCHY) {
			jrand = (int)cauchyDistribution(generator);
			jrand = abs(jrand) % NP;
		} else if(dist == LOGISTIC) {
			jrand = (int)secondDistribution(generator);
			jrand = logisticMap(jrand, 4.0);
			jrand = abs(jrand) % NP;
		} else if(dist == KENT) {
			jrand = (int)secondDistribution(generator);
			jrand = kentMap(jrand, 0.7);
			jrand = abs(jrand) % NP;
		}

		for(j = 0; j < population[i].x.size(); j++) {


			if(dist == UNIFORM) {
				rand = uniformDistribution(generator);
			} else if(dist == GAUSSIAN) {
				rand = normalDistribution(generator);
				rand = (fmod(fabs(rand),RAND_MAX));
			} else if(dist == CAUCHY) {
				rand = cauchyDistribution(generator);
				rand = (fmod(fabs(rand),RAND_MAX));
			} else if(dist == LOGISTIC) {
				rand = uniformDistribution(generator);
				rand = logisticMap(rand, 4.0);
				rand = (fmod(fabs(rand),RAND_MAX));
			} else if(dist == KENT) {
				rand = uniformDistribution(generator);
				rand = kentMap(rand, 0.7);
				rand = (fmod(fabs(rand),RAND_MAX));
			}

			if(rand <= population[i].CR || j == jrand) {
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
		if(rosenbrockFunction(population[i].trialVector) < rosenbrockFunction(population[i].x)) {
			population[i].x.clear();
			population[i].x = population[i].trialVector;
		} else {
			population[i].x = population[i].x;
		}
	}
}

double jDE::sphereFunction(vector<double> x) {

	/* Defining fitness function - Sphere's function
	   Dimensions: d
	   Domain: x ∈ [-5.12, 5.12]
	   Global minimum: (0, ..., 0)
	*/

	int numDim = x.size();
	double sum1 = 0, result;

	for(int i = 0; i < numDim; i++) {
		sum1 += pow(x[i], 2);
	}
	result = sum1;
   	return result;
}

double jDE::ackleyFunction(vector<double> x) {

	/* Defining fitness function - Ackley's function 
       Dimensions: d
       Domain: x ∈ [-32.768, 32.768]
       Global minimum: (0, ..., 0)
	*/

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

double jDE::rastriginFunction(vector<double> x) {

	/* Defining fitness function - Rastrigin's function 
	   Dimensions: d
	   Domain: x ∈ [-5.12, 5.12]
	   Global minimum: (0, ..., 0)
	*/

	int numDim = x.size(), i;
	double result = 0;

	for(i = 0; i < numDim; i++) {
        result += pow(x[i], 2.0) - 10 * cos(2 * M_PI * x[i]) + 10;
    }

    return result;
}

double jDE::rosenbrockFunction(vector<double> x) {
	/* Defining fitness function - Rosenbrock's function 
	   Dimensions: d
	   Domain: x ∈ [-5, 10]
	   Global minimum: (1, ..., 1)
	*/	

	double result = 0.0;
	int numDim = x.size(), i;

	for(i = 0; i < numDim-1; i++) {
        result += 100.0 * pow((x[i+1] - pow(x[i],2.0)),2) + pow((1.0 - x[i]),2);
    }

    return result;
}

double jDE::schafferFunction(vector<double> x) {

	/* Defining fitness function - Schaffer's function
	   Dimensions: d
	   Domain: x ∈ [-100.0, 100.0]
	   Global minimum: (0, ..., 0)
	*/

	int numDim = x.size();
	double sum1 = 0, aux, aux1;

	aux1 = 0;
    for(int j = 0; j < numDim; j++) {
    	aux = aux + (pow(x[j],(double)2));
    }

    aux = pow(aux,(double)0.25);

    for(int j = 0; j < numDim; j++) {
    	aux1 = aux1 + (pow(x[j],(double)2));
    }

    aux1 = pow(aux1,(double)0.1);
    aux1 = pow(sin(50*aux1),(double)2) + 1.0;

    return aux*aux1;
}

double jDE::griewankFunction(vector<double> x) {

	/* Defining fitness function - Griewank's function
	   Dimensions: d
	   Domain: x ∈ [-600.0, 600.0]
	   Global minimum: (0, ..., 0)
	*/

	int numDim = x.size();
	double aux, aux1, aux2 = 0;

    aux = aux1 = 0;
    aux2 = 1;
    for(int j = 0; j < numDim; j++) {
    	aux1 = aux1 + pow((x[j]),(double)2);
    	aux2 = aux2 * cos((((x[j])/sqrt((double)(j+1)))*M_PI)/180);
    }
    aux = (1/(double)4000) * aux1 - aux2 + 1;

    return aux;
}

double jDE::eggholderFunction(vector<double> x) {

	/* Defining fitness function - Eggholder's function 
	   Dimensions: 2
	   Domain: x ∈ [-512, 512]
	   Global minimum: (512, 404.2319)
	*/

	double result = 0.0;
	int numDim = x.size(), i;

	for(i = 0; i < numDim-1; i++) {
		result += -(x[i+1] + 47.0) * sin(sqrt(fabs(x[i+1] + x[i] * 0.5 + 47.0))) + sin(sqrt(fabs(x[i] - (x[i+1] + 47.0)))) * (-x[i]);
	}
	  
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
