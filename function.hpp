#pragma once
#include <iostream>
//#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
//#include <lapack.h>
#include <random>
#include <ctime>
#include <cmath>


using namespace std;
using namespace arma;

bool est_p_voisin(int i, int j, int Dim,const dmat &Position);
double potentiel_chimique_dichotomie(vec valp, int Nb_val, double beta);
double potentiel_chimique_meth (vec valp, int Nb_val, double beta);
double mu_max (double T, int meth);