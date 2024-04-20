#pragma once 
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void save_magnetic_data(const vec& mx, const vec& my, const vec& mz, const string& filename);
