#include <iostream>
#include <armadillo>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;
using namespace arma;

void plot_magnetic_moments(const vec& mx, const vec& my, int Dim);