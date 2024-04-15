#include <iostream>
#include <armadillo>
#include "matplotlibcpp.h"
#include "plot_magnetic.hpp"

namespace plt = matplotlibcpp;
using namespace std;
using namespace arma;

void plot_magnetic_moments(const vec& mx, const vec& my, int Dim) {
    // Prepare data for the quiver plot
    std::vector<double> x, y, u, v;
    for (int i = 0; i < Dim; i++) {
        for (int j = 0; j < Dim; j++) {
            x.push_back(j);
            y.push_back(i);
            u.push_back(mx(i * Dim + j)); // Convert linear index to 2D
            v.push_back(my(i * Dim + j)); // Convert linear index to 2D
        }
    }

    // Set up plot limits and labels
    plt::quiver(x, y, u, v);
    plt::xlim(-1, Dim);
    plt::ylim(-1, Dim);
    plt::title("Magnetic Moment Directions");
    plt::show();
    
}