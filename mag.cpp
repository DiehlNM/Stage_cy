#include <iostream>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;

void plot_magnetic_moments(const vector<double>& mx, const vector<double>& my, int Dim) {
    // Prepare data for the quiver plot
    std::vector<double> x, y, u, v;
    for (int i = 0; i < Dim; i++) {
        for (int j = 0; j < Dim; j++) {
            x.push_back(j);
            y.push_back(i);
            u.push_back(mx[i * Dim + j]); // Convert linear index to 2D
            v.push_back(my[i * Dim + j]); // Convert linear index to 2D
        }
    }

    // Set up plot limits and labels
    plt::quiver(x, y, u, v);
    plt::xlim(-1, Dim);
    plt::ylim(-1, Dim);
    plt::title("Magnetic Moment Directions");
    plt::show();
}

int main() {
    int Dim = 4; // Dimension of the lattice
    vector<double> mx(Dim * Dim, 0); // Example vector for x components
    vector<double> my(Dim * Dim, 0); // Example vector for y components

    // Initialize vectors with some pattern
    // Uncomment to simulate a ferromagnetic state
    // fill(mx.begin(), mx.end(), 1);
    // fill(my.begin(), my.end(), 1);

    // Uncomment for antiferromagnetic pattern
    for (int i = 0; i < Dim; ++i) {
        for (int j = 0; j < Dim; ++j) {
            if ((i + j) % 2 == 0) {
                mx[i * Dim + j] = 1;
                my[i * Dim + j] = 1;
            } else {
                mx[i * Dim + j] = -1;
                my[i * Dim + j] = -1;
            }
        }
    }

    plot_magnetic_moments(mx, my, Dim);

    return 0;
}
