#include <iostream>
#include <armadillo>
#include "matplotlibcpp.h"
#include "plot_magnetic.hpp"

using namespace matplotlibcpp;
using namespace std;
using namespace arma;

void plot_magnetic_moments(const vec& mx, const vec& my, const vec& mz, int Dim) {
    // Préparation des données pour les trois plans
    std::vector<double> x, y, z, u, v, w;
    for (int i = 0; i < Dim; i++) {
        for (int j = 0; j < Dim; j++) {
            x.push_back(i);  // Coordonnée x
            y.push_back(j);  // Coordonnée y
            z.push_back(j);  // Coordonnée z utilisée pour le plan yz
            u.push_back(mx(i * Dim + j));  // Composante u du vecteur
            v.push_back(my(i * Dim + j));  // Composante v du vecteur
            w.push_back(mz(i * Dim + j));  // Composante w du vecteur
        }
    }

    // Configuration des subplots
    figure();
    quiver(x, y, u, v);
    xlim(-1, Dim);
    ylim(-1, Dim);
    title("Magnetic Moments in XY Plane");
    //show();

    figure ();
    quiver(x, z, v, w);
    xlim(-1, Dim);
    ylim(-1, Dim);
    title("Magnetic Moments in YZ Plane");
    //show();

    figure ();
    quiver(x, z, u, w);
    xlim(-1, Dim);
    ylim(-1, Dim);
    title("Magnetic Moments in XZ Plane");
    show();

}
