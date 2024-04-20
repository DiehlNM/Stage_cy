#include <fstream>
#include <iostream>
#include <armadillo>
#include "save_data.hpp"


using namespace std;
using namespace arma;

void save_magnetic_data(const vec& mx, const vec& my, const vec& mz, const string& filename) {
    ofstream file(filename);
    if (file.is_open()) {
        file << "mx,my,mz\n"; // En-tête du CSV
        for (size_t i = 0; i < mx.size(); ++i) {
            file << mx(i) << "," << my(i) << "," << mz(i) << "\n";
        }
        file.close();
    } else {
        cout << "Unable to open file";
    }
}
