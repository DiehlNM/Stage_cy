#include <fstream>
#include <iostream>
#include <filesystem>
#include <armadillo>
#include "save_data.hpp"


using namespace std;
using namespace arma;
//namespace fs = std::filesystem;


void save_magnetic_data(const vec& mx, const vec& my, const vec& mz, int Dim, double density, double T) {
    
    //fs::create_directory("Image");
    //std::filesystem::create_directory("Image");
    string filename = "MF_Data/Data/magnetic_data_dim_" + to_string(Dim) + "_density_" + to_string(density) + "_temperature" + to_string(T) + ".csv" ;
    //string filename2 = "Results"
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
