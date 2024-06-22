#include <iostream>
#include <cstdlib>
//#define ARMA_DONT_USE_WRAPPER
//#include "matplotlibcpp.h"
#include <armadillo>
//#include <lapack.h>
#include <random>
#include <fstream>
#include <thread>
#include <ctime>
#include <cmath>
#include "function2.hpp"


using namespace std;
using namespace arma;
using namespace std::chrono;

//using namespace matplotlibcpp;



int main(){

    //double mu;
    string fichier("mu_max_T.dat");
    std::ofstream flux(fichier.c_str());

    int N=100;
    int meth;

    //cout << "Choisir methode 1 ou 2" << endl;
    //cin >> meth ;
    //cout << endl;
    double density=1, a;
    auto start = high_resolution_clock::now();
    auto last_time = start;
    //cout << "test\n" ;

    /*while (density>=0.7){
        
        a = mu_max(0.000001, 1, 20, density);
        density -= 0.05;
        cout << "Density " << density << endl;

    }*/

    /*a = mu_max(0.000001, 1, 4, 1);
    a = mu_max(0.000001, 1, 5, 1);
    a = mu_max(0.000001, 1, 6, 1);
    a = mu_max(0.000001, 1, 7, 1);
    a = mu_max(0.000001, 1, 8, 1);
    a = mu_max(0.000001, 1, 9, 1);
    a = mu_max(0.000001, 1, 10, 1)*/;

    string filename = "MF_Data/Plot/moment_max_mf_data_dim_" + to_string(4) + ".csv" ;
    //string filename2 = "Results"
    ofstream file(filename);
    if (file.is_open()) {
        file << "mu_max,T\n"; // En-tête du CSV
    } else {
        cout << "Unable to open file";
    }


    const int num_values = 1000;
    
    // Vecteur pour stocker les valeurs
    std::vector<double> values;
    values.reserve(num_values);
    vector<double> mu(num_values),T(num_values);
    
    // Générer les valeurs équitablement espacées
    for (int i = 0; i < num_values; ++i) {
        values.push_back(static_cast<double>(i) / (num_values - 1));
    }

    /*for (int i=0; i<num_values; i++){
        if (i==0){
            mu[i] = mu_max(0.000001, 1, 4, 1);
            T[i] = pow(10,-i);
            file << mu[i] << "," << T[i] << endl;
        } else {
            mu[i] = mu_max(values[i], 1, 4, 1);
            T[i] = values[i];
            file << mu[i] << "," << T[i] << endl;
        }
        
    }*/ file.close();

    //a = mu_max(0.000001, 1, 21, 1);

    for(int k =0; k<5; k++){
        for (int j=4; j<13; j++ ){
            a = mu_max(0.000001, 1, j, 1);
        }
    }


    auto now = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(now - last_time);
    
    auto elapsed = duration_cast<seconds>(now - start);
    auto hour = duration_cast<hours>(elapsed).count();
    auto minute = duration_cast<minutes>(elapsed).count() % 60;
    auto second = elapsed.count() % 60;
    cout << endl << "Temps total d'execution " << hour << "h " << minute << "m " << second << "s" << endl;
    last_time = now;
    

    cout << "Done" << endl;
        
    //double a =mu_max(0.000001, meth, 4, 1);

    /*int result = system("/bin/python3 /media/diehlm/Data/Linux/Stage_CY/plot_moment.py /media/diehlm/Data/Linux/Stage_CY/criter.py");
    if (result != 0) {
        std::cerr << "Erreur lors de l'exécution du script Python." << std::endl;
        return 1; // Retourne une erreur
    }

    std::cout << "Le script Python a été exécuté avec succès." << std::endl;*/

    return 0;

}

