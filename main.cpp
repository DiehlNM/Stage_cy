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

    double mu;
    string fichier("mu_max_T.dat");
    std::ofstream flux(fichier.c_str());

    int N=200;
    int meth;
    vector<double> x(N),y(N);
    //cout << "Choisir methode 1 ou 2" << endl;
    //cin >> meth ;
    //cout << endl;
    double density=1, a;
    auto start = high_resolution_clock::now();
    auto last_time = start;
    //cout << "test\n" ;

    while (density>=0.5){
        
        a = mu_max(0.000001, 1, 16, density);
        density -= 0.1;
        cout << "Density " << density << endl;

    }

    auto now = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(now - last_time);
    
    auto elapsed = duration_cast<seconds>(now - start);
    auto hour = duration_cast<hours>(elapsed).count();
    auto minute = duration_cast<minutes>(elapsed).count() % 60;
    auto second = elapsed.count() % 60;
    cout << "Temps total d'execution " << hour << "h " << minute << "m " << second << "s" << endl;
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

