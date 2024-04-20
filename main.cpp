#include <iostream>
#include <cstdlib>
//#define ARMA_DONT_USE_WRAPPER
//#include "matplotlibcpp.h"
#include <armadillo>
//#include <lapack.h>
#include <random>
#include <fstream>
#include <ctime>
#include <cmath>
#include "function2.hpp"


using namespace std;
using namespace arma;
//using namespace matplotlibcpp;


//bool est_p_voisin(int i, int j, int Dim,const dmat &Position);
//double potentiel_chimique_dichotomie(vec valp, int Nb_val, double beta);
//double potentiel_chimique_meth (vec valp, int Nb_val, double beta);



int main(){

    double mu;
    string fichier("mu_max_T.dat");
    std::ofstream flux(fichier.c_str());

    int N=200;
    int meth;
    vector<double> x(N),y(N);
    cout << "Choisir methode 1 ou 2" << endl;
    cin >> meth ;
    cout << endl;
    
    //for (int i=0; i<N; i++){

        //mu=mu_max(i);
        //flux << i << "\t" << mu << endl;
        //x.at(i)=i;
        //cout << x[i] << "\t";
        
        double a =mu_max(0.000001, meth, 4);

        int result = system("/bin/python3 /media/diehlm/Data/Linux/Stage_CY/plot_moment.py");
        if (result != 0) {
            std::cerr << "Erreur lors de l'exécution du script Python." << std::endl;
            return 1; // Retourne une erreur
        }

    std::cout << "Le script Python a été exécuté avec succès." << std::endl;
        //cout << y[i] << endl;
        //flux << i << "\t" << y[i] << endl;

    //}

        /*cout << " xsize " << x.size() << " ysize  " << y.size() << endl;


        //plot("mu_max_T.dat");
        figure();
        cout << "Test 1 " << endl;
        plot(x,y);
        cout << "Test 2" << endl;
        xlim(0, 199);
        cout << "Test 3" << endl;
        ylim(0, 2);
        cout << "Test 4" << endl;
        show();
        cout << "Test 5" << endl;
        savefig("mu_t.png");
        savefig("mu_T.pdf");*/
    
    /*vector<double> x(200),y(200);

    for(int i=0; i<200; i++){

        x.at(i)=i;
        y.at(i)=i*i;

    }

    figure();
    plot(5);
    cout << " xsize " << x.size() << " ysize  " << y.size() << endl;
    plot(x,y);
    savefig("test.png");
    show();*/

    return 0;

}



/*bool est_p_voisin (int i, int j, int Dim, const dmat &Position){

    double x = Position.at(i,1)-Position.at(j,1); // Différence des abcisses
    double y = Position.at(i,2)-Position.at(j,2); // Différence des ordonnées

    if ( fabs(x)==0 && (fabs(y)==1 || fabs(y)==Dim-1 )) {

        return true; 

    } else if ( fabs(y)==0 && (fabs(x)==1 || fabs(x)==Dim-1) ) {

        return true; 

    } else {

        return false;

    }
    

}


double potentiel_chimique_dichotomie (vec valp, int Nb_val, double beta){

    double mu_1, mu_2, mu_val;
    int const N= valp.size();
    mu_1 = valp(0);
    mu_2 = valp(N-1);
    cout << " mu1 " << mu_1 << " mu2 " << mu_2 << endl;

    double Ne_1=0, Ne_2=0, Ne_val=0;
    int ater=0;

    cout << " abs " << abs(Ne_val-Nb_val) << "\t";

    do {

        Ne_val = Ne_1 = Ne_2 = 0;
        ater ++;
        cout << "ater\t" << ater << endl;
        
        for (int i=0; i<N; i++){

            Ne_1 += 1/(1+exp(beta*(valp(i)-mu_1)));
            Ne_2 += 1/(1+exp(beta*(valp(i)-mu_2)));

        }

        mu_val = mu_1 * (Ne_2-Nb_val)/(Ne_2-Ne_1) + mu_2 * (Nb_val-Ne_1)/(Ne_2-Ne_1);

        for(int i=0; i<N; i++){

            Ne_val += 1/(1+exp(beta*(valp(i)-mu_val)));

        } 

        if (Ne_val<Nb_val){

            mu_1 = mu_val;

        } else {

            mu_2 = mu_val;

        }


    } while (abs(Ne_val-Nb_val) < 10e-10 && ater <500);

    cout << " Nbval " << Nb_val << " Ne1 " << Ne_1 << " Ne2 " << Ne_2 << " Iter " << ater << endl;
    cout << " mu1 " << mu_1 << " mu2 " << mu_2 << endl;

    return mu_val;


}

double potentiel_chimique_meth (vec valp, int Nb_val, double beta){


    double potentiel_chimique = valp(0);
    double step = ((valp(2*Nb_val-1)-valp(0))/Nb_val)/1000;
    double nb_elec =0, w;

    while (nb_elec<Nb_val){

        potentiel_chimique += step;
        nb_elec=0;

        for (int i=0; i<2*Nb_val-1; i++){

            w = exp(beta*(valp(i)-potentiel_chimique));

            nb_elec += 1/(1+w);

        }

    }

    return potentiel_chimique; 


}*/ 