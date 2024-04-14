#include <iostream>
//#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
//#include <lapack.h>
#include <random>
#include <ctime>
#include <cmath>
#include "function2.hpp"
#include "matplotlibcpp.h"


using namespace std;
using namespace arma;
using namespace matplotlibcpp;

bool est_p_voisin (int i, int j, int Dim, const dmat &Position){

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
    //cout << " mu1 " << mu_1 << " mu2 " << mu_2 << endl;

    double Ne_1=0, Ne_2=0, Ne_val=0;
    int ater=0;

    //cout << " abs " << abs(Ne_val-Nb_val) << "\t";

    do {

        Ne_val = Ne_1 = Ne_2 = 0;
        ater ++;
        //cout << "ater\t" << ater << endl;
        
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

    //cout << " Nbval " << Nb_val << " Ne1 " << Ne_1 << " Ne2 " << Ne_2 << " Iter " << ater << endl;
    //cout << " mu1 " << mu_1 << " mu2 " << mu_2 << endl;

    cout << "\n ATER " << ater << "\n";

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


}

double mu_max (double T, int meth, int Dim){

    std::cout << "Implémentation des variables " << std::endl ;

    // Implémentation des variables  

    //int Dim=4;           // Dimension du réseau 
    int Nb_at = Dim*Dim;    // Nombre d'atomes = Dim*Dim
    double a=1;          // Paramètre de maille 
    double t=1;          // Integrale de transfert
    double E_zero=0;     // Energie ground 

    int k=0;
    
    mat Coord(Dim, Dim); // Tableau avec les coordonnées
    mat Position (Nb_at,3); // Tableau contenant la position de chaque atome numéroté de 0 à Nb_at-1 
                            // et les coordonnées i,j associées, avec l'origine placée en l'atome 0

    // Remplissage du tableau 

    for (int i=0; i<Dim; i++){

        for (int j=0; j<Dim; j++){

            Position(k,0)=k;
            Position(k,1)=i;
            Position(k,2)=j;
            k++;

        }
    }

    //cout << "Test affichage mat " << endl ;
    //cout << Position << endl;

    // Définition du hamiltonien 

    mat H(Nb_at,Nb_at, fill::zeros); //Déclaration d'une matrice nulle

    for (int i=0; i<Nb_at; i++){

        for(int j=0; j<Nb_at; j++){
            
            if (i==j) {

                H(i,j)=E_zero;

            } else if ( est_p_voisin(i,j, Dim, Position) ) {

                H(i,j)=t;

            }

        }
    }

    // Diagonalisation du Hamiltonien 

    double Det_H = det(H);
    
    //cout << "Determinant de H est " << Det_H << endl;

    //cout << endl;
    //cout << "Affichage de l'hamiltonien  \n" << H << endl;

    cx_vec H_eigenvalue;
    cx_mat H_eigenvector;

    eig_gen(H_eigenvalue, H_eigenvector, H);

    //cout << "Valeurs propres " << H_eigenvalue << endl;

    //cout << "Vecteurs propres " << H_eigenvector << endl;


    // Vérifions l'erreur sur le calcul H*Vec_propre=Val_propre*Vec_propre
    
    cx_vec prod1, prod2;

    //for (int i=0; i<Nb_at; i++) {

        //prod1= H*H_eigenvector.col(i);
        //prod2= H_eigenvalue(i)*H_eigenvector.col(i);


        //cout << "Erreur " << i << "= " << prod1 - prod2 << endl;

        
    //}

    
    // Calcul de l'energie moyenne

    cx_mat transposee = H_eigenvector.t();
    cx_double valeur_moyenne , valeur_moyenne2;

    for (int i=0; i<Nb_at; i++){

        for(int j=0; j<Nb_at; j++){

            valeur_moyenne = as_scalar(transposee.row(j)*(H*H_eigenvector.col(i)));

            //cout << "Valeur moyenne i= " << i << " j= " << j << " m= " << valeur_moyenne << endl;
            
            valeur_moyenne2 = as_scalar(transposee.row(j)*(H_eigenvalue(i)*H_eigenvector.col(i)));
            
            //cout << "Valeur moyenne2 i= " << i << " j= " << j << " m= " << valeur_moyenne << endl;



        }
    }


    // Initialisation de valeurs aléatoire pour le valeurs moyennes 
    //de densité de spin de spin
    // L'initialisation est faitede façon à se rapprocher de la solution (antiferro) attendue

    vec val_m_eup(Nb_at, fill::zeros);
    vec val_m_edown(Nb_at, fill::zeros);
    cx_vec si_p (Nb_at, fill::zeros);
    cx_vec si_m ;

    double valeurs;
    std::default_random_engine re(time(0));
    std::uniform_real_distribution<double> distr(0,0.5);
    

    for( int i; i<Nb_at; i++){

        valeurs = distr(re);

        if (i%2==0){

            val_m_eup(i)= 0.5 + valeurs;

        }else {

            val_m_eup(i)= 0.5 - valeurs;

        }

        val_m_edown(i)= 1 - val_m_eup(i);

    }


    vec mx(Nb_at, fill::zeros), my(Nb_at, fill::zeros);

    std::uniform_real_distribution<double> distri(0,1);


    for( int i=0; i<Nb_at; i++){

        valeurs = distri(re);

        mx(i)= valeurs;

    }

    for( int i=0; i<Nb_at; i++){

        valeurs = distri(re);

        my(i)= valeurs;

    }

    si_p.set_real(mx);
    si_p.set_imag(my);

    si_m=conj(si_p);

    //cout << "Si m " << si_m << endl;
    //cout << " sip " << si_p << endl; 
    //cout << "Moyenne electron up " << val_m_eup << endl;
    //cout << "Moyenne des electrons down" << val_m_edown << endl;

    // Calcul auto cohérent des Ni


    double criter_up=1.0 , criter_down=1.0, U=0.1, criter_m=1.0, criter_p=1.0;
    double /*T=500,*/ kB=0.0000861328, eV= 1.9e-19;
    double beta= 1/(kB*T);
    //cout << "beta " << beta << endl;
    int iter=0;
    double mu;

    double energie_totale, energie_fermi, gap;
    mat Val_m_edown, Val_m_eup;
    cx_mat Si_m, Si_p;
    mat H_eup, H_edown;
    cx_vec val_prop_up, val_prop_down, val_prop_tot;
    cx_mat vec_prop_up, vec_propre_down, vec_prop_tot;
    cx_mat H_total(2*Nb_at, 2*Nb_at, fill::zeros);
    vec val_prop;

    const double p=0.3; 
    //double potentiel_chimique;// Potentiel chimique 
    //double step;
    //double nb_elec; // Nombre d'electron 
    //double w_up, w_down; //Variable pour les fonctions de répartition 
    //double fermi_dirac;


    while ( (criter_down > 0.00001) && (criter_up > 0.00001) /*&& (criter_m>0.00001) && (criter_p > 0.00001)*/ ){

        double potentiel_chimique=0;// Potentiel chimique 
        double step=0;
        double nb_elec=0; // Nombre d'electron 
        double w_up=0, w_down=0; //Variable pour les fonctions de répartition 
        double fermi_dirac=0;

        iter +=1;

        if (iter==1000){
            break;
        }

        energie_totale = 0;

        Val_m_eup = diagmat(val_m_edown); // Transformation du vecteur en matrice diagonale
        //cout << "mat up \n" << Val_m_eup << endl; 
        H_eup = H + (U*Val_m_eup);
        //eig_gen( val_prop_up, vec_prop_up, H_eup);



        //for (int i=0; i<Nb_at; i++) {

            //prod1= H_eup*vec_prop_up.col(i);
            //prod2= val_prop_up(i)*vec_prop_up.col(i);


            //cout << "Erreur " << i << "= " << prod1 - prod2 << endl;

        
        //}


        Val_m_edown = diagmat (val_m_eup);

        //cout << " Mat down" << Val_m_edown << endl;
        H_edown = H + (U*Val_m_edown);
        //eig_gen( val_prop_down, vec_propre_down, H_edown); // Verifier les vecteurs propres 

        Si_m = U*diagmat(si_p);
        Si_p = U*diagmat(si_m);

        //Si_m.zeros();
        //Si_p.zeros();

        for (int i=0; i<Nb_at ; i++){

            for (int j=0; j<Nb_at ; j++){

                H_total(i,j) = H_eup(i,j);
                H_total(i+Nb_at, j+Nb_at) = H_edown(i,j);
                H_total(i, j+Nb_at) = (-Si_p(i,j));
                H_total(i+Nb_at, j) = (-Si_m(i,j));

            }

            
        } //cout << " H total " << H_total << endl;



        /*for (int i=0; i<2*Nb_at; i++){

            for (int j=0; j<2*Nb_at; j++){

                if(j==Nb_at){
                    cout << "\t";
                }
                if(H_total(i,j).real()==0){
                    cout << " 0";
                } else {
                    cout << " 1";
                }

            }   if(i==Nb_at-1){
                    cout << endl;
                }
                cout << endl;
        } cout << endl;*/



        //for (int i=0; i<Nb_at; i++) {

            //prod1= H_edown*vec_propre_down.col(i);
            //prod2= val_prop_down(i)*vec_propre_down.col(i);


            //cout << "Erreur " << i << "= " << prod1 - prod2 << endl;

        
        //}

        //cout << "Val propre up " << H_total;
        //
        //cout << "val prpre down" << val_prop_down;

        //val_prop = join_cols(real(val_prop_up),real(val_prop_down));
        //val_prop = sort(val_prop); // Les valeurs propres sont triées par ordre croissant
        //cout << "Val prop triées " << val_prop << endl;

        eig_gen (val_prop_tot, vec_prop_tot, H_total);
        vec val_prop_tot_real = sort(real(val_prop_tot));
        
        //cout << "test1 " << endl;
        //cout << " Vec prop tot " << vec_prop_tot;
        
        /*for (int i=0; i<2*Nb_at; i++){

            for (int j=0; j<2*Nb_at; j++){

                if(j==Nb_at){
                    cout << "\t";
                }
                if(vec_prop_tot(i,j).real()==0){
                    cout << " 0";
                } else {
                    cout << " 1";
                }

            }   if(i==Nb_at-1){
                    cout << endl;
                }
                cout << endl;
        } cout << endl;*/


        // Calcul de l'energie de Fermi 


        energie_fermi = (val_prop_tot_real(Nb_at-1));
        
        gap = real(val_prop_tot_real(Nb_at)-val_prop_tot_real(Nb_at-1));

        //cout << "Energie de fermi " << energie_fermi << "\t" ;
        //cout << "Energie du Gap " << gap << "\t";

        vec nval_moyen_eu(Nb_at, fill::zeros), nval_moyen_ed(Nb_at, fill::zeros);
        cx_vec nsi_m(Nb_at, fill::zeros), nsi_p(Nb_at, fill::zeros);

        //Calcul du potentiel chimique

        if (meth == 1){

            potentiel_chimique = potentiel_chimique_meth(val_prop_tot_real,Nb_at, beta);

        } else if (meth == 2){

            potentiel_chimique = potentiel_chimique_dichotomie (val_prop_tot_real, Nb_at, beta);

        }
        //potentiel_chimique +=- step/2; // ?????????? just -=
        //cout << "Potentiel chimique " << potentiel_chimique << "\t";
        //cout << "Nombre electrons " << nb_elec << "\t";

        //cout << "test2 " << endl;
        /*cout << "nbat " << Nb_at << endl;
        cout << "val " << val_prop_tot.size() << endl;
        cout << "vec " << vec_prop_tot.size() << endl;*/
        

        for (int i=0; i<Nb_at; i++){  //Boucle sur les orbitales 

            for (int j=0; j<2*Nb_at; j++){ // Boucle sur les états propres 


                // Somme sur toutes les énergies infèrieures à Ef 
                // Somme sur toute les énergie * Fermi-Dirac
                
                //cout << "test3 " << endl;
                //cout << "i " << i << " j " << j << endl;

                w_up = exp(beta*(real(val_prop_tot(j))-potentiel_chimique));

                fermi_dirac = 1/(1+w_up);
        //cout << "test31 " << i << " " << j << endl;

                nval_moyen_eu(i) += (real(conj(vec_prop_tot(i,j))*vec_prop_tot(i,j))*fermi_dirac);
        //cout << "test32 " << endl;
                //nsi_p(i) += (conj(vec_prop_tot(j,i))*vec_prop_tot(j,i+Nb_at))*fermi_dirac;


                if(i==0) {

                    energie_totale += real(val_prop_tot(j))*fermi_dirac;

                }



                w_down = exp(beta*(real(val_prop_tot(j))-potentiel_chimique));

        //cout << "val2 " << val_prop_down(j) << " pot " << potentiel_chimique << endl;
                fermi_dirac = 1/(1+w_down);
                //cout << "Fermi2 " << i << " " << j << " " << fermi_dirac << endl;

                nval_moyen_ed(i) += (real(conj(vec_prop_tot(i+Nb_at,j))*vec_prop_tot(i+Nb_at,j))*fermi_dirac);
                
                /*if(i==0) {

                    energie_totale += real(val_prop_tot_real(j))*fermi_dirac;

                }*/
                //cout << "test4 " << endl;

                //nsi_m(i) += (conj(vec_prop_tot(j,i+Nb_at))*vec_prop_tot(j,i))*fermi_dirac;

                //cout << "test5 " << endl;

            } 
                    //cout << " test nval i= " << i << " " << nval_moyen_ed(i) + nval_moyen_eu(i) << "\t";

        } //cout << endl;

            

        //cout << "test " << endl;
        //cout << " i " << i << " j " << j << endl;

        
        for (int i=0; i<Nb_at; i++){

            energie_totale += U*nval_moyen_eu(i)*nval_moyen_ed(i);
        
        }

        //cout << " Ni up " << val_m_eup << endl;
        //cout << " Ni down " << val_m_edown << endl;
        
        criter_up = max (abs(val_m_eup-nval_moyen_eu));

        criter_down = max (abs(val_m_edown-nval_moyen_ed));

        criter_m = max (abs(si_m-nsi_m));

        criter_p = max (abs(si_p-nsi_p));

        //cout << "New Ni " << nval_moyen_ed + nval_moyen_eu << endl;

        val_m_eup = ((1-p)*val_m_eup) + (p*nval_moyen_eu);
        val_m_edown= ((1-p)*val_m_edown) + (p*nval_moyen_ed);
        si_m = ((1-p)*si_m) + (p*nsi_m);
        si_p = ((1-p)*si_p) + (p*nsi_p);

        // Nullifier les Si 
        //si_m.zeros();
        //si_p.zeros();

        //cout << "Iter = " << iter << " criter up " << criter_up << endl; 
        //cout << " criter down" << criter_down << endl;
        //cout << "Ni = " << val_m_eup + val_m_edown << "\t";


        

    }
    
    double mu_maxi;
    mu_maxi = max(val_m_eup-val_m_edown)/2;
    mx = real(si_p);
    my = imag(si_p);
    vector<double> moment(Dim,Dim) ; 
    mat moment_mag (Dim, Dim, arma::fill::zeros);
    

    for (int i=0; i<Dim; i++){

        for( int j=0; j<Dim; j++){

            moment_mag(i,j)=(val_m_eup(i*Dim+j)-val_m_edown(i*Dim+j));

        }
    }

    std::vector<std::vector<double>> data(moment_mag.n_rows, std::vector<double>(moment_mag.n_cols));

    for(size_t i = 0; i < moment_mag.n_rows; ++i) {
        for(size_t j = 0; j < moment_mag.n_cols; ++j) {
            data[i][j] = moment_mag(i, j);
        }
    }    

    imshow(data);
    colorbar();
    show();

    cout << "\n \nFin d'itération \n\n\n" ;
    cout << "Energie de fermi " << energie_fermi << endl;
    cout << "Energie totale " << energie_totale << endl;
    cout << "Valeur du gap " << gap << endl;
    cout << "Nombre ITERATION " << iter << endl;
    cout << "Criter up " << criter_up << endl;
    cout << "Criter down " << criter_down << endl;
    cout << "Criter m " << criter_m << endl;
    cout << "Criter p " << criter_p << endl; 
    cout << "energie totale " << energie_totale << endl;
    cout << "Temperature\t" << T << endl;
    //cout << " Si m\t" << si_m;
    //cout << "Si p\t" << si_p;
    //cout << " Ni up " << val_m_eup << endl;
    //cout << " Ni down " << val_m_edown << endl;
    //cout << " Ni " << val_m_edown + val_m_eup << endl;
    //cout << " si_m " << si_m << endl;
    cout << " Mx " << mx << endl;
    cout << " My " << my << endl; 
    cout << "Moment " << (val_m_eup - val_m_edown)/2 << endl;




    // Display of the magnetic configuration 



 

    //cout << "Moment " << val_m_eup - val_m_edown << endl;
    //cout << " " <<  << endl;
    //cout << " " <<  << endl;

    return abs(mu_maxi);

}
