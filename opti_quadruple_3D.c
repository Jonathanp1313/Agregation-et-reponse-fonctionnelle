
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <omp.h>


#ifndef pi
#define pi 3.14159265358979323846
#endif

////// Le code réalise automatiquement les 4 configurations d'agrégation ////////
// Tous solitaires 
// Proies seulement solitaires
// Prédateurs seulement solitaires 
// Proies + Prédateurs grégaires


/////////////////////////////////  Dimension du Tore  ///////////////////////////

double torex = 600.0;
double torey = 600.0;
double torez = 600.0;




////////////////////////////////  Paramètres de base  ///////////////////////////

double delta_t  = 0.1;
int    Nmax     = 3000; // Nombre max de proies
int    P        = 20; // Nombre de prédateur
int    T        = 60; // Temps de simulation
int    Tstart   = 50; // Début de l'interaction proie-prédateur




////////////////////////////////////  Vitesses  /////////////////////////////////

double vinertie = 20.0; // Poids de l'inertie

double vproie = 15.0;
double vpred  = 40.0;




/////////////////////////////// Paramètres Von mises ////////////////////////////

double kappa_proie_ref = 200.0; // Valeur que l'on aurait pour delta_t = 1
double kappa_pred_ref = 200.0;




/////////////////////////// Paramètres de banc des proies ///////////////////////

double Rproie_att = 100.0;   double Cproie_att = 20.0;
double Rproie_al  =  30.0;   double Cproie_al  =  5.0;
double Rproie_rep =  5.0;   double Cproie_rep =  10.0;

#define Rproie_att2 (Rproie_att * Rproie_att)
#define Rproie_al2  (Rproie_al  * Rproie_al)
#define Rproie_rep2 (Rproie_rep * Rproie_rep)





///////////////////////// Paramètres de banc des prédateurs //////////////////////

double Rpred_att = 200.0;    double Cpred_att = 20.0;
double Rpred_al  =  80.0;    double Cpred_al  =  5.0;
double Rpred_rep =  30.0;    double Cpred_rep =  10.0;

#define Rpred_att2 (Rpred_att * Rpred_att)
#define Rpred_al2  (Rpred_al  * Rpred_al)
#define Rpred_rep2 (Rpred_rep * Rpred_rep)




//////////// Décroissance exponentielle de l'attraction pour les bancs ///////////

double alpha_att = 0.05;     // Décroissance exponentielle de l'attraction





////////////////////////// Interactions proies-prédateurs ////////////////////////

double Rfuite  = 200.0;   double Cfuite      = 300.0;   double beta_fuite   = 0.025;
double Rchasse = 200.0;   double Cchasse     = 250.0;   double gamma_chasse = 0.055;
double Rcapt   =   5.0;

#define Rcapt2 (Rcapt * Rcapt)





//////////////////////////// Temps de manip et estomac //////////////////////////

double estomacmax           = 2.0;
double energie_proie        = 1.0;
double decroissance_estomac = 0.95;
int    Tmanip               = 1;





///////////////////////////////////// Divers ////////////////////////////////////

double epsilon   = 0.01;

static volatile unsigned int global_call_counter = 0;


typedef struct { double x, y, z; }              vec3;
typedef struct { double theta, phi; }            dir2;
typedef struct { double x, y, z, theta, phi; }  vec5;

///////////////////////// Construction des loi uniformes /////////////////////////

// Sur [0,1] 

static double runif_c(unsigned int *seed) {
    return (rand_r(seed) + 0.5) / ((double)RAND_MAX + 1.0);
}


// Sur l'ensemble du tore


vec3 runif_ctore(double tx, double ty, double tz, unsigned int *seed) {
    double x = runif_c(seed) * 2.0 * tx - tx;
    double y = runif_c(seed) * 2.0 * ty - ty;
    double z = runif_c(seed) * 2.0 * tz - tz;
    vec3 v = {x, y, z};
    return v;
}

dir2 runif_sphere(unsigned int *seed) {
    double phi   = runif_c(seed) * 2.0 * pi - pi;
    double theta = acos(2.0 * runif_c(seed) - 1.0);
    dir2 d = {theta, phi};
    return d;
}




///////////////////////// Construction loi de Von mises /////////////////////////

double von_mises(double mu, double kappa, unsigned int *seed) {
    double tau   = 1.0 + sqrt(1.0 + 4.0 * kappa * kappa);
    double rho   = (tau - sqrt(2.0 * tau)) / (2.0 * kappa);
    double r     = (1.0 + rho * rho) / (2.0 * rho);
    double z, f, c, theta;
    while (1) {
        double u1 = runif_c(seed);
        double u2 = runif_c(seed);
        double u3 = runif_c(seed);
        z = cos(pi * u1);
        f = (1.0 + r * z) / (r + z);
        c = kappa * (r - f);
        if (c * (2.0 - c) > u2 || log(c / u2) + 1.0 - c >= 0.0) {
            theta = (u3 > 0.5) ? mu + acos(f) : mu - acos(f);
            break;
        }
    }
    return theta;
}




///////////////// Construction de la grille (pour l'optimisation) /////////////////

#define GRID_W 6
#define GRID_H 6
#define GRID_D 6
#define MAX_PER_CELL 200


typedef struct {
    int indices[MAX_PER_CELL];
    int count;
} Cell;

// Conversion coordonnée → cellule (avec wrap toroïdal)
static inline int cell_col(double x) {
    int c = (int)((x + torex) / (2.0 * torex / GRID_W));
    return ((c % GRID_W) + GRID_W) % GRID_W;
}
static inline int cell_row(double y) {
    int r = (int)((y + torey) / (2.0 * torey / GRID_H));
    return ((r % GRID_H) + GRID_H) % GRID_H;
}
static inline int cell_dep(double z) {
    int d = (int)((z + torez) / (2.0 * torez / GRID_D));
    return ((d % GRID_D) + GRID_D) % GRID_D;
}

// Construction de la grille
void build_grid(double *x, double *y, double *z, int N,
                Cell grid[GRID_D][GRID_H][GRID_W]) {

    for (int d = 0; d < GRID_D; d++)
        for (int r = 0; r < GRID_H; r++)
            for (int c = 0; c < GRID_W; c++)
                grid[d][r][c].count = 0;

    for (int i = 0; i < N; i++) {
        if (isnan(x[i])) continue;
        int d = cell_dep(z[i]);
        int r = cell_row(y[i]);
        int c = cell_col(x[i]);
        Cell *cell = &grid[d][r][c];
        if (cell->count < MAX_PER_CELL)
            cell->indices[cell->count++] = i;
    }
}


/////////////////////////// Conditions aux bords du tore //////////////////////////

vec3 veriftore(double xt, double yt, double zt) {
    if (xt >  torex) xt -= 2.0 * torex;
    if (xt < -torex) xt += 2.0 * torex;
    if (yt >  torey) yt -= 2.0 * torey;
    if (yt < -torey) yt += 2.0 * torey;
    if (zt >  torez) zt -= 2.0 * torez;
    if (zt < -torez) zt += 2.0 * torez;
    vec3 v = {xt, yt, zt};
    return v;
}

// Calcule delta_x, delta_y, delta_z en tenant compte du tore




vec3 dxdy_tore(double x1, double y1, double z1,
               double x2, double y2, double z2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    if (dx >  torex) dx -= 2.0 * torex;
    if (dx < -torex) dx += 2.0 * torex;
    if (dy >  torey) dy -= 2.0 * torey;
    if (dy < -torey) dy += 2.0 * torey;
    if (dz >  torez) dz -= 2.0 * torez;
    if (dz < -torez) dz += 2.0 * torez;
    vec3 v = {dx, dy, dz};
    return v;
}



//////////////////////////////// Formation des bancs ///////////////////////////////



// Fonction d'interaction entre congénères lorsqu'il y formation de banc

typedef struct {
    bool   inter;
    double vx_att, vy_att, vz_att;
    double vx_al,  vy_al,  vz_al;
    double vx_rep, vy_rep, vz_rep;
} interaction_banc_resultats;


interaction_banc_resultats interaction_banc_grid(
        double *xmtn, double *ymtn, double *zmtn,
        double *theta_mtn, double *phi_mtn,
        int indiv,
        double Ratt2, double Ral2, double Rrep2, double Catt, double Cal, double Crep,
        Cell grid[GRID_D][GRID_H][GRID_W]) {

    interaction_banc_resultats res = {false,
                                      0.0,0.0,0.0,
                                      0.0,0.0,0.0,
                                      0.0,0.0,0.0};

    int dep0 = cell_dep(zmtn[indiv]);
    int row0 = cell_row(ymtn[indiv]);
    int col0 = cell_col(xmtn[indiv]);

    for (int dd = -1; dd <= 1; dd++) {
        for (int dr = -1; dr <= 1; dr++) {
            for (int dc = -1; dc <= 1; dc++) {

                int dep = ((dep0 + dd) % GRID_D + GRID_D) % GRID_D;
                int r   = ((row0 + dr) % GRID_H + GRID_H) % GRID_H;
                int c   = ((col0 + dc) % GRID_W + GRID_W) % GRID_W;
                Cell *cell = &grid[dep][r][c];

                for (int ci = 0; ci < cell->count; ci++) {
                    int conj = cell->indices[ci];

                    if (conj == indiv)       continue;
                    if (isnan(xmtn[conj]))   continue;

                    vec3 d = dxdy_tore(xmtn[indiv], ymtn[indiv], zmtn[indiv],
                                       xmtn[conj],  ymtn[conj],  zmtn[conj]);

                    double distcongenere2 = d.x*d.x + d.y*d.y + d.z*d.z;

                    double dist = -1.0;


                    ///////////////////// Attraction //////////////////////

                    if (distcongenere2 < Ratt2 && distcongenere2 > Ral2) {
                        if (dist < 0.0) dist = sqrt(distcongenere2);
                        double poids  = Catt * exp(-alpha_att * dist);
                        res.vx_att   += poids * d.x;
                        res.vy_att   += poids * d.y;
                        res.vz_att   += poids * d.z;
                        res.inter     = true;
                    }


                    ///////////////////// Alignement //////////////////////

                    if (distcongenere2 < Ral2
                    &&  distcongenere2 > Rrep2
                    &&  !isnan(theta_mtn[conj])
                    &&  !isnan(phi_mtn[conj])) {

                        double tc = theta_mtn[conj];
                        double pc = phi_mtn[conj];

                        // Vecteur directeur 3D du congénère
                        res.vx_al += Cal * sin(tc) * cos(pc);
                        res.vy_al += Cal * sin(tc) * sin(pc);
                        res.vz_al += Cal * cos(tc);
                        res.inter  = true;
                    }


                    ///////////////////// Répulsion //////////////////////

                    if (distcongenere2 < Rrep2) {
                        if (dist < 0.0) dist = sqrt(distcongenere2);
                        double poids  = Crep / (dist + epsilon);
                        res.vx_rep   -= poids * d.x;
                        res.vy_rep   -= poids * d.y;
                        res.vz_rep   -= poids * d.z;
                        res.inter     = true;
                    }
                }
            }
        }
    }

    return res;
}

//////////////////////////////// Déplacement des proies ///////////////////////////////


vec5 deplacement_proie(int i, bool predateurs_actifs,
                       double *xproie_mtn, double *yproie_mtn, double *zproie_mtn,
                       double *xpred_mtn,  double *ypred_mtn,  double *zpred_mtn,
                       double *thetaproie_mtn, double *phiproie_mtn,
                       Cell grid_proie[GRID_D][GRID_H][GRID_W],
                       Cell grid_pred [GRID_D][GRID_H][GRID_W],
                       int N, int t, bool banc_proie, double kappa_proie, unsigned int *seed) {


    //////////////////////// 1) Inertie ////////////////////////

    double theta_prev = thetaproie_mtn[i];
    double phi_prev   = phiproie_mtn[i];

    double vx_inertie = vinertie * sin(theta_prev) * cos(phi_prev);
    double vy_inertie = vinertie * sin(theta_prev) * sin(phi_prev);
    double vz_inertie = vinertie * cos(theta_prev);




    //////////////////////// 2) Fuite ///////////////////////////

    double vx_fuite = 0.0, vy_fuite = 0.0, vz_fuite = 0.0;

    if (predateurs_actifs) {
        int dep0 = cell_dep(zproie_mtn[i]);
        int row0 = cell_row(yproie_mtn[i]);
        int col0 = cell_col(xproie_mtn[i]);

        for (int dd = -1; dd <= 1; dd++) {
            for (int dr = -1; dr <= 1; dr++) {
                for (int dc = -1; dc <= 1; dc++) {

                    int dep = ((dep0 + dd) % GRID_D + GRID_D) % GRID_D;
                    int r   = ((row0 + dr) % GRID_H + GRID_H) % GRID_H;
                    int c   = ((col0 + dc) % GRID_W + GRID_W) % GRID_W;
                    Cell *cell = &grid_pred[dep][r][c];

                    for (int ci = 0; ci < cell->count; ci++) {
                        int k = cell->indices[ci];
                        vec3 d = dxdy_tore(xproie_mtn[i], yproie_mtn[i], zproie_mtn[i],
                                           xpred_mtn[k],  ypred_mtn[k],  zpred_mtn[k]);
                        double dist2 = d.x*d.x + d.y*d.y + d.z*d.z;
                        if (dist2 < Rfuite * Rfuite) {
                            double poids = Cfuite * exp(-beta_fuite * sqrt(dist2));
                            vx_fuite -= poids * d.x;
                            vy_fuite -= poids * d.y;
                            vz_fuite -= poids * d.z;
                        }
                    }
                }
            }
        }
    }




    /////////////////// 3) Interaction de banc ////////////////////

    double vx_banc = 0.0, vy_banc = 0.0, vz_banc = 0.0;

    if (banc_proie) {
        interaction_banc_resultats res = interaction_banc_grid(
            xproie_mtn, yproie_mtn, zproie_mtn,
            thetaproie_mtn, phiproie_mtn,
            i, Rproie_att2, Rproie_al2, Rproie_rep2,
            Cproie_att, Cproie_al, Cproie_rep,
            grid_proie);

        ///////////////// Attraction /////////////////

        vx_banc += res.vx_att;
        vy_banc += res.vy_att;
        vz_banc += res.vz_att;


        ///////////////// Alignement /////////////////

        if (res.vx_al != 0.0 || res.vy_al != 0.0 || res.vz_al != 0.0) {
            vx_banc += res.vx_al;
            vy_banc += res.vy_al;
            vz_banc += res.vz_al;
        }


        ///////////////// Répulsion //////////////////

        vx_banc += res.vx_rep;
        vy_banc += res.vy_rep;
        vz_banc += res.vz_rep;
    }




    //////////////////// 4) Somme vectorielle ////////////////////
    ////////////////////// Extraction angles /////////////////////
    //////////// Application d'une vitesse constante /////////////

    double vx_tot = vx_inertie + vx_banc + vx_fuite;
    double vy_tot = vy_inertie + vy_banc + vy_fuite;
    double vz_tot = vz_inertie + vz_banc + vz_fuite;

    // Extraction des angles
    double norme = sqrt(vx_tot*vx_tot + vy_tot*vy_tot + vz_tot*vz_tot);
    double theta_tot, phi_tot;

    if (norme > epsilon) {
        theta_tot = acos(fmax(-1.0, fmin(1.0, vz_tot / norme)));
        double rho_xy = sqrt(vx_tot*vx_tot + vy_tot*vy_tot);
        phi_tot = (rho_xy > epsilon) ? atan2(vy_tot, vx_tot) : phi_prev;
    } else {
        // Vecteur nul : on conserve la direction précédente
        theta_tot = theta_prev;
        phi_tot   = phi_prev;
    }

    // Bruit directionnel Von Mises
    theta_tot = von_mises(theta_tot, kappa_proie, seed);
    phi_tot   = von_mises(phi_tot, kappa_proie, seed);

    // Mise à jour de la position
    double xn = xproie_mtn[i] + vproie * sin(theta_tot) * cos(phi_tot) * delta_t;
    double yn = yproie_mtn[i] + vproie * sin(theta_tot) * sin(phi_tot) * delta_t;
    double zn = zproie_mtn[i] + vproie * cos(theta_tot)                * delta_t;

    vec3 pos = veriftore(xn, yn, zn);

    vec5 result = {pos.x, pos.y, pos.z, theta_tot, phi_tot};
    return result;
}

/////////////////////////////// Déplacement des prédateurs //////////////////////////////


vec5 deplacement_pred(int k, bool predateurs_actifs,
                      double *xproie_mtn, double *yproie_mtn, double *zproie_mtn,
                      double *xpred_mtn,  double *ypred_mtn,  double *zpred_mtn,
                      double *thetapred_mtn, double *phipred_mtn,
                      double *estomac, double *tempsmange,
                      Cell grid_proie[GRID_D][GRID_H][GRID_W],
                      Cell grid_pred [GRID_D][GRID_H][GRID_W],
                      int N, int t, bool banc_pred, double kappa_pred, unsigned int *seed) {


    //////////////////////// 1) Inertie ////////////////////////

    double theta_prev = thetapred_mtn[k];
    double phi_prev   = phipred_mtn[k];

    double vx_inertie = vinertie * sin(theta_prev) * cos(phi_prev);
    double vy_inertie = vinertie * sin(theta_prev) * sin(phi_prev);
    double vz_inertie = vinertie * cos(theta_prev);




    //////////////////////// 2) Chasse ///////////////////////////

    double vx_chasse = 0.0, vy_chasse = 0.0, vz_chasse = 0.0;

    if (predateurs_actifs && estomac[k] < estomacmax && tempsmange[k] >= Tmanip) {
        int dep0 = cell_dep(zpred_mtn[k]);
        int row0 = cell_row(ypred_mtn[k]);
        int col0 = cell_col(xpred_mtn[k]);

        for (int dd = -1; dd <= 1; dd++) {
            for (int dr = -1; dr <= 1; dr++) {
                for (int dc = -1; dc <= 1; dc++) {

                    int dep = ((dep0 + dd) % GRID_D + GRID_D) % GRID_D;
                    int r   = ((row0 + dr) % GRID_H + GRID_H) % GRID_H;
                    int c   = ((col0 + dc) % GRID_W + GRID_W) % GRID_W;
                    Cell *cell = &grid_proie[dep][r][c];

                    for (int ci = 0; ci < cell->count; ci++) {
                        int i = cell->indices[ci];
                        if (isnan(xproie_mtn[i])) continue;

                        vec3 d = dxdy_tore(xpred_mtn[k], ypred_mtn[k], zpred_mtn[k],
                                           xproie_mtn[i], yproie_mtn[i], zproie_mtn[i]);
                        double dist2 = d.x*d.x + d.y*d.y + d.z*d.z;
                        if (dist2 < Rchasse * Rchasse) {
                            double poids = Cchasse * exp(-gamma_chasse * sqrt(dist2));
                            vx_chasse += poids * d.x;
                            vy_chasse += poids * d.y;
                            vz_chasse += poids * d.z;
                        }
                    }
                }
            }
        }
    }




    /////////////////// 3) Interaction de banc ////////////////////

    double vx_banc = 0.0, vy_banc = 0.0, vz_banc = 0.0;

    if (banc_pred) {
        interaction_banc_resultats res = interaction_banc_grid(
            xpred_mtn, ypred_mtn, zpred_mtn,
            thetapred_mtn, phipred_mtn,
            k, Rpred_att2, Rpred_al2, Rpred_rep2,
            Cpred_att, Cpred_al, Cpred_rep,
            grid_pred);

        vx_banc = res.vx_att + res.vx_al + res.vx_rep;
        vy_banc = res.vy_att + res.vy_al + res.vy_rep;
        vz_banc = res.vz_att + res.vz_al + res.vz_rep;
    }




    //////////////////// 4) Somme vectorielle ////////////////////
    ////////////////////// Extraction angles /////////////////////
    //////////// Application d'une vitesse constante /////////////

    double vx_tot = vx_inertie + vx_banc + vx_chasse;
    double vy_tot = vy_inertie + vy_banc + vy_chasse;
    double vz_tot = vz_inertie + vz_banc + vz_chasse;

    // Extraction des angles
    double norme = sqrt(vx_tot*vx_tot + vy_tot*vy_tot + vz_tot*vz_tot);
    double theta_tot, phi_tot;

    if (norme > epsilon) {
        theta_tot = acos(fmax(-1.0, fmin(1.0, vz_tot / norme)));
        double rho_xy = sqrt(vx_tot*vx_tot + vy_tot*vy_tot);
        phi_tot = (rho_xy > epsilon) ? atan2(vy_tot, vx_tot) : phi_prev;
    } else {
        theta_tot = theta_prev;
        phi_tot   = phi_prev;
    }

    // Bruit directionnel Von Mises-Fisher

    theta_tot = von_mises(theta_tot, kappa_pred, seed);
    phi_tot   = von_mises(phi_tot, kappa_pred, seed);

    // Mise à jour de la position
    double xn = xpred_mtn[k] + vpred * sin(theta_tot) * cos(phi_tot) * delta_t;
    double yn = ypred_mtn[k] + vpred * sin(theta_tot) * sin(phi_tot) * delta_t;
    double zn = zpred_mtn[k] + vpred * cos(theta_tot)                * delta_t;

    vec3 pos = veriftore(xn, yn, zn);

    vec5 result = {pos.x, pos.y, pos.z, theta_tot, phi_tot};
    return result;
}



/////////////////////////// Fonction principale de simulation //////////////////////////

double simul(int n_proie, int n_pred, int t_max, bool banc_proie, bool banc_pred) {

    unsigned int call_id;
    #pragma omp atomic capture
    call_id = global_call_counter++;

    unsigned int seed = call_id * 2246822519u ^ (call_id * 1664525u + 1013904223u);

    double kappa_proie = kappa_proie_ref / delta_t;
    double kappa_pred  = kappa_pred_ref  / delta_t;

    int Ploc = n_pred;
    int Tloc = t_max;

    int Nb_sous_pas = (int)round(1.0 / delta_t);

    int    capture    = 0;
    double *estomac    = (double*)calloc(Ploc, sizeof(double));
    double *tempsmange = (double*)malloc(Ploc * sizeof(double));
    for (int k = 0; k < Ploc; k++) tempsmange[k] = (double)Tmanip;


    /////////////////// Allocation des tableaux t-1 ///////////////////

    double *xproie_tm1     = (double*)malloc(n_proie * sizeof(double));
    double *yproie_tm1     = (double*)malloc(n_proie * sizeof(double));
    double *zproie_tm1     = (double*)malloc(n_proie * sizeof(double));
    double *thetaproie_tm1 = (double*)malloc(n_proie * sizeof(double));
    double *phiproie_tm1   = (double*)malloc(n_proie * sizeof(double));

    double *xpred_tm1      = (double*)malloc(Ploc * sizeof(double));
    double *ypred_tm1      = (double*)malloc(Ploc * sizeof(double));
    double *zpred_tm1      = (double*)malloc(Ploc * sizeof(double));
    double *thetapred_tm1  = (double*)malloc(Ploc * sizeof(double));
    double *phipred_tm1    = (double*)malloc(Ploc * sizeof(double));

    Cell grid_proie[GRID_D][GRID_H][GRID_W];
    Cell grid_pred [GRID_D][GRID_H][GRID_W];


    //////////////////// Initialisation //////////////////////

    for (int i = 0; i < n_proie; i++) {
        vec3 pos = runif_ctore(torex, torey, torez, &seed);
        dir2 dir = runif_sphere(&seed);
        xproie_tm1[i]     = pos.x;
        yproie_tm1[i]     = pos.y;
        zproie_tm1[i]     = pos.z;
        thetaproie_tm1[i] = dir.theta;
        phiproie_tm1[i]   = dir.phi;
    }

    for (int k = 0; k < Ploc; k++) {
        vec3 pos = runif_ctore(torex, torey, torez, &seed);
        dir2 dir = runif_sphere(&seed);
        xpred_tm1[k]     = pos.x;
        ypred_tm1[k]     = pos.y;
        zpred_tm1[k]     = pos.z;
        thetapred_tm1[k] = dir.theta;
        phipred_tm1[k]   = dir.phi;
    }


    /////////////////// Allocation des buffers mtn/apres ///////////////////

    double *xproie_mtn      = (double*)malloc(n_proie * sizeof(double));
    double *yproie_mtn      = (double*)malloc(n_proie * sizeof(double));
    double *zproie_mtn      = (double*)malloc(n_proie * sizeof(double));
    double *thetaproie_mtn  = (double*)malloc(n_proie * sizeof(double));
    double *phiproie_mtn    = (double*)malloc(n_proie * sizeof(double));

    double *xpred_mtn       = (double*)malloc(Ploc * sizeof(double));
    double *ypred_mtn       = (double*)malloc(Ploc * sizeof(double));
    double *zpred_mtn       = (double*)malloc(Ploc * sizeof(double));
    double *thetapred_mtn   = (double*)malloc(Ploc * sizeof(double));
    double *phipred_mtn     = (double*)malloc(Ploc * sizeof(double));

    double *xproie_apres     = (double*)malloc(n_proie * sizeof(double));
    double *yproie_apres     = (double*)malloc(n_proie * sizeof(double));
    double *zproie_apres     = (double*)malloc(n_proie * sizeof(double));
    double *thetaproie_apres = (double*)malloc(n_proie * sizeof(double));
    double *phiproie_apres   = (double*)malloc(n_proie * sizeof(double));

    double *xpred_apres      = (double*)malloc(Ploc * sizeof(double));
    double *ypred_apres      = (double*)malloc(Ploc * sizeof(double));
    double *zpred_apres      = (double*)malloc(Ploc * sizeof(double));
    double *thetapred_apres  = (double*)malloc(Ploc * sizeof(double));
    double *phipred_apres    = (double*)malloc(Ploc * sizeof(double));

    


    ////////////////////////////// Boucle temporelle //////////////////////////////


    double Rfuite2        = Rfuite  * Rfuite;
    double Rchasse2       = Rchasse * Rchasse;
    double estomac_decay  = pow(decroissance_estomac, delta_t);
    
    for (int t = 1; t < Tloc; t++) {

        bool pred_actifs = (t >= Tstart);

        memcpy(xproie_mtn,     xproie_tm1,     n_proie * sizeof(double));
        memcpy(yproie_mtn,     yproie_tm1,     n_proie * sizeof(double));
        memcpy(zproie_mtn,     zproie_tm1,     n_proie * sizeof(double));
        memcpy(thetaproie_mtn, thetaproie_tm1, n_proie * sizeof(double));
        memcpy(phiproie_mtn,   phiproie_tm1,   n_proie * sizeof(double));

        memcpy(xpred_mtn,     xpred_tm1,     Ploc * sizeof(double));
        memcpy(ypred_mtn,     ypred_tm1,     Ploc * sizeof(double));
        memcpy(zpred_mtn,     zpred_tm1,     Ploc * sizeof(double));
        memcpy(thetapred_mtn, thetapred_tm1, Ploc * sizeof(double));
        memcpy(phipred_mtn,   phipred_tm1,   Ploc * sizeof(double));

        for (int s = 0; s < Nb_sous_pas; s++) {

            build_grid(xproie_mtn, yproie_mtn, zproie_mtn, n_proie, grid_proie);
            build_grid(xpred_mtn,  ypred_mtn,  zpred_mtn,  Ploc,    grid_pred);


            ////////////// Mouvements des proies ///////////////

            for (int i = 0; i < n_proie; i++) {

                if (isnan(xproie_mtn[i])) {
                    xproie_apres[i]     = NAN;
                    yproie_apres[i]     = NAN;
                    zproie_apres[i]     = NAN;
                    thetaproie_apres[i] = NAN;
                    phiproie_apres[i]   = NAN;
                    thetaproie_mtn[i]   = NAN;
                    phiproie_mtn[i]     = NAN;
                    continue;
                }

                vec5 res = deplacement_proie(i, pred_actifs,
                                xproie_mtn, yproie_mtn, zproie_mtn,
                                xpred_mtn,  ypred_mtn,  zpred_mtn,
                                thetaproie_mtn, phiproie_mtn,
                                grid_proie, grid_pred,
                                n_proie, t, banc_proie, kappa_proie, &seed);
                xproie_apres[i]     = res.x;
                yproie_apres[i]     = res.y;
                zproie_apres[i]     = res.z;
                thetaproie_apres[i] = res.theta;
                phiproie_apres[i]   = res.phi;
            }


            ////////////// Mouvements des prédateurs ///////////////

            for (int k = 0; k < Ploc; k++) {
                vec5 res = deplacement_pred(k, pred_actifs,
                                xproie_mtn, yproie_mtn, zproie_mtn,
                                xpred_mtn,  ypred_mtn,  zpred_mtn,
                                thetapred_mtn, phipred_mtn,
                                estomac, tempsmange,
                                grid_proie, grid_pred,
                                n_proie, t, banc_pred, kappa_pred, &seed);
                xpred_apres[k]     = res.x;
                ypred_apres[k]     = res.y;
                zpred_apres[k]     = res.z;
                thetapred_apres[k] = res.theta;
                phipred_apres[k]   = res.phi;
            }


            /////////////// Recensement des captures ///////////////

            if (pred_actifs) {

                build_grid(xproie_apres, yproie_apres, zproie_apres, n_proie, grid_proie);

                for (int k = 0; k < Ploc; k++) {
                    if (tempsmange[k] >= (double)Tmanip && estomac[k] < estomacmax) {

                        int dep0 = cell_dep(zpred_apres[k]);
                        int row0 = cell_row(ypred_apres[k]);
                        int col0 = cell_col(xpred_apres[k]);

                        double best_d2 = Rcapt2;
                        int    best_i  = -1;

                        for (int dd = -1; dd <= 1; dd++) {
                            for (int dr = -1; dr <= 1; dr++) {
                                for (int dc = -1; dc <= 1; dc++) {

                                    int dep = ((dep0 + dd) % GRID_D + GRID_D) % GRID_D;
                                    int r   = ((row0 + dr) % GRID_H + GRID_H) % GRID_H;
                                    int c   = ((col0 + dc) % GRID_W + GRID_W) % GRID_W;
                                    Cell *cell = &grid_proie[dep][r][c];

                                    for (int ci = 0; ci < cell->count; ci++) {
                                        int i = cell->indices[ci];
                                        if (!isnan(xproie_apres[i])) {
                                            vec3 v = dxdy_tore(xpred_apres[k], ypred_apres[k], zpred_apres[k],
                                                               xproie_apres[i], yproie_apres[i], zproie_apres[i]);
                                            double d2 = v.x*v.x + v.y*v.y + v.z*v.z;
                                            if (d2 < best_d2) {
                                                best_d2 = d2;
                                                best_i  = i;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (best_i >= 0) {
                            xproie_apres[best_i] = NAN;
                            yproie_apres[best_i] = NAN;
                            zproie_apres[best_i] = NAN;
                            xproie_mtn[best_i]   = NAN;
                            yproie_mtn[best_i]   = NAN;
                            zproie_mtn[best_i]   = NAN;
                            capture++;
                            estomac[k]   += energie_proie;
                            tempsmange[k] = 0.0;
                        }
                    }
                }

                for (int k = 0; k < Ploc; k++) tempsmange[k] += delta_t;
                for (int k = 0; k < Ploc; k++) estomac[k] *= estomac_decay;
            }


            ///////////////// Rotation des buffers /////////////////

            double *tmp;
            tmp = xproie_mtn;     xproie_mtn     = xproie_apres;     xproie_apres     = tmp;
            tmp = yproie_mtn;     yproie_mtn     = yproie_apres;     yproie_apres     = tmp;
            tmp = zproie_mtn;     zproie_mtn     = zproie_apres;     zproie_apres     = tmp;
            tmp = thetaproie_mtn; thetaproie_mtn = thetaproie_apres; thetaproie_apres = tmp;
            tmp = phiproie_mtn;   phiproie_mtn   = phiproie_apres;   phiproie_apres   = tmp;

            tmp = xpred_mtn;     xpred_mtn     = xpred_apres;     xpred_apres     = tmp;
            tmp = ypred_mtn;     ypred_mtn     = ypred_apres;     ypred_apres     = tmp;
            tmp = zpred_mtn;     zpred_mtn     = zpred_apres;     zpred_apres     = tmp;
            tmp = thetapred_mtn; thetapred_mtn = thetapred_apres; thetapred_apres = tmp;
            tmp = phipred_mtn;   phipred_mtn   = phipred_apres;   phipred_apres   = tmp;
        }

        memcpy(xproie_tm1,     xproie_mtn,     n_proie * sizeof(double));
        memcpy(yproie_tm1,     yproie_mtn,     n_proie * sizeof(double));
        memcpy(zproie_tm1,     zproie_mtn,     n_proie * sizeof(double));
        memcpy(thetaproie_tm1, thetaproie_mtn, n_proie * sizeof(double));
        memcpy(phiproie_tm1,   phiproie_mtn,   n_proie * sizeof(double));

        memcpy(xpred_tm1,     xpred_mtn,     Ploc * sizeof(double));
        memcpy(ypred_tm1,     ypred_mtn,     Ploc * sizeof(double));
        memcpy(zpred_tm1,     zpred_mtn,     Ploc * sizeof(double));
        memcpy(thetapred_tm1, thetapred_mtn, Ploc * sizeof(double));
        memcpy(phipred_tm1,   phipred_mtn,   Ploc * sizeof(double));
    }


    /////////////////////////// Libération mémoire ///////////////////////////

    free(xproie_mtn);      free(yproie_mtn);      free(zproie_mtn);
    free(thetaproie_mtn);  free(phiproie_mtn);
    free(xproie_apres);    free(yproie_apres);    free(zproie_apres);
    free(thetaproie_apres);free(phiproie_apres);
    free(xproie_tm1);      free(yproie_tm1);      free(zproie_tm1);
    free(thetaproie_tm1);  free(phiproie_tm1);

    free(xpred_mtn);       free(ypred_mtn);       free(zpred_mtn);
    free(thetapred_mtn);   free(phipred_mtn);
    free(xpred_apres);     free(ypred_apres);     free(zpred_apres);
    free(thetapred_apres); free(phipred_apres);
    free(xpred_tm1);       free(ypred_tm1);       free(zpred_tm1);
    free(thetapred_tm1);   free(phipred_tm1);

    free(estomac);
    free(tempsmange);

    return (double)capture / ((double)Ploc * (double)(Tloc - Tstart));
}




int rep_fonctionnelle(int nmax, int n_pred, int t_max, bool banc_proie, bool banc_pred, const char *filename) {
    if (nmax <= 0 || n_pred <= 0 || t_max <= 0 || filename == NULL) {
        fprintf(stderr, "rep_fonctionnelle: arguments invalides\n");
        return -1;
    }

    int nb_rep = 10;
    int total_taches = nmax * nb_rep;

    FILE *f = fopen(filename, "w");
    if (!f) { perror("fopen"); return -1; }
    fprintf(f, "# N capture\n");

    double *sommes    = (double*)calloc(nmax + 1, sizeof(double));
    double *resultats = (double*)malloc((nmax + 1) * sizeof(double));
    if (!sommes || !resultats) {
        fprintf(stderr, "rep_fonctionnelle: allocation échouée\n");
        free(sommes); free(resultats); fclose(f);
        return -1;
    }

    int taches_terminees = 0;

    #pragma omp parallel for schedule(dynamic, 1)
    for (int task = 0; task < total_taches; task++) {
        int N = (task / nb_rep) + 1;

        double val = simul(N, n_pred, t_max, banc_proie, banc_pred);

        #pragma omp atomic
        sommes[N] += val;

        #pragma omp atomic
        taches_terminees++;

        // N est terminé seulement quand toutes ses répétitions sont faites

        #pragma omp critical
        {
            int N_termines = taches_terminees / nb_rep;
            printf("\rSimulation en cours : N = %d / %d (%.1f%%)",
                   N_termines, nmax,
                   100.0 * N_termines / nmax);
            fflush(stdout);
        }
    }

    for (int N = 1; N <= nmax; N++) {
        resultats[N] = sommes[N] / nb_rep;
        fprintf(f, "%d %.6f\n", N, resultats[N]);
    }

    free(sommes);
    free(resultats);
    printf("\nSimulations terminées.\n");
    fclose(f);
    return 0;
}


int main(void) {

    printf("=== Tous solitaires ===\n");
    if (rep_fonctionnelle(Nmax, P, T, false, false, "solitaires.txt") != 0)
        fprintf(stderr, "Erreur lors de la création de solitaires.txt\n");

    printf("=== Proies en banc, prédateurs solitaires ===\n");
    if (rep_fonctionnelle(Nmax, P, T, true, false, "proie_banc.txt") != 0)
        fprintf(stderr, "Erreur lors de la création de proie_banc.txt\n");

    printf("=== Proies solitaires, prédateurs en banc ===\n");
    if (rep_fonctionnelle(Nmax, P, T, false, true, "pred_banc.txt") != 0)
        fprintf(stderr, "Erreur lors de la création de pred_banc.txt\n");

    printf("=== Proies + prédateurs en banc ===\n");
    if (rep_fonctionnelle(Nmax, P, T, true, true, "banc_banc.txt") != 0)
        fprintf(stderr, "Erreur lors de la création de banc_banc.txt\n");

    return 0;
}



// A rentrer dans le terminale pour lancer le code : 

// 1) gcc -O3 -march=native -fopenmp -flto -fno-math-errno -fno-trapping-math -o simulation opti_quadruple_3D.c -lm


// 2) ./simulation