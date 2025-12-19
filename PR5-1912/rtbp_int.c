// rtbp_int.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "flux_b.h"
#include "rtbp.h"

int main() {
    double prm[1] = {3.040357143e-6};//mu
    double T = 3.05177804298575;
    double x[6] = {-0.9889514941464619, 0, 0.003249420704787417,0, -0.9992369544112847, 0};

    int m = 3;  // RTBP: 3 graus de libertat
    double t = 0.0;
    int n = 2*m;
    
    // Parámetros integración
    double h = 0.01;
    double hmin = 1e-10;
    double hmax = 0.1;
    double tol = 1e-13;
    int npasmax = 100000;
    
    if(flux(&t,x,&h,T,hmin,hmax,tol,npasmax,n,rtbp,prm,1)!=0){
        return -1;
    }
    //llamar rtbp_int.exe >orb.txt
    return 0;
}
