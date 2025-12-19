// opham_main_pendol.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "opham_b.h"
#include "pendol.h"
#include "QRsolve_b.h"
#include "Newton_b.h"
#include "rk45_b.h"
#include "flux_b.h"

int main(int argc, char *argv[]) {
    /* Consola*/
    double tol, cg[3];
    // Coeficientes del hiperplano: g(x) = p = 0 -> cg[3] = {0.0, 1.0, 0.0};
    int maxit;
    if(argc<6
       ||sscanf(argv[1],"%lf",&tol)!=1
       ||sscanf(argv[2],"%d",&maxit)!=1
       ||sscanf(argv[3],"%lf",&cg[0])!=1
       ||sscanf(argv[4],"%lf",&cg[1])!=1
       ||sscanf(argv[5],"%lf",&cg[2])!=1){
       fprintf(stderr,"%s tol maxit cg0 cg1 cg2 \n",argv[0]);
    return -1;
  }
 /*Fi consola*/
 
    int m = 1;  // péndul: 1 grau de llibertat
    //int n = 2*m; // n = 2
    //int dimF=2+n, dimX=1+n;//dimF=4, dimX=3
    // Parámetres
    double hh = -cos(0.5);  // energía h = H(0.5, 0)
    //pendol_h(xx+1,&hh,NULL); // xx+1 = x0
    double xx[3]={2 * M_PI, 0.5, 0.0};
    //xx=(T,x) R^(n+1)
    //double T=xx[0];
    //double *x0=xx+1;//x0 te n components
 
    
    // Parámetres per flux()
    double h0rk = 0.01;
    double hminrk = 1e-10;
    double hmaxrk = 0.1;
    double tolrk = 1e-13;
    int npasmxrk = 10000;  

    printf("Energia h = %f\n", hh);
    printf("Aproximacio inicial: T = %f, x = (%f, %f)\n", xx[0], xx[1], xx[2]);

    int result = opham(m, hh, xx, tol, maxit, cg, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, pendol_h, pendol, NULL);

    if (result >= 0) {
        printf("Periode T = %.15f\n", xx[0]);
        printf("Condició inicial: x = (%.15f, %.15f) \n", xx[1], xx[2]);
        printf("T0-Tfin=%e\n", fabs(2 * M_PI-xx[0]));
        printf("q0-qfin=%e\n", fabs(0.5-xx[1]));
        printf("p0-pfin=%e\n", fabs(0.0-xx[2]));
        double H;
        pendol_h(xx+1,&H,NULL); 
        printf("H(x)-h=%e\n", H-hh);
        printf("g(x)=%e\n", cg[0]*xx[1]+cg[1]*xx[2]+cg[2]);
    } else {
        printf("\nNo convergencia despres de %d iteracions\n", maxit);
    }

    return 0;
}
