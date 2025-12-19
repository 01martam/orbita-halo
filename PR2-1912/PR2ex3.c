#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Newton_b.h"
#include "TestFunction.h"
#include "QRsolve_b.h"

/* ARRELS en coordenades reals:
   (1,0), (-1/2, sqrt(3)/2), (-1/2, -sqrt(3)/2)
*/
#define MIN(x,y) ((x)<(y)? (x):(y))

double dist(double *x, double *r)
{
  return hypot(x[0] - r[0], x[1] - r[1]);
  //return sqrt((x[0] - r[0]) * (x[0] - r[0]) + (x[1] - r[1]) * (x[1] - r[1]));
}


int main(int argc, char *argv[]) {
    int maxit;
    double tol;
    if(argc<3
         ||sscanf(argv[1],"%lf",&tol)!=1|| tol <= 0
         ||sscanf(argv[2],"%d",&maxit)!=1|| maxit <= 0){
         fprintf(stderr,"%s tol maxit \n",argv[0]);
      return -1;
    }
    /*Fi consola*/
    
    // Arrels del polinomi
    double sqrt3_2 = sqrt(3.0) / 2.0;
    double arrel1[2] = {1.0, 0.0};
    double arrel2[2] = {-0.5, sqrt3_2};
    double arrel3[2] = {-0.5, -sqrt3_2};
    
    // Obrir arxius de sortida
    FILE *g1 = fopen("g1.txt", "w");
    if (g1==NULL){
      fprintf(stderr, "ERROR: No es pot obrir g1.txt\n");
      return -1;
    }
    FILE *g2 = fopen("g2.txt", "w");
    if (g2==NULL){
      fprintf(stderr, "ERROR: No es pot obrir g2.txt\n");
      fclose(g1);
      return -1;
    }
    FILE *g3 = fopen("g3.txt", "w");
    if (g3==NULL){
      fprintf(stderr, "ERROR: No es pot obrir g3.txt\n");
      fclose(g1);fclose(g2);
      return -1;
    }
    //malla [-1,1] x [-1,1]
    int N=500;// Resolució de la malla
    int i,j,punts1=0,punts2=0,punts3=0;
    double x[2],z[2];
    double step= 2.0 / (N - 1);
    

    // Bucle principal sobre la malla
    for(i=0;i<N;i++){
      for (j=0;j<N;j++){
        //Punt inicial de la malla
        x[0]=-1.0+i*step;
        x[1]=-1.0+j*step;
        
        //Copiem aixo a un nou vector per no modificarlo.
        //double z[2]={x[0],x[1]};
        z[0]=x[0];
        z[1]=x[1];
        
        //Apliquem Newton
        newton(2,2,z,funcio_fractal,maxit,tol,0,NULL);
        //Distancies a les arrels
        double d1=dist(z,arrel1);
        double d2=dist(z,arrel2);
        double d3=dist(z,arrel3);
        double minim=MIN(d1,MIN(d2,d3));
        
        //Guardem el punt inicial
        
        if(minim==d1){
          fprintf(g1, "%lf %lf\n", x[0], x[1]);
          punts1++;
        }
        else if(minim==d2){
          fprintf(g2, "%lf %lf\n", x[0], x[1]);
          punts2++;
        }
        else{
          fprintf(g3, "%lf %lf\n", x[0], x[1]);
          punts3++;
        }
        
      }
    }

    // Tancar fitxers
    fclose(g1); fclose(g2); fclose(g3);

    // Estadístiques
    printf("\n=== Resultats del fractal complet ===\n");
    printf("Punts que convergeixen a l'arrel (z=1): %d\n", punts1);
    printf("Punts que convergeixen a l'arrel (z=-0.5+i√3/2): %d\n", punts2);
    printf("Punts que convergeixen a l'arrel (z=-0.5-i√3/2): %d\n", punts3);
    printf("Total de punts procesats: %d\n", N*N);

    printf("\nFitxers generats: g1.txt, g2.txt, g3.txt\n");
    //gnuplot> plot \"g1.txt\" u 1:2 w dots lt rgb \"red\"
                  //\"g2.txt\" u 1:2 w dots lt rgb \"blue\"
                  //\"g3.txt\" u 1:2 w dots lt rgb \"green\"
                  
    return 0;
}
