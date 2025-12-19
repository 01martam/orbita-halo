#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Newton_b.h"
#include "TestFunction.h"
#include "QRsolve_b.h"


int main(int argc, char *argv[]){
    int m=3,n=3,passos,maxit;
    double x[n],tol;
    //tol=1e-5;
    //maxit=1e3;
 /*
 * Línia de comandes:
  double tol double maxit int passos double x0[n]
 */
     if(argc<7
          ||sscanf(argv[1],"%lf",&tol)!=1|| tol <= 0
          ||sscanf(argv[2],"%d",&maxit)!=1|| maxit <= 0
          ||sscanf(argv[3],"%d",&passos)!=1
          ||sscanf(argv[4],"%lf",&x[0])!=1
          ||sscanf(argv[5],"%lf",&x[1])!=1
          ||sscanf(argv[6],"%lf",&x[2])!=1){
          
          fprintf(stderr,"%s tol maxit passos x0[0] x0[1] x0[2] \n",argv[0]);
          return -1;
     }
/*
 * Fi línia de comandes
 */




    int result=newton(m,n,x,function1,maxit,tol,passos,NULL);
    if (result>=0){
        printf("Solucio trobada:\n x=");
        imprimir_matriu(1,n,x);
        double f[m],df[m*n];
        printf("Si x3=pi/2+pi*k , k qualsevol enter on pi/2 es aproximadament %.8lf. \n",M_PI/2);
        printf("solucions analitiques són: (0 0 x3) o (-1/x3 1/x3^2 x3)\n");
        function1(m,n,x,f,df,NULL);
        printf("F(x)=");
        imprimir_matriu(1,m,f);
        printf("DF=");
        imprimir_matriu(m,n,df);
        return 0;
    }
    else if(result==-1){
        printf("El metode no convergeix.\n");
    }
    else {
        fprintf(stderr, "Error en el codi de Newton\n");
        return -1;
    }

}
