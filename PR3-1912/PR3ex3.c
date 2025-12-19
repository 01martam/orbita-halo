#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rk45_b.h"

int EDO_lorenz(int n, double t, double *x, double *f, double *prm){
  if (n!=3){
        fprintf(stderr,"ERROR en la dimensió (n=3).\n");
        return -1;}
  double sigma=prm[0];
  double rho=prm[1];
  double beta=prm[2];
  //x=x[0] y=x[1] z=x[2]
  f[0]=sigma*(x[1]-x[0]);
  f[1]=rho*x[0]-x[1]-x[0]*x[2];
  f[2]=-beta*x[2]+x[0]*x[1];
  return 0;
}


int main(int argc, char *argv[]){
  int n=3;
  double x[n],prm[3]={10.0,28.0,8./3.};//parametre sigma, rho, beta
  double T=50,t=0,hmin,hmax,h,tol;

  
  /*linea de comandes */
  if(argc<8||sscanf(argv[1],"%lf",&tol)!=1||tol<=0
       ||sscanf(argv[2],"%lf",&x[0])!=1
       ||sscanf(argv[3],"%lf",&x[1])!=1
       ||sscanf(argv[4],"%lf",&x[2])!=1
       ||sscanf(argv[5],"%lf",&h)!=1
       ||sscanf(argv[6],"%lf",&hmin)!=1
       ||sscanf(argv[7],"%lf",&hmax)!=1)
  {
    fprintf(stderr,"%s tol x0 y0 z0 h hmin hmax \n",argv[0]);
    return -1;
  }
  
  /*Fi de linia de comandes*/

   FILE *F1 = fopen("Lorenz.txt","w");
    if(F1==NULL){
        fprintf(stderr,"No s'ha pogut obrir Lorenz.txt\n");
        return -1;
    }
  printf("t\t x\t y\t z\n");
  fprintf(F1,"%.6f\t %.6f\t %.6f\t %.6f\n",t,x[0],x[1],x[2]);
   while(fabs(t)<fabs(T)){fprintf(F1,"%.6f\t %.6f\t %.6f\t %.6f\n",t,x[0],x[1],x[2]);
        if(fabs(t+h)>fabs(T)){
            h=T-t;//asegurarnos que no sortim del interval
        }
        //Fer un pas del rk45
        if(rk45(&t,x,&h,hmin,hmax,tol,n,EDO_lorenz,prm,0)!=0){
            fprintf(stderr, "Error en rk45 \n");
            fclose(F1);
            return -1;
        }
        fprintf(F1,"%.6f\t %.6f\t %.6f\t %.6f\n",t,x[0],x[1],x[2]);

   }
    fclose(F1);
 return 0;
}
