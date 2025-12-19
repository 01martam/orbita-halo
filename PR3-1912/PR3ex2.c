#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rk45_b.h"

int camp(int n, double t, double *x, double *f, double *prm){
  if (n!=2){return -1;}
  double a=prm[0];
  f[0]=a*x[0]-x[1];
  f[1]=x[0]+a*x[1];
  return 0;
}

int funcio_exacta(double t,double *f, double *x0,double *prm){
  double a=prm[0];
  double expo=exp(a*t), cost=cos(t), sint=sin(t);//evitar calcular la mateixa cosa diversos cops.
  f[0]=expo*(x0[0]*cost-x0[1]*sint);
  f[1]=expo*(x0[0]*sint+x0[1]*cost);
  return 0;
}

int main(int argc, char *argv[]){
  /*linea de comandes */
  double t=0, tf,h,hmax,hmin,prm[1],x_inicial[2];
  int miss;//tf=T
  if(argc<8||sscanf(argv[1],"%lf",&prm[0])!=1
       ||sscanf(argv[2],"%lf",&tf)!=1
       ||sscanf(argv[3],"%lf",&x_inicial[0])!=1
       ||sscanf(argv[4],"%lf",&x_inicial[1])!=1
       ||sscanf(argv[5],"%lf",&h)!=1
       ||sscanf(argv[6],"%lf",&hmin)!=1
       ||sscanf(argv[7],"%lf",&hmax)!=1
       ||sscanf(argv[8],"%d",&miss)!=1)
  {
    fprintf(stderr,"%s a T x0 y0 h hmin hmax missatge(si:1 , No:0)\n",argv[0]);
    return -1;
  }

  /*Fi de linia de comandes*/
  FILE *archiu= fopen("arxiu.txt","w");
  if(archiu==NULL){
    fprintf(stderr,"No s'ha pogut obrir arxiu.txt\n");
    return -1;
  }
  double x[2],x_exac[2];
  t=0;
  x[0]=x_inicial[0];//x0
  x[1]=x_inicial[1];//y0
  if(tf<0){h=-fabs(h);}//asegurar-nos de que fem el pas cap enrere.
  printf("t\t x_exacta\t y_exacta\t x_aprox\t y_aprox\t diff_x\t diff_y\n");
  funcio_exacta(t,x_exac,x_inicial,prm);
  
  fprintf(archiu,"%.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.2e\t %.2e\n",t,x_exac[0],x_exac[1],
            x[0],x[1],fabs(x[0]-x_exac[0]),fabs(x[1]-x_exac[1]));

  while(fabs(t)<fabs(tf)){

    if(fabs(t+h)>fabs(tf)){
      h=tf-t;//asegurarnos que no sortim del interval
    }

    //fer un pas
    if(rk45(&t,x,&h,hmin,hmax,1e-12,2,camp,prm,miss)!=0){
      fprintf(stderr, "Error en rk45\n");
      fclose(archiu);
      return -1;
    }
    funcio_exacta(t,x_exac,x_inicial,prm);
    fprintf(archiu,"%.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.2e\t %.2e\n",t,x_exac[0],x_exac[1],
            x[0],x[1],fabs(x[0]-x_exac[0]),fabs(x[1]-x_exac[1]));

  }
  fclose(archiu);
  return 0;
}
