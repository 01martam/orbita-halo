#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "rk45_b.h"




void imprimir_vector(int n,double *a) {
      //printf("(");
        for (int j = 0; j < n; j++) {
            printf("%.15lf\t ", a[j]);
        }
        printf("\n");
}

/*
 * hmin, hmax, tol són resp. el pas mı́nim, pas màxim i la tolerància per passara rk45()
 * El paràmetre npasmax és el nombre màxim de crides a rk45() que permetem.
 * Els paràmetre n és el nombre de variables dependents (dimensió del sistema d’EDO)
 * camp és un punter a funció que implementa les equacions que volem integrar, i prm són paràmetres per passar al camp.
 * miss si es 1 imprimeix cada valor t i x calculat entremitg.
 *
 *Retorna
 *0 si ha acabat amb èxit,
 *-1 si rk45() ha tornat un error,
 * -2 si s’han superat npasmax passos d’integració numèrica.
 *




 */
int flux2(double *t, double *x, double *h, double T, double hmin,
          double hmax, double tol, int npasmax, int n,
          int (*camp)(int n, double t, double *x, double *f, double *prm),
          double *prm,int miss)
{
    //Assegurar que h té el signe correcte respecte T.
    *h = (T > 0) ? fabs(*h) : -fabs(*h);

    double tf=*t + T;
    int pasos = 0;
    if(T>0){
      while (*t<tf) {
          if (miss==1){
              printf("%.6lf\t",*t);
              imprimir_vector(n,x);
          }
  
          if(*t+*h>tf){
              *h=tf-*t;//asegurarnos de no sortir del interval
          }
  
  
          if (rk45(t, x, h, hmin, hmax,tol, n, camp, prm,0)!= 0) {
              fprintf(stderr, "Error intern en rk45()\n");
              return -1;
          }
          //Se ha integrat per tant podem contar un pas.
          pasos++;
  
          if (pasos > npasmax) {
              fprintf(stderr, "Superat el nombre de passos permessos\n");
              return -2;  // Número máximo de pasos excedido
          }
      }
    return 0;
    }


    else if(T<0){
      while (*t>tf) {
          if (miss==1){
              printf("%.6lf\t",*t);
              imprimir_vector(n,x);
          }
  
          if(*t+*h<tf){
              *h=tf-*t;//asegurarnos de no sortir del interval
          }
  
  
          if (rk45(t, x, h, hmin, hmax,tol, n, camp, prm,0)!= 0) {
              fprintf(stderr, "Error intern en rk45()\n");
              return -1;
          }
          //Se ha integrat per tant podem contar un pas.
          pasos++;
  
          if (pasos > npasmax) {
              fprintf(stderr, "Superat el nombre de passos permessos\n");
              return -2;  // Número máximo de pasos excedido
          }
      }
    return 0;
    }
    return 6666;
}

int flux(double *t, double *x, double *h, double T, double hmin,
          double hmax, double tol, int npasmax, int n,
          int (*camp)(int n, double t, double *x, double *f, double *prm),
          double *prm,int miss)
{
    //Assegurar que h té el signe correcte respecte T.
    *h = (T > 0) ? fabs(*h) : -fabs(*h);

    double tf=*t + T;
    int pasos = 0;
    while ((T>0&& *t<tf)||(T<0&& *t>tf)) {
        if (miss==1){
            printf("%.6lf\t",*t);
            imprimir_vector(n,x);
        }

        if((T>0&&*t+*h>tf)||(T<0&&*t+*h<tf)){
            *h=tf-*t;//asegurarnos de no sortir del interval
        }

        if (rk45(t, x, h, hmin, hmax,tol, n, camp, prm,0)!= 0) {
            fprintf(stderr, "Error intern en rk45()\n");
            return -1;
        }
        //Se ha integrat per tant podem contar un pas.
        pasos++;

        if (pasos > npasmax) {
            fprintf(stderr, "Superat el nombre de passos permessos\n");
            return -2;  // Número máximo de pasos excedido
        }
    }
    return 0;

}
