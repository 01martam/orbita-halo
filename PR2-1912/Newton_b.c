#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QRsolve_b.h"

/* Funció que implementa el mètode de newton en diverses variables
m:nombre de components de F
n: nombre de components de la solució
x:En entrar x[i=0,..n-1] ha de ser la llavor x_0
  En sortir ha de ser la solució x_k proposada pel mètode de newton
fdf:Funció que calcula F i DF
maxit: el nombre maxim d'iterats que farà el mètode
tol:Tolerancia per decidir si el mètode troba solució o no.
prm: parametres necesaris en la funcio fdf()
passos: Imprimir calculs intermitjos
    si passos ==0 llavors no s'imprimeix res.
    Si passos ==1 llavors s'imprimeix norma de f,norma de y
    Si passos ==2 llavors s'imprimeix x,f,df,y
La funcio retorna
    (el nombre de iterats si)>0 tot va bé
    -1 si no hi ha convergencia
    -2 Si hi ha hagut algun problema en resoldre algun QRsolve
    -3 si no s'ha pogut avaluar la funcio fdf
    -4 si no s'ha pogut asignar memòria
*/
double norma_2(int m, double *f){
  double norm=0;
  for(int i=0; i<m;i++){
    norm=norm+f[i]*f[i];
  }
  return sqrt(norm);
}


int newton(int m, int n, double *x, int (*fdf)(int m,int n,double *x,double *f,double *df, void *prm),int maxit, double tol,int passos, void *prm){
  //Fixemnos que el que hi ha dintre de fdf no l'estem donan nombre si no que diem quantes variables crida i quin tipus, perfectament podriem treure els noms.
  int iterats,i;
  double norm_f,norm_y,*y,*f,*df,*tau;
  //Reserva memoria
  f=(double *)malloc(m*sizeof(double));
  df=(double *)malloc(m*n*sizeof(double));
  tau=(double *)malloc(n*sizeof(double));
  y=(double *)malloc(n*sizeof(double));
  
  if (y == NULL || f == NULL || df == NULL || tau == NULL){
    fprintf(stderr,"Error en asignar memoria.\n");
    free(y); free(f); free(df); free(tau);
    return -4;
  }
  
  for (iterats=0;iterats<maxit;iterats++){
    //Avaluar F(x_k) i DF(x_k)
    if(fdf(m,n,x,f,df,prm)!=0){
      fprintf(stderr,"Error en la avaluació de la funció en l'iterat%d\n",iterats);
      free(y); free(f); free(df); free(tau);
      return -3;
    }
    //Calcular la norma de F(x)
    norm_f=norma_2(m,f);
    
    //Verifiquem la convergencia
    if(norm_f<tol){
      if(passos==1){
        printf("it %d norm_f %e \n",iterats,norm_f);
      }
      iterats++;//+1 si contem 0,1,2,3,4,...
      printf("Newton convergeix despres de %d iteracions\n",iterats);
      free(y); free(f); free(df); free(tau);
      
      return iterats;
    }
    //Resoldre el sistema DF(x_k)y_k=F(x_k)
    if (QRsolve(m, n, df, f, tau, y, tol,0) != 0) {
      fprintf(stderr, "Error resolent el sistema en l'iterat %d\n", iterats);
      free(y); free(f); free(df); free(tau);
      return -2;
    }
    if(passos==2){
      printf("\nIteracio %d\n",iterats);
      printf("x_{%d}=",iterats);
      imprimir_matriu(1,n,x);
      printf("f(x_{%d})=",iterats);
      imprimir_matriu(1,m,f);
      printf("Df(x_{%d})=",iterats);
      imprimir_matriu(m,n,df);
      printf("y_{%d}=",iterats);
      imprimir_matriu(1,n,y);
    }
    if(passos==1){
      //Calcular la norma de y
      norm_y=norma_2(n,y);
      printf("it %d norm_f %e norm_y %e\n",iterats,norm_f,norm_y);
    }
    
    //Actualitza el valor de x
    for(i=0;i<n;i++){
      x[i]=x[i]-y[i];
    }
  }
  //Newton no convergeix
  fprintf(stderr,"Error : Newton no ha convergit amb %d iteracions.\n",maxit);
  free(y); free(f); free(df); free(tau);
  return -1;
}
