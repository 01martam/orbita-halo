#include <stdio.h>
#include <math.h>
#include "QRsolve_b.h"
#include "flux_b.h"
#include "rk45_b.h"


// #define A(i, j) a[(i) + (j) * m]
// void imprimir_matriu2(int m, int n,double *a) {
//   for (int i = 0; i < m; i++) {
//     printf("|");
//     for (int j = 0; j < n; j++) {
//       printf("%g ", A(i, j));
//     }
//     printf("|\n");
//   }
// }
// #undef A

// Sistema: ẋ = α(1-r²)x - y, ẏ = x + α(1-r²)y
int EDO(int n, double t, double *x, double *f, double *prm){
  double alpha=prm[0],r2;
  if (n!=2&&n!=6){
    fprintf(stderr,"ERROR en la dimensió (n=2 per flux, n=6 per la variacionals).\n");
    return -1;
  }
  else if (n==2){
    r2=x[0]*x[0]+x[1]*x[1];
    f[0]=alpha*(1-r2)*x[0]-x[1];
    f[1]=x[0]+alpha*(1-r2)*x[1];
    return 0;
  }
  else{
    n=2;
    #define A(i,j) x[n+(j)*n+(i)]
    r2=x[0]*x[0]+x[1]*x[1];
    f[0]=alpha*(1-r2)*x[0]-x[1];
    f[1]=x[0]+alpha*(1-r2)*x[1];

    //Calculem Jacobia
    double J[2][2];
    J[0][0]=alpha*(-2*x[0]*x[0]+(1-r2));//df0/dx
    J[1][0]=1-2*alpha*x[0]*x[1];//df1/dx
    J[0][1]=-1-2*alpha*x[0]*x[1];//df0/dy
    J[1][1]=alpha*(-2*x[1]*x[1]+(1-r2));//df1/dy
    //variacional
    f[2]=J[0][0]*A(0,0)+J[0][1]*A(1,0);
    f[3]=J[1][0]*A(0,0)+J[1][1]*A(1,0);
    f[4]=J[0][0]*A(0,1)+J[0][1]*A(1,1);
    f[5]=J[1][0]*A(0,1)+J[1][1]*A(1,1);
    #undef A
    return 0;

  }

}

//Versio 1
// Función para calcular Dϕ usando diferencias finitas
int derivada_numerica(double t0, double *x0, double T, double delta,double h,double hmin,double hmax, double tol, int npasmax,
                      int (*camp)(int n, double t, double *x, double *f, double *prm),
                      double *prm, double *Dphi) {
  
    double xp[2],xm[2],t,haux;
    for (int j = 0; j < 2; j++) {
        // Perturbació
        for (int i = 0; i < 2; i++) {
            xp[i] = x0[i];
            xm[i] = x0[i];
        }
        xp[j] += delta;//x0 + δ·e_j
        xm[j] -= delta;//x0 - δ·e_j

        // Integrar flux per x0 + δ·e_j
        t = t0;
        haux = h;
        
        if(flux(&t, xp, &haux, T, hmin, hmax, tol, npasmax, 2, camp, prm,0)<0){
          return -1;
        }

        // Integrar flux per x0 - δ·e_j
        t = t0;
        haux = h; 
        
        if(flux(&t, xm, &haux, T, hmin, hmax, tol, npasmax, 2, camp, prm,0)<0){
          return -1;
        }
        
        // Diferencia finita centrada per ∂phi/∂x_j
        for (int i = 0; i < 2; i++) {
          Dphi[(i) + (j)*2] = (xp[i] - xm[i]) / (2.0 * delta);
        }
    }
    return 0;
}



// Función para calcular Dϕ usando ecuaciones variacionales
int derivada_variacional(double t0, double *x0, double T,double h,double hmin,double hmax, double tol, int npasmax,
                         int (*camp)(int n, double t, double *x, double *f, double *prm),
                         double *prm, double *Dphi,double *x) 
{
   //x + A (matriz n×n)
    //Concicions inicials
    x[0]=x0[0];
    x[1]=x0[1];
    
    x[2]=1.0;//A00
    x[3]=0.0;//A10
    x[4]=0.0;//A01
    x[5]=1.0;//A11
    
    // Integrar el sistema extendido
    double t = t0;

    int r = flux(&t, x, &h, T, hmin, hmax,tol, npasmax, 6, camp, prm,0);
    if (r<0){
      return r;
    }
    //Extreure matriu Dphi
    Dphi[0] = x[2];  //A00 
    Dphi[1] = x[3];  // A01 
    Dphi[2] = x[4];  // A10 
    Dphi[3] = x[5];  // A11 
    return 0;
}

//versio 2

#define A(i,j) x[n+(j)*n+(i)]
#define K(i,j) Dphi[(j)*n+(i)]

int dernum(int n,double t0, double *x0, double T, double delta,double h,double hmin,double hmax, double tol, int npasmax,
           int (*camp)(int n, double t, double *x, double *f, double *prm),
           double *prm, double *Dphi) {
  
  double xp[n],xm[n];
  double t,haux;
  for (int j = 0; j < n; j++) {
    // Perturbació
    for (int i = 0; i < n; i++) {
      xp[i] = x0[i];
      xm[i] = x0[i];
    }
    xp[j] += delta;//x0 + delta·e_j
    xm[j] -= delta;//x0 - delta·e_j
    
    // Integrar flux per x0 + delta·e_j
    t = t0;
    haux = h;
    
    if(flux(&t, xp, &haux, T, hmin, hmax, tol, npasmax, n, camp, prm,0)<0){
      return -2;
    }
    
    // Integrar flux per x0 - delta·e_j
    t = t0;
    haux = h; 
    
    if(flux(&t, xm, &haux, T, hmin, hmax, tol, npasmax, n, camp, prm,0)<0){
      return -2;
    }
    
    // Diferencia finita centrada per ∂phi/∂x_j
    for (int i = 0; i < n; i++) {
      K(i,j) = (xp[i] - xm[i]) / (2.0 * delta);
    }
  }
  return 0;
}

int dervar(int n, double t0, double *x0, double T,double h,double hmin,double hmax, double tol, int npasmax,
           int (*camp)(int n, double t, double *x, double *f, double *prm),
           double *prm, double *Dphi,double *x) 
{
  int ntotal=n*(n+1);
  double t = t0;
  //x + A (matriz n×n)
  //Concicions inicials
  for (int i=0;i<n;i++){
    x[i]=x0[i];
  }
  for (int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A(i,j)=(i==j)?1.0:0.0;
    }
  }
  
  
  // Integrar el sistema extendido
  int r = flux(&t, x, &h, T, hmin, hmax,tol, npasmax, ntotal, camp, prm,0);
  if (r<0){
    return r;
  }
  for (int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      K(i,j)=A(i,j);
    }
  }
  return 0;
}

#undef A
#undef K


int main() {
    double prm[1] = {0.5}, x0[2] = {0.5, 0.3},X[6];//alpha
    int npasmax=100000;
    double t0 = 0.0,T = 0.5;
    double h=1e-3,hmin=1e-8,hmax=1.0,tol=1e-10,delta = 1e-6;
    double Dnum[4];  //  diferencies finites
    double Dvar[4];  // per equacions variacionals

    // 1. Calcular Dphi amb variacionals
    if(derivada_variacional(t0, x0, T,h,hmin,hmax,tol,npasmax, EDO, prm, Dvar,X)!=0){return -1;}
    
    printf("phi_T(x0) = (%.15f, %.15f)\n", X[0], X[1]);
    imprimir_matriu(1, 6, X);
    
    printf("\n Dphi(T) per variacionals:\n");
    imprimir_matriu(2, 2, Dvar);
    
    // 1. Calcular Dphi amb variacionals
    if(dervar(2,t0, x0, T,h,hmin,hmax,tol,npasmax, EDO, prm, Dvar,X)!=0){return -1;}

    printf("\n Dphi(T) per variacionals:\n");
    imprimir_matriu(2, 2, Dvar);
    
    // 2. Calcular Dphi amb diferencies finites
    if (derivada_numerica(t0, x0, T, delta, h, hmin, hmax,tol,npasmax,EDO,prm, Dnum)!=0){return -1;}
    printf("\nDphi(T) per diferencies finites(ordre 1):\n");
    imprimir_matriu(2, 2, Dnum);
    
    // 2. Calcular Dphi amb diferencies finites
    if (dernum(2,t0, x0, T, delta, h, hmin, hmax,tol,npasmax,EDO,prm, Dnum)!=0){return -1;}
    printf("\nDphi(T) per diferencies finites(ordre 1):\n");
    imprimir_matriu(2, 2, Dnum);
    
    return 0;
}
