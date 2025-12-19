#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "opham_b.h"
#include "pendol.h"
#include "QRsolve_b.h"
#include "Newton_b.h"
#include "rk45_b.h"
#include "flux_b.h"


#define DF(i,j) df_num[(i) + dimF*(j)]
int main() {
  int m = 1;  // péndul: 1 grau de llibertat
  int n = 2*m; // n = 2
  int dimF=2+n, dimX=1+n;//dimF=4, dimX=3
  // Parámetres
  double hh = -cos(0.5);  // energía h = H(0.5, 0)
  //pendol_h(xx+1,&hh,NULL); // xx+1 = x0
  double xx[3]={0.1,0.5,0.0};
  //xx=(T,x) R^(n+1)
  //double T=xx[0];
  //double *x0=xx+1;//x0 te n components
  double cg[3]={0.0,1.0,0.0};// Coeficients del hiperplà g(x) = p 
  
  // Parámetres per flux()
  double h0rk = 0.01;
  double hminrk = 1e-10;
  double hmaxrk = 0.1;
  double tolrk = 1e-13;
  int npasmxrk = 10000;
  
// 1. Calcular F y DF amb opham_fdf
  double f_o[4];      // F(X)
  double df_o[(4)*(3)]={0}; // DF(X)
  
  printf("=== Calcul opham ===\n");
  
  if (opham_fdf(m, hh, xx, cg, f_o, df_o, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, pendol_h, pendol, NULL) != 0) {
    fprintf(stderr, "Error en opham_fdf\n");
    return -1;
  }
  
  //Imprimir
  printf("F(X) = ");
  imprimir_matriu(1,dimF,f_o);
  printf("DF(X) = ");
  imprimir_matriu(dimF,dimX,df_o);
  
// 2. Calcular DF numéricamente con diferencias finitas
  printf("\n=== Calcul numeric (diferencies finites) ===\n");
  
  double delta = 1e-6;
  double df_num[(4)*(3)] = {0};
  // Perturbar +delta
  double xp[3];
  double xm[3];
  double f_p[4];
  double f_m[4];
  
  /*Derivació numèrica*/
  for (int j = 0; j < dimX; j++) {
    for (int k = 0; k < dimX; k++) {
      xp[k] = xx[k];//xp=(T,x0)
      xm[k] = xx[k];//xm=(T,x0)
    }
    
    xp[j] += delta;
    xm[j] -= delta;
    
    // Calcular F(X+delta) y F(X-delta)
    
    // F(X+delta) 
    if (opham_fdf(m, hh, xp, cg, f_p, NULL,h0rk, hminrk, hmaxrk, tolrk, npasmxrk, pendol_h, pendol, NULL) != 0) {
      printf("Error en opham_fdf para X+delta\n");
      return -1;
    }
    
    // F(X-delta) 
    if (opham_fdf(m, hh, xm, cg, f_m, NULL, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, pendol_h, pendol, NULL) != 0) {
      printf("Error en opham_fdf per X-delta\n");
      return -1;
    }
    
    // Diferencia finita centrada: ∂F_i/∂X_j ≈ [F_i(X+δ) - F_i(X-δ)] / (2δ)
    for (int i = 0; i < dimF; i++) {
      DF(i,j) = (f_p[i] - f_m[i]) / (2.0 * delta);
    }
  }
  
  imprimir_matriu(dimF,dimX,df_num);
  return 0;
}
#undef DF