#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "QRsolve_b.h"
#include "Newton_b.h"
#include "rk45_b.h"
#include "flux_b.h"

/*La funció opham_fdf() ha de tornar dins *f el valor F (X), i, si dff!=NULL, ha de
 tornar dins *dff la diferencial DF (X), guardada per columnes (tal com l’espera
 QRsolve()).
 Si dff==NULL, opham_fdf() no ha de fer servir dff.
 .

 on
 -m dimensio,
 -hh = h (nivell d’energia)
 -xx = X=(T,x) ∈ R1+n (perı́ode i condició incial),
 -*cg conté els coeficients (c0 , . . . , cn−1 , cn) de l’hiperplà g(x) = c0 x0 + c1 x1 + · · · + cn−1xn−1 + cn,
 -ham és un punter a funció que implementa el Hamiltonià H (l’argument *x de ham() conté les coordenades del punt x on es vol trobar el valor del Hamiltonià, que s’ha de tornar dins *h)
 -camp és un punter a una funció que implementa el camp f i les seves variacionals primeres, tal com ho fèiem a la Pràctica 4.
 -f sortira F(x)
 -df sortira DF(x)
 -parametres flux(): hminrk, hmaxrk, tolrk, npasmxrk,h0rk. (El paràmetre h0rk s’ha de prendre com a pas inicial per a la integració numèrica)

 La funcio retorna
    0 si tot anat bé
    -1 si falla en guarda memoria,
    -2 si no s'ha evaluat la funció ham() o camp()
    -3 si no s'ha evaluat la funció flux()


 */


/*
  x es Rⁿ
  f:Rⁿ->Rⁿ
  H:R^m×R^m->R
  
  n=2*m
  x=(q,p)
  f(x)=(dH/dp, -dH/dq )
  X=(T,x) es R^(n+1)
  F:R^(1+n)->R^(2+n)
  F(X)=[H(x)-h]
  [g(x)]
  [phiT(x)-x]
  g Rⁿ->R     
 
  xx R^(1+n)
  xx=(T,x) x Rⁿ
  x=(q,p)  q,p R^m
  f Rⁿ
  H R
  F R^(n+2)
  g R
 */

#define A(i,j) x_var[n+(j)*n+(i)]
#define DF(i,j) df[(i) + dimF*(j)]

int opham_fdf (
    int m, double hh, double *xx, double *cg, double *f, double *df,
    double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,
    int (*ham)(double *x, double *h, double *prm),
    int (*camp)(int n, double t, double *x, double *f, double *prm),
    double *prm)
{

  int i,j;
  int n=2*m;//dimensió del sistema hamiltonià
// 0. Condicions inicials
  //xx=(T,x) R^(n+1)
  double T=xx[0];
  double *x=xx+1;//x te n components

//1.Avaluar H(x)-h
  double Hx;// el que dona la funcio H(x)
  if(ham(x,&Hx,prm)!=0){
    fprintf(stderr,"Error en l'evaluacio de la funcio ham()");
    return -2;
  }
  f[0]=Hx-hh;//primera coordenada de F(x)

//2.Avaluar g(x)=c_n + sum_{i=0}^{n-1} c_i*x_i
  f[1]=cg[n];//terme independent (cg te n+1 elements)
  for (i=0;i<n;i++){
    f[1]+=cg[i]*x[i];// "segona component" de F(x)
  }

//3.Bifurcacio ¿df=NULL?
  double t1 = 0.0;
  double h1 = h0rk;
  if(df==NULL){//Nomes necesitem F
    // 4. Calcular phi_T(x) fent integracio
      double *x_int = (double *)malloc(n * sizeof(double));
      if (x_int == NULL) {
        fprintf(stderr, "Error en guardar memoria de x_int\n");
        return -1;
      }
      
      // condicions inicials
      for (i = 0; i < n; i++) {
        x_int[i] = x[i];
      }
      
      // Integrar sense variacionals
      if (flux(&t1, x_int, &h1, T, hminrk, hmaxrk, tolrk, npasmxrk, n, camp, prm,0) != 0) {
        fprintf(stderr, "Error en flux() per calcular phi_T(x)\n");
        free(x_int);
        return -3;
      }
      //x_int=phi_T(x)
      // 5. Calcular phi_T(x)-x i omplir en F
      for (i = 0; i < n; i++) {
        f[2 + i] = x_int[i] - x[i];//"tercera component" de F(x)
      }
      free(x_int);
      return 0;
  }
  else{
    //Necesitem F I DF
    // 4. Calcular phi_T(x) i Dphi_T(x) amb variacionals fent integracio
    // Dimensió per variacionals: n (sistema) + n*n (matriu A)
      int nva = n + n * n;//dimensió del sistema amb variacionals
      double *x_var = (double *)malloc(nva * sizeof(double));//x_var=(x+A)
      if (x_var == NULL) {
        fprintf(stderr, "Error en guardar memoria de x_var\n");
        return -1;
      }
      
      // Copiar x0 a les primeres n components
      for (i = 0; i < n; i++) {
        x_var[i] = x[i];
      }
      
      // Inicializar A(0) = I (identidad)
      for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
          A(i,j) = (i == j) ? 1.0 : 0.0;
        }
      }
      
      // Integrar el PVI (x, A)
      if (flux(&t1, x_var, &h1, T, hminrk, hmaxrk, tolrk, npasmxrk, nva, camp, prm,0) != 0) {
        fprintf(stderr, "Error en flux() per calcular variacionals\n");
        free(x_var);
        return -3;
      }
      
      // 5. Calcular phi_T(x)-x i omplir en F
      for (i = 0; i < n; i++) {
        f[2 + i] = x_var[i] - x[i];//"tercera component" de F(x)
      }
      
      // 6. Construir DF
      // DF es una matriu de (2+n) x (1+n)
      int dimF = 2 + n;
      //int dimX = 1 + n;
      
      //a) posar els 0's
      DF(0,0)=0.0;
      DF(1,0)=0.0;
      
      //// Inicializar amb zeros
      //for (i = 0; i < dimF * dimX; i++) {
       // df[i] = 0.0;
      //}
      
      
      //b)Calcular f(phi_T(x)) para d/dT Columna 0 I files 2 a 2+n-1
      double *f_phi = (double *)malloc(n * sizeof(double));
      if (f_phi == NULL) {
        fprintf(stderr, "Error en guardar memoria de f_phi\n");
        free(x_var);
        return -1;
      }
      
      if (camp(n, 0.0, x_var, f_phi, prm) != 0) {
        fprintf(stderr, "Error en l'evaluacio de la funcio camp() per phi_T\n");
        free(x_var);
        free(f_phi);
        return -2;
      }
      
      // d/dT(phi_T(x)) = f(phi_T(x))=f(x_int)=f_phi
      for (i = 0; i < n; i++) {
        DF(2 + i, 0) = f_phi[i];
      }
      
      //Fila 0 I 1 Columna 1 a n
      //f_phi=f(x)
      if (camp(n, 0.0, x, f_phi, prm) != 0) {
        fprintf(stderr, "Error en l'evaluacio de la funcio camp() per x\n");
        free(x_var);
        free(f_phi);
        return -2;
      }
      
      //c) DH(x) = (-f_m, ..., -f_{n-1}, f_0, ..., f_{m-1})
      // Fila 0, columnes 1..n (derivadas respecto a x)
      for (i = 0; i < m; i++) {
        DF(0, 1 + i) = -f_phi[m + i];      // respecto a q_i: -f_{m+i}
        DF(0, 1 + i + m) = f_phi[i];     // respecto a p_i: f_i
      }
      
      //Fila 1, columnes 1..n
      //d) Dg(x) = (c_0, c_1, ..., c_{n-1})
      for (i = 0; i < n; i++) {
        DF(1, 1 + i) = cg[i];
      }
      
      // Dphi_T(x) - I
      //e) A(i,j) contiene Dphi_T(x)[i][j]
      for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
          DF(2 + i, 1 + j) = A(i,j) - ((i == j) ? 1.0 : 0.0);
        }
      }
      
      // Liberar memoria
      free(x_var);
      free(f_phi);
      
      return 0;
  }
  
}

#undef A
#undef DF






/*La funció opham ha de trobar, a partir de X (0) ∈ R1+n, el perı́ode i una condició
 incial dins l’hiperplà definit pels coeficients c0 , . . . , cn de l’òrbita periòdica d’energia
 h del sistema Hamiltonià definit per ham() i camp(). Per a això, ha de fer servir el
 programa newton().

 En cas d’èxit, ha de tornar X (k) dins *xx i el nombre d’iterats
 emprats via return. En cas de no convergència, ha de tornar el darrer iterat dins xx i −1 via return.
 on
 -m = m,
 -hh = h (nivell d’energia),
 -xx = X (0) ∈ R1+n (perı́ode i condició incial aproximats),
 -*cg conté els coeficients de l’hiperplà (c0 , . . . , cn−1, cn),
 -ham és un punter a funció que implementa el Hamiltonià H (l’argument *x de ham() conté les coordenades del punt x on es vol avaluar el Hamiltonià), i camp és un punter a
     una funció que implementa el camp f i les seves variacionals primeres, de manera
     anàloga a la primera pràctica.
 -Els arguments tol i maxit són els mateixos que a newton().
 
 La funcio retorna
 (el nombre de iterats si)>0 tot va bé
 -1 si no hi ha convergencia
 -2 Si hi ha hagut algun problema en resoldre algun QRsolve
 -3 si no s'ha pogut avaluar la funcio fdf
 -4 si no s'ha pogut asignar memòria

 A cada iterat de Newton, escriu el nombre d’iterat, ∥F(X(k))∥2 i ∥∆X(k)∥2 .
 Per exemple:
 opham(): it 0 nf 0.00040621 nc 1.08646E-06
 opham(): it 1 nf 1.67035E-07 nc 1.08125E-06
 opham(): it 2 nf 1.04319E-09 nc 2.36645E-11
 opham(): it 3 nf 1.33993E-12
 */


//Versio 1 programar el Newton dintre de opham
int opham2 (int m, double hh, double *xx, double tol, int maxit, double *cg,
           double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,
           int (*ham)(double *x, double *h, double *prm),
           int (*camp)(int n, double t, double *x, double *f, double *prm),
           double *prm)
{
  //Fixemnos que el que hi ha dintre de fdf no l'estem donan nombre si no que diem quantes variables crida i quin tipus, perfectament podriem treure els noms.
  int iterats,i;
  int n=2*m,dimF=(2+n),dimX=(1+n);
  double norm_f,norm_y,*y,*f,*df,*tau;
  //Reserva memoria
  f=(double *)malloc(dimF*sizeof(double));
  df=(double *)malloc(dimF*dimX*sizeof(double));
  tau=(double *)malloc(dimX*sizeof(double));
  y=(double *)malloc(dimX*sizeof(double));
  
  if (y == NULL || f == NULL || df == NULL || tau == NULL){
    fprintf(stderr,"Error en asignar memoria.\n");
    free(y); free(f); free(df); free(tau);
    return -4;
  }
  
  for (iterats=0;iterats<maxit;iterats++){
    //Avaluar F(x_k) i DF(x_k)
    if(opham_fdf (m, hh, xx, cg, f, df,h0rk, hminrk, hmaxrk, tolrk, npasmxrk,ham,camp,prm)!=0){
      fprintf(stderr,"Error en la avaluació de la funció en l'iterat%d\n",iterats);
      free(y); free(f); free(df); free(tau);
      return -3;
    }
    //Calcular la norma de F(x)
    norm_f=norma_2(dimF,f);
    
    //Verifiquem la convergencia
    if(norm_f<tol){
      printf("opham():it %d norm_f %e \n",iterats,norm_f);
      iterats++;
      printf("Newton convergeix despres de %d iteracions\n",iterats);
      free(y); free(f); free(df); free(tau);
      return iterats;//+1 si contem 0,1,2,3,4,...
    }
    
    //Resoldre el sistema DF(x_k)y_k=F(x_k)
    if (QRsolve(dimF, dimX, df, f, tau, y, tol,0) != 0) {
      fprintf(stderr, "Error resolent el sistema en l'iterat %d\n", iterats);
      free(y); free(f); free(df); free(tau);
      return -2;
    }
    //Calcular la norma de y
    norm_y=norma_2(dimX,y);
    printf("opham():it %d norm_f %e norm_y %e\n",iterats,norm_f,norm_y);
    
    
    //Actualitza el valor de x
    for(i=0;i<dimX;i++){
      xx[i]=xx[i]-y[i];
    }
  }
  //Newton no convergeix
  fprintf(stderr,"Error : Newton no ha convergit amb %d iteracions.\n",maxit);
  free(y); free(f); free(df); free(tau);
  return -1;
}
  
//Versio 2 Aprofitant la programació de newton
typedef struct {
    int m;
    double hh;
    double *cg;
    double h0rk;
    double hminrk;
    double hmaxrk;
    double tolrk;
    int npasmxrk;
    int (*ham)(double *x, double *h, double *prm);
    int (*camp)(int n, double t, double *x, double *f, double *prm);
    double *prm;
  } opham_params;

  
  
int fdf(int m, int n, double *x, double *f, double *df, void *params)
  {
    // Hacer un casting del puntero void a la estructura esperada (opham_params_t)
    opham_params *p = (opham_params *)params;
    
    // Usar los parámetros pasados en la estructura
    return opham_fdf(p->m, p->hh, x, p->cg, f, df, p->h0rk, p->hminrk, p->hmaxrk, 
                     p->tolrk, p->npasmxrk, p->ham, p->camp, p->prm);
}

int opham (int m, double hh, double *xx, double tol, int maxit, double *cg,
             double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,
             int (*ham)(double *x, double *h, double *prm),
             int (*camp)(int n, double t, double *x, double *f, double *prm),
             double *prm)
{
  int n=2*m,dimF=(2+n),dimX=(1+n);

  opham_params params = {m,hh, cg, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, ham, camp, prm};
  

  return newton(dimF, dimX, xx, fdf ,maxit, tol,1,&params);

}

  
  


