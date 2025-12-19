#include <math.h>

int pendol (int n, double t, double *x, double *f, double *prm) {
   if(n!=2&&n!=6){return-1;}
   f[0]=x[1];
   f[1]=-sin(x[0]);
   if (n==6) {
      f[2]=x[3];
      f[3]=-x[2]*cos(x[0]);
      f[4]=x[5];
      f[5]=-x[4]*cos(x[0]);
   }
   return 0;
}

int pendol_h(double x[], double *h, double *prm) {
    //h[0] 
    *h= 0.5 * x[1]*x[1] - cos(x[0]);   /* H(q,p) */
    
    ///* Gradiente del Hamiltoniano */
    //h[1] = sin(x[0]);   /* ∂H/∂q */
    //h[2] = x[1];        /* ∂H/∂p */
    return 0;
}

