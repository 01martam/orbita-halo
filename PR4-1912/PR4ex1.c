#include <stdio.h>
#include <math.h>
#include "flux_b.h"
#include "rk45_b.h"

int camp(int n, double t, double *x, double *f, double *prm){
  double a=prm[0];
  if (n!=2&&n!=6){
    fprintf(stderr,"ERROR en la dimensió (n=2 per flux, n=6 per la variacionals).\n");
    return -1;
  }
  else if(n==2){
    f[0]=a*x[0]-x[1];
    f[1]=x[0]+a*x[1];
    return 0;
  }
  else {
    n=2;
    #define A(i,j) x[n+(j)*n+(i)]
    f[0]=a*x[0]-x[1];
    f[1]=x[0]+a*x[1];
    // Calcula Jacobiana
    double J[2][2];
    J[0][0]=a;//df0/dx
    J[1][0]=1;//df1/dx
    
    J[0][1]=-1.0;//df0/dy
    J[1][1]=a;//df1/dy
    
    //variacional
    f[2]=J[0][0]*A(0,0)+J[0][1]*A(1,0);
    f[3]=J[1][0]*A(0,0)+J[1][1]*A(1,0);
    f[4]=J[0][0]*A(0,1)+J[0][1]*A(1,1);
    f[5]=J[1][0]*A(0,1)+J[1][1]*A(1,1);
    #undef A
    return 0;
    
    
  }
  
}

int funcio_exacta(double t,double *f, double *x0,double *prm){
  double a=prm[0];
  double expo=exp(a*t), cost=cos(t), sint=sin(t);//evitar calcular la mateixa cosa diversos cops.
  f[0]=expo*(x0[0]*cost-x0[1]*sint);
  f[1]=expo*(x0[0]*sint+x0[1]*cost);
  return 0;
}


int main(int argc, char *argv[]) {
/*linea de comandes */
  double t=0,T,h,hmax,hmin,tol=1e-12,prm[1],x_inicial[2];//parametre a
  int n,npasmax;
  if(argc<10||sscanf(argv[1],"%lf",&prm[0])!=1
       ||sscanf(argv[2],"%lf",&T)!=1
       ||sscanf(argv[3],"%lf",&x_inicial[0])!=1
       ||sscanf(argv[4],"%lf",&x_inicial[1])!=1
       ||sscanf(argv[5],"%lf",&h)!=1
       ||sscanf(argv[6],"%lf",&hmin)!=1
       ||sscanf(argv[7],"%lf",&hmax)!=1
       ||(sscanf(argv[8],"%d",&n)!=1||(n!=2 && n!=6))
       ||sscanf(argv[9],"%d",&npasmax)!=1)
  {
    fprintf(stderr,"%s a T x0 y0 h hmin hmax n npassmax\n",argv[0]);
    return -1;
  }
  
  /*Fi de linia de comandes*/
    double x[n];
    x[0]=x_inicial[0];
    x[1]=x_inicial[1];
    #define A(i,j) x[2+(j)*2+(i)]
    //Inicialitzar La identitat (condicio inicial de variacionals)
    if (n==6){//n(n+1)
        for (int i=0;i<2;i++){
            for(int j=0;j<2;j++){
                A(i,j)=(i==j)?1.0:0.0;
            }
        }
    }

    int resultat =flux(&t,x,&h,T,hmin,hmax,tol,npasmax,n,camp,prm,1);
    if (resultat== 0) {
        printf("t_final = %f\n", t);
        printf("x_final = (%.15lf, %.15lf)\n", x[0], x[1]);
        double f[2];
        funcio_exacta(t,f, x_inicial, prm);
        printf("x_exacta = (%.15lf, %.15lf)\n", f[0], f[1]);
        if (n==6){//n(n+1)
            printf("A=");
            for (int i=0;i<2;i++){
                for(int j=0;j<2;j++){
                    printf("%.6lf\t",A(i,j));
                }
                printf("\n");
            }
        }
    } else {
        printf("Error en l'integració: %d\n", resultat);
    }
    #undef A
    return 0;
}
