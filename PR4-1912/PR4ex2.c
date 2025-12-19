#include <stdio.h>
#include <math.h>
#include "flux_b.h"
#include "rk45_b.h"

int EDO_lorenz(int n, double t, double *x, double *f, double *prm){
  double sigma=prm[0];
  double rho=prm[1];
  double beta=prm[2];
  if (n!=3&&n!=12){
    fprintf(stderr,"ERROR en la dimensió (n=3 per flux, n=9 per la variacionals).\n");
    return -1;
  }
  else if (n==3){
    //x=x[0] y=x[1] z=x[2]
    f[0]=sigma*(x[1]-x[0]);
    f[1]=rho*x[0]-x[1]-x[0]*x[2];
    f[2]=-beta*x[2]+x[0]*x[1];
    return 0;
  }
  else{
    n=3;
    //x=x[0] y=x[1] z=x[2]
    #define A(i,j) x[n+(j)*n+(i)]
    f[0]=sigma*(x[1]-x[0]);
    f[1]=rho*x[0]-x[1]-x[0]*x[2];
    f[2]=-beta*x[2]+x[0]*x[1];
    // Calcula Jacobiana
    double J[3][3];
    J[0][0]=-sigma;//df0/dx
    J[1][0]=rho-x[2];//df1/dx
    J[2][0]=x[1];//df2/dx

    J[0][1]=sigma;//df0/dy
    J[1][1]=-1.0;//df1/dy
    J[2][1]=x[0];//df2/dy

    J[0][2]=0;//df0/dz
    J[1][2]=-x[0];//df1/dz
    J[2][2]=-beta;//df2/dz
    //Variacional
    f[3]=J[0][0]*A(0,0)+J[0][1]*A(1,0)+J[0][2]*A(2,0);
    f[4]=J[1][0]*A(0,0)+J[1][1]*A(1,0)+J[1][2]*A(2,0);
    f[5]=J[2][0]*A(0,0)+J[2][1]*A(1,0)+J[2][2]*A(2,0);

    f[6]=J[0][0]*A(0,1)+J[0][1]*A(1,1)+J[0][2]*A(2,1);
    f[7]=J[1][0]*A(0,1)+J[1][1]*A(1,1)+J[1][2]*A(2,1);
    f[8]=J[2][0]*A(0,1)+J[2][1]*A(1,1)+J[2][2]*A(2,1);

    f[9]=J[0][0]*A(0,2)+J[0][1]*A(1,2)+J[0][2]*A(2,2);
    f[10]=J[1][0]*A(0,2)+J[1][1]*A(1,2)+J[1][2]*A(2,2);
    f[11]=J[2][0]*A(0,2)+J[2][1]*A(1,2)+J[2][2]*A(2,2);

    #undef A
    return 0;
  }
}


int main(int argc, char *argv[]) {
/*linea de comandes */
  double t,T,h,hmax,hmin,tol,prm[3]={10.0,28.0,8./3.},x_inicial[3];//parametre sigma, rho, beta
  int n,npasmax;


  if(argc<12
       ||sscanf(argv[1],"%lf",&T)!=1
       ||sscanf(argv[2],"%lf",&t)!=1//t0
       ||sscanf(argv[3],"%lf",&x_inicial[0])!=1//x0
       ||sscanf(argv[4],"%lf",&x_inicial[1])!=1//y0
       ||sscanf(argv[5],"%lf",&x_inicial[2])!=1//z0
       ||sscanf(argv[6],"%lf",&h)!=1
       ||sscanf(argv[7],"%lf",&hmin)!=1
       ||sscanf(argv[8],"%lf",&hmax)!=1
       ||(sscanf(argv[9],"%d",&n)!=1||(n!=3 && n!=12))
       ||sscanf(argv[10],"%lf",&tol)!=1
       ||sscanf(argv[11],"%d",&npasmax)!=1){
       fprintf(stderr,"%s T t0 x0 y0 z0 h hmin hmax n tol npasmax\n",argv[0]);
    return -1;
  }

  /*Fi de linia de comandes*/

    double x[n];
    x[0]=x_inicial[0];
    x[1]=x_inicial[1];
    x[2]=x_inicial[2];
    #define A(i,j) x[3+(j)*3+(i)]
    //Inicialitzar La identitat (condicio inicial de variacionals)
    if (n==12){//n(n+1)
        for (int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                A(i,j)=(i==j)?1.0:0.0;
            }
        }
    }

    int resultat =flux(&t,x,&h,T,hmin,hmax,tol,npasmax,n,EDO_lorenz,prm,1);
    if (resultat== 0) {
        printf("t_final = %f\n", t);
        printf("x_final = (%lf, %lf, %lf)\n", x[0], x[1], x[2]);
        if (n==12){//n(n+1)
            printf("A=\n");
            for (int i=0;i<3;i++){
                printf("[");
                for(int j=0;j<3;j++){
                    printf("%.15lf\t",A(i,j));
                }
                printf("]\n");
            }
        }
    } else {
        printf("Error en l'integració: %d\n", resultat);
    }
 
    #undef A
    return 0;
}
