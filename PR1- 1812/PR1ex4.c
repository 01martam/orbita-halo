#include <stdio.h> //print f
#include <stdlib.h>// reserva de memoria
#include <math.h>
#include <time.h>

#include "QRsolve_b.h" // On esta el QRsolve i imprimir_matriu i matriu_x_vector i matriu_x_matriu

#define A(i, j) a[(i) + (j) * m] /* i = fila j= Columna */


int main(int argc, char *argv[]){
    int n;
 /*
 * Línia de comandes:
    int n
 */
    //Mirar si s'ha introduit un valor i si esta mirar que es un enter positiu.
    if (argc!=2|| sscanf(argv[1],"%d",&n)!=1|| n <= 0) {
      fprintf(stderr,"%s n (enter positiu)\n",argv[0]);
      return -1;
   }

    srand(time(0));
    int m=n,i,j;
    //printf("Introdueix un valor (n):");
    //scanf("%d", &n);
    int kmax = (m == n)? n - 1 : n;

    double *a,*b, *x,*tau;
    double tol= 1e-12,error,max_error=0;

    //asignem memoria
    tau=(double*)malloc(kmax*sizeof(double));
    if (tau==NULL){
        fprintf(stderr,"Error en asignar memoria per tau.\n");
        return -1;
    }

    a=(double*)malloc(m*n*sizeof(double));
    if (a==NULL){
        fprintf(stderr,"Error en asignar memoria per la matriu A.\n");
        return -1;
    }

    b=(double*)malloc(m*sizeof(double));
    if (b==NULL){
        fprintf(stderr,"Error en asignar memoria per el vector b.\n");
        return -1;
    }
    x=(double*)malloc(n*sizeof(double));
    if (x==NULL){
        fprintf(stderr,"Error en asignar memoria per x.\n");
        return -1;
    }
    //Generar matriu aleatoria i guardar copia
    for (j=0;j<n;j++){
        for(i=0;i<m;i++){
            A(i,j)=((double)rand()/RAND_MAX);
        }
    }
    printf("Tenim la matriu\n");
    imprimir_matriu(m,n,a);


    //Generar b tal que la solució sigui (1,1,...,1)
    for(i=0;i<m;i++){
        b[i]=0;
        for (j=0;j<n;j++){
            b[i]=b[i]+A(i,j)*1.0;
        }
    }

    printf("I el vector b\n");
    imprimir_matriu(m, 1,b);
    printf("\n\n");

    //Resoldre el sistema creat
    QRsolve(m,n,a,b,tau,x,tol,0);

    printf("\nLa solucio del sistema es :\n");
    printf("x=");
    imprimir_matriu(1, n, x);

    //calcular error maxim
    for(i=0;i<n;i++){
        error=fabs(x[i]-1.0);
        if(error>max_error){
            max_error=error;
        }
    }
    printf("error maxim:%e\n",max_error);

    free(tau);free(x);free(a);free(b);
    return 0;
}
#undef A
