#include <stdio.h> //print f
#include <stdlib.h>// reserva de memoria
#include <math.h>
#include "QRsolve_b.h"// On esta el QRsolve i imprimir_matriu i matriu_x_vector i matriu_x_matriu

//#define A(i, j) a[(i) + (j) * m] /* i = fila j= Columna */


int main(){
    int m = 3, n = 3;
    int kmax = (m == n)? n - 1 : n;
    double b[]={1,2,3},x[n],tau[kmax], tol= 1e-12;
    double a[] = {
        12, 6, -4,   // columna 0
        -51, 167, 24,   // columna 1
        4, -68, -41,
    };

    printf("Tenim la matriu\n");
    imprimir_matriu(m,n,a);

    printf("I el vector b\n");
    imprimir_matriu(m, 1,b);
    printf("\n\n");

    QRsolve(m,n,a,b,tau,x,tol,0);

    printf("\nLa solucio del sistema es :\n");
    printf("x=");
    imprimir_matriu(1, n, x);

    printf("\nR=");
    imprimir_matriu(m,n,a);

    printf("\nQb=");
    imprimir_matriu(m, 1,b);


    return 0;
}
//#undef A
