#include <stdio.h> //print f
#include <stdlib.h>// reserva de memoria
#include <math.h>
#include "QRsolve_b.h" // On esta el QRsolve i imprimir_matriu i matriu_x_vector i matriu_x_matriu

//#define A(i, j) a[(i) + (j) * m] /* i = fila j= Columna */

int main(){
    int m = 3, n = 2;
    int kmax = (m == n)? n - 1 : n;
    double b[]={1,3,2},x[n],tau[kmax],tol= 1e-12;
    double a[] = {
        0, 0, 5,   // columna 0
        -4, 0, -2,   // columna 1
    };


    printf("Tenim la matriu\n");
    imprimir_matriu(m,n,a);

    printf("I el vector b\n");
    imprimir_matriu(m, 1,b);
    printf("\n\n");

    QRsolve(m,n,a,b,tau,x,tol,1);


    printf("\nLa solucio del sistema es :\n");
    printf("x=");
    imprimir_matriu(1, n, x);

    return 0;
}
//#undef A
