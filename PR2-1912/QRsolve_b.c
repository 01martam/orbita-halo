#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Reescribim la qualsevol matriu en un vector que guarda cada columna una radera de l'altre*/
#define A(i, j) a[(i) + (j) * m] /* i = fila j= Columna m,p=nombre de files*/

/*Funcio que serveix per imprimir les "matrius" de la forma a[(i) + (j) * m] per consola.
m: nombre de files
n:nombre de columnes
*/
void imprimir_matriu(int m, int n,double *a) {
    for (int i = 0; i < m; i++) {
         printf("|");
        for (int j = 0; j < n; j++) {
            printf("%g\t", A(i, j));
        }
        printf("|\n");
    }
}

/*Funcio QRsolve que serveix per reslodre sistemes.
m:nombre de files.
n:nombre de columnes. m>=n
a->A(i,j): En entrar:A(i,j) és una matriu m*n del sistema lineal. i=0,...,m-1, j=0,..,n-1
           En sortir:Dins A(i,j) per j>=i hi ha la matriu R de la descomposició QR, dins A(i,k) i=k+1,...,m-1
           hi ha les coordenades de la k-essima reflexió de Householder que s'ha aplicat a k=0,..,n-1,
           és a dir els vectors u_k que hem calculat que són normats, els quals tenen la primera component 1 per construcció i per això no es guarda aquesta.
b:En entrar:Terme independent del sistema lineal sobredeterminat Ax=b.
  En sortir: Q^T*b
  Si b!=NULL no es fa servir.
tau:En sortir conte les tau_k a la posició k-essima per utilitzar l'algoritme.
x:En sortir conté la solució del sistema sobredeterminat . Si b!=NULL  no es fa servir.
tol: tolerancia per determinar la singularitat de A, usualment s'utilitza valors molts propers a 1.
passos: Si introdueixes qualsevol valor enter que no sigui 0, sortiran misatges on es mostrara pasos intermitjos del calcul.
        Si introdueixes 0, no imprimeix res.
La funció retorna
    0 si tot anat bé.
    -1 si la matriu A és singular.
    -2 Si la matriu A no és sobredeterminada.
*/
int QRsolve(int m, int n, double *a, double *b, double *tau, double *x, double tol,int pasos) {
    //asegurarnos de que la matriu es de la mida correcta.
    if(m<n){
        fprintf(stderr,"Error: ha de ser m>=n.\n");
        return -2;}

    int i, j, k;
    double sk, beta;
    // evitar fer un pas extra quan la matriu es determinada en lloc de sobredeterminada.
    int kmax = (m == n)? n - 1 : n;

    if(pasos!=0){
        printf("Descomposicio de A:\n");
    }

// Algoritme 1:Descomposició QR amb remplaçament
    for (k = 0; k < kmax; k++) {
        if(pasos!=0){
            printf("iterat k=%d:\n",k);
        }

        sk = 0;
        // Càlcul de ||a^(k)||_2 ^2
        for (i = k; i < m; i++) {
            sk = sk + A(i,k) * A(i,k);
        }

        // Detecció de matriu singular
        if (sk < tol) {
            return -1;
        }

        // Càlcul de s_k
        sk = sqrt(sk);//Treure el quadrat de la norma
        if (A(k,k) > 0) {
            sk = -sk;// -sign(a_1^(k))
        }

        //Utilitzem tau com a auxiliar
        tau[k] = A(k,k) - sk;

        // Normalizació del vector u^(k)
        for (i = k + 1; i < m; i++) {
            A(i,k) = A(i,k) / tau[k];
        }
        //Càlcul del multiplicador tau_k=-(a_{k,k}-s_k)/s_k
        tau[k] = -tau[k] / sk;

        //s_k*e_1
        A(k,k) = sk;

        //Arreglem les altres columnes
        for (j = k + 1; j < n; j++) {
            beta = A(k,j);

            for (i = k + 1; i < m; i++) {
                beta = beta + A(i,k) * A(i,j);
            }

            beta = beta * tau[k];
            A(k,j) = A(k,j) - beta;

            for (i = k + 1; i < m; i++) {
                A(i,j) = A(i,j) - beta * A(i,k);
            }
        }
        if(pasos!=0){
            printf("Matriu A despres del pas %d:\n", k);
            printf("A^(%d)=\n",k+1);
            imprimir_matriu(m, n, a);
        }


// Algoritme 2: Càlcul de Qb a partir de u^(k)
        if (b != NULL) {
            beta = b[k];

            for (i = k + 1; i < m; i++) {
                beta = beta + A(i,k) * b[i];
            }

            beta = beta * tau[k];
            b[k] = b[k] - beta;

            for (i = k + 1; i < m; i++) {
                b[i] = b[i] - beta * A(i,k);
            }
            if(pasos!=0){
                printf("vector b despres del pas %d:\n", k);
                printf("b^(%d)=\n",k+1);
                imprimir_matriu(m, 1, b);
            }
        }
    }
    if(pasos!=0){
        printf("tau=");
        imprimir_matriu(1, kmax, tau);
    }

// Resolució del sistema Rx = Qb
    if (b != NULL && x != NULL) {
        for (i = n - 1; i >= 0; i--) {
            x[i] = b[i];
            for (j = i + 1; j < n; j++) {
                x[i] = x[i]-A(i, j) * x[j];
            }
            x[i] = x[i]/A(i,i);
        }
    }
    return 0;
}


#undef A

