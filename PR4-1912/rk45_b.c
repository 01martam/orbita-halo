#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*Implementació del metode runge-kutta de 6 etapes d'ordres 4 i 5 d'un pas adaptatiu
de manera que permeti resoldre qualsevol sistema d'equacion de dimensio n
t:En entrar temps t_{k}
  En sortir temps t_{k+1}
x:En entrar,x[i=0,...,n-1] han de ser les entrades x_{k}
  En sortir  les entrades de x_{k+1}
h:En entrar, pas h_{k} pas proposat
  En sortir pas h_{k+1} h nou
hmax: El maxim valor que pot pendre h (valor absolut )
hmin:  El minim valor que pot pendre h (valor absolut )
tol: tolerancia absoluta epsilon (epsilon*)
n:dimensio del camp n
camp: funció que calcula f(t,x) i guarda al vector f.
prm: llista de parametres del camp f (si escau)
passos:text que mostra els passos intermitjos fins obtenir l'h.
    1  mostra  misatge
    qualsevol altre no hi ha misatges

en retorna:
    0 si el pas s'acepta,
    1 si no s'acepta el pas,
    -1 si falla memoria en guarda memoria,
    -2 no s'evalua ve el camp
    -3 Un dels calculs es NaN(Not a number)

  */
//mingw32-make

#define A(i,j) (((i) <= (j)) ? 0.0 : (a[(j) + (i)*((i)-1)/2]))//fer que si esta per sobre de la diagonal son 0 si no son els valors de a
#define K(i,j) k[(i) + n*(j)]// i index de la variable, j els diferents kappas k_j
#define MIN(x,y) ((x)<(y)? (x):(y))
#define MAX(a,b) ((a)>(b)? (a):(b))
#define SIGNE(h)  ((h) > 0 ? 1 : -1)

/*Coeficients RK45 (estatics)*/
	static double c[6] = {0., 1./4., 3./8., 12./13., 1., 1./2.};
	static double a[15] = {1./4., 3./32., 9./32., 1932./2197.,
		-7200./2197.,7296./2197., 439./216., -8., 3680./513.,
		-845./4104., -8./27., 2., -3544./2565., 1859./4104., -11./40.};//matriu triangular superior sense guardar els 0 guardada per files
	static double b[6] = {25./216., 0., 1408./2565., 2197./4104., -1./5.,0.};//ordre 4
	static double bp[6] = {16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.};//ordre 5
	static double E[6] = {-1./360., 0., 128./4275., 2197./75240., -1./50., -2./55.};//diferencia de l'ordre 5 amb el 4



//tol=epsilon*
//err = epsilon o e_k
int rk45(double *t, double *x,
         double *h, double hmin, double hmax,
         double tol, int n,
         int (*camp)(int n, double t, double *x, double *f, double *prm), double *prm,int passos){

    //assegurem que hmin i hmax siguin positius.
    hmin=fabs(hmin);
    hmax=fabs(hmax);

    //Asegurem que l'h està dintre del interval, però conservant el signe per conservar si l'integració es foward o backward [hmin,hmax]
    (*h)=SIGNE(*h)*MIN(hmax,MAX(hmin,fabs(*h)));
    // Reservar dinámicament les kappa's
    /*Els indexos de les kappas en lloc d'anar de 1 a 6 va de 0 a 5*/
    double *k = malloc(6* n * sizeof(double));
    if (k == NULL) {
    fprintf(stderr, "Error al reservar memoria.\n");
    return -1;
    }
    /*si pensem com si fos una matriu llavors les diferents kappas estan guardades per columnes es a dir |k1[0,..,n-1],k2[0,..,n-1],...,k6[0,..,n-1]|
    | |  | ...  ||
    |k1 k2 ... k6|
    | |  | ...  ||
    //K(componente, kappa)
    */

    //Comencem calculan la primera kappa on no canvia t i x
    /*κ_1*/
    if (camp(n, *t, x, &K(0,0), prm)!=0) {free(k);return -2;}

    //Reservar memoria per una Variable auxiliar
    // calcular les x que utilitzarem per calcular kappa
    double *extraXk = (double *)malloc(n * sizeof(double));
    if (extraXk==NULL) {
            free(k);
            fprintf(stderr,"Error al reservar memoria.\n");
            return -1;
    }
    while(1){
        // Bucle per calcular les etapes κ_2..κ_6
        for (int etapa=1; etapa<6; etapa++) {
            // Construir extraXk = x + h * sum_j A(etapa,j) * kappa_j
            for (int i=0; i<n; i++) {// i=component del vector
                extraXk[i] = x[i]; // extraXk = x
                for (int j=0; j<etapa; j++) {//j=
                  //K(component, kappa)
                    extraXk[i] += (*h) * A(etapa,j) * K(i,j);//Construir extraXk = x + h * sum_j A(etapa,j) * kappa_j
                }
            }


            // Avaluar el camp en aquesta etapa
            if (camp(n, *t + c[etapa]*(*h), extraXk, &(K(0,etapa)), prm) != 0) {free(k);free(extraXk);return -2;}//k_{etapa}
        }

        //Anem a calcular la nostra x5
        double norm=0.0, tolrel = 0.0, err = 0.0;
        /*
        tol error màxim acceptable que volem
        tolrel error tolerable relatiu adaptada al pas
        err Error estimat de l'aproximació obtinguda
        norm indicador de la magnitud dels increments de Runge-Kutta
                   no controla el pas adaptatiu, nomes serveix com una mena de alarma per detectar que la solució s'ha tornat numéricament inválida.
        */
        //Reutilitzarem extraXk per guardar la nostra prediccio
        for (int i=0; i<n; i++) {
            //Calcul de x5
            extraXk[i] = 0.0;
            double diff = 0.0;
            for (int j=0; j<6; j++) {
              //K(componente, kappa)
            extraXk[i] += bp[j] * K(i,j);// sum_j bp[j] * kappa_j
            diff += E[j] * K(i,j);//sum_j (bp[j]-b[j]) * kappa_j
            }
            norm += fabs(extraXk[i]);//||sum_j bp[j] * kappa_j||_1 (suma del valor absolut dels component d'un vector)
            err+=fabs(diff);// ||sum_j (bp[j]-b[j]) * kappa_j||_1
            extraXk[i] = x[i] + (*h) * extraXk[i];//x5 proposat
            tolrel += fabs(extraXk[i]);//||x5||_1
        }

        if (isnan(norm) || isnan(err)||isnan(tolrel)) {free(k);free(extraXk);return -3;}//Comproba si es NaN(Not a number), es a dir valors que poden apareixer per  divisions per 0, arrels de nombres negatius,etc.
        err = fabs(*h) * (err);
         tolrel = tol * (1.0 + tolrel/100.0);//||x7-x8||_1
         //tolrel = tol;//no canvia gairé si agafo un o el altre
        //Veiem si aquest pas es valid i el modifiquem el pas
//Error d'aproximació aceptable
         if (err < tolrel) {
            //Pass bó
            *t += (*h);
             for (int i=0; i < n; i++){x[i] = extraXk[i];}//X5 nou
            //Error molt petit
             if (err <= ldexp(tolrel, -5)) {
                //e_k+1<<2^5*epsilon*
                 *h = SIGNE(*h) * MIN(hmax, fabs(*h) * 1.8);
                if (passos==1){printf("(1)h=%lf\n",*h);}
                free(k); free(extraXk);
                return 0;
                //multipliquem per 2*rho si no es pasa de hmax
            }else{
            //Error no massa petit
             *h = SIGNE(*h)*MIN(hmax, MAX(hmin, fabs(*h)*0.9 * pow(tolrel / err,0.2)));
             free(k);free(extraXk);
              if (passos==1){printf("(2)h=%lf\n",*h);}
             return 0;}
         }
//Error d'aproximacio no acceptable
        //pass masa petit
         if (fabs(*h) <= hmin) {
             fprintf(stderr, "h més petita que la minima.\n");
            *t += (*h);
            *h = SIGNE(*h)*hmin;
            for (int i=0; i < n; i++){x[i] = extraXk[i];}//X5 nou
             if (passos==1){printf("(3)h=%lf\n",*h);}
             free(k);free(extraXk);
             return 0;
         }
         *h = SIGNE(*h)*MAX(hmin, fabs(*h)*0.9 * pow(tolrel / err,0.2));
         if (passos==1){printf("(4)h=%lf\n",*h);}
    }
    if (passos==1){printf("(5)h=%lf\n",*h);}
    free(k);free(extraXk);
    return 0;
}




#undef A
#undef K
#undef MAX
#undef MIN
