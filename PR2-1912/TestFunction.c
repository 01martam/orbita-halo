#include <stdio.h>
#include <math.h>

int function1(int m, int n, double *x, double *f, double *df, void *prm)
{
	if (m!=3||n!=3){
		printf("function1(): Error en les dimensions de la funció F.\n");
		return -1;
	}
	#define X1 x[0]
	#define X2 x[1]
	#define X3 x[2]
	#define F1 f[0]
	#define F2 f[1]
	#define F3 f[2]
	#define DF(i,j) df[(i)+(m)*(j)] /* i=fila, j=columna */
	F1 = X1+X2*X3;
	F2 = X2-X1*X1;
	F3 = cos(X3);
	DF(0,0) = 1;//d/dx1 f1
	DF(0,1) = X3;//d/dx2 f1
	DF(0,2) = X2;//d/dx3 f1
	
	DF(1,0) = -2*X1;//d/dx1 f2
	DF(1,1) = 1;//d/dx2 f2
	DF(1,2) = 0;//d/dx3 f2
	
	DF(2,0) = 0;//d/dx1 f3
	DF(2,1) = 0;//d/dx2 f3
	DF(2,2) = -sin(X3);//d/dx3 f3
	#undef X1
	#undef X2
	#undef X3
	#undef F1
	#undef F2
	#undef F3
	#undef DF
	return 0;
}

int funcio_fractal(int m, int n, double *x, double *f, double *df, void *prm){
    /*
    Escribim la funcio complexa com una funcio a R^2
    x[0]=Re(z=:x_0), x[1]=Im(z)=:x_1
    p(z)=z^3-1=(x_0+ix_1)^3-1
    parts real de p:(x_0)^3-3*(x_0)*(x_1)^2-1
    part imaginaria de p:3*(x_0)^2*(x_1)-(x_3)^3
    */
    if (m!=2||n!=2){
		printf("function1(): Error en les dimensions de la funció F.\n");
		return -1;
	}
  #define X x[0]
	#define Y x[1]
	#define DF(i,j) df[(i)+(m)*(j)]
  double x2=X*X, y2=Y*Y;
	f[0]=X*x2-3*X*y2-1.0;
	f[1]=3*x2*Y-Y*y2;
  
  double aux1= 3*x2-3*y2, aux2=6*X*Y;
	DF(0,0)=aux1;//d/dx f0
	DF(0,1)=-aux2;//d/dy f0
	
	DF(1,0)=aux2;//d/dx f1
	DF(1,1)=aux1;//d/dy f1
	#undef X
	#undef Y
	#undef DF
	return 0;

}
