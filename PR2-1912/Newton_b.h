double norma_2(int m, double *f);
//int newton(int m, int n, double *x, int (*fdf)(int m, int n,double *x,double *f,double *df),int maxit, double tol,int passos);
int newton(int m, int n, double *x, int (*fdf)(int m, int n,double *x,double *f,double *df, void *prm),int maxit, double tol,int passos,void *prm);
