#ifndef LNSRCHMP_H_
#define LNSRCHMP_H_

#include "matrix.h"

#define EPSILON 2.7182818284590452353602874713526624977572



#ifdef __cplusplus
extern "C" {
#endif

double macheps();

double pmax(double a, double b);

double pmin(double a, double b);

double signx(double x);

double l2norm(double *vec, int N);

int stopcheck_mt(double fx, int N, double *xc, double *xf, double *jac, double *dx, double fsval, double gtol, double stol, int retval);

int stopcheck2_mt(double fx, int N, double fo, double *jac, double *dx, double eps,double stoptol, double functol, int retval);

void grad_fd(double(*funcpt)(double *, int), double *x, int N, double *dx, double eps, double *f);

void grad_cd(double(*funcpt)(double *, int), double *x, int N, double *dx, double eps, double *f);

int cstep(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy, double *stp, double *fp, double *dp, int *brackt,
	double  stpmin, double stpmax);

int cvsrch(double(*funcpt)(double *, int), double *x, double *f, double *g, double *stp, double *s, int N, double *dx, double maxstep,
	int MAXITER,double eps2,double ftol, double gtol, double xtol);

int lnsrchmt(double(*funcpt)(double *, int), double *xi,double *f, double *jac,double *alpha, double *p, int N, double *dx, double maxstep, int MAXITER,
		double eps2,double ftol, double gtol, double xtol, double *x);

#ifdef __cplusplus
}
#endif

#endif /* LNSRCHMP_H_ */
