/*
 * optimc.h
 *
 *  Created on: Mar 16, 2014
 *      Author: HOME
 */

#ifndef OPTIMC_H_
#define OPTIMC_H_

#include "neldermead.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct opt_set* opt_object;

opt_object opt_init(int N);

struct opt_set{
	int N;
	double objfunc;
	double eps;
	double gtol;
	double stol;
	double ftol;
	double xtol;
	int MaxIter;
	int Iter;
	int Method;
	int retval;
	char MethodName[50];
	double xopt[1];
};

void summary(opt_object obj);

void setMaxIter(opt_object obj,int MaxIter);

void setTOL(opt_object obj,double gtol,double stol,double ftol,double xtol);

int fminsearch(double (*funcpt)(double *,int),int N,double *xi,double *xf);

double fminbnd(double (*funcpt)(double),double a, double b);

int fminunc(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,int method,double *xf);

int fminqn(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,
		double delta,int method,double *dx,double fsval,double *xf);

double brent_local_min(double (*funcpt)(double ),double a, double b, double t, double eps, double *x);

void optimize(opt_object obj,double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,
		int method);

void free_opt(opt_object object);

#ifdef __cplusplus
}
#endif

#endif /* OPTIMC_H_ */
