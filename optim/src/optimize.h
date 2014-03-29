/*
 * optimize.h
 *
 *  Created on: Mar 16, 2014
 *      Author: HOME
 */

#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

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

int fminsearch(double (*funcpt)(double *,int),int N,double *xi,double *xf);

double fminbnd(double (*funcpt)(double),double a, double b);

int fminunc(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,int method,double *xf);

int fminqn(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,
		double delta,int method,double *dx,double fsval,double *xf);

void optimize(opt_object obj,double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,
		int method);

#ifdef __cplusplus
}
#endif

#endif /* OPTIMIZE_H_ */
