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


int fminsearch(double (*funcpt)(double *,int),int N,double *xi,double *xf);

int fminunc(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,int method,double *xf);

#ifdef __cplusplus
}
#endif

#endif /* OPTIMIZE_H_ */
