/*
 * neldormead.h
 *
 *  Created on: Jan 5, 2014
 *      Author: HOME
 */

#ifndef NELDORMEAD_H_
#define NELDORMEAD_H_

#include "brent.h"

#ifdef __cplusplus
extern "C" {
#endif

int nel_min(double (*funcpt)(double *,int),double *xc,int N,double *dx,double fsval,int MAXITER,int *niter,
		double eps,double *xf);

#ifdef __cplusplus
}
#endif


#endif /* NELDORMEAD_H_ */
