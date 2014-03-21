/*
 * optimize.c
 *
 *  Created on: Mar 16, 2014
 *      Author: Rafat Hussain
 */

#include "optimize.h"

int fminsearch(double (*funcpt)(double *,int),int N,double *xi,double *xf) {
	int i,retval,MAXITER,niter;
	double fsval,eps;
	double *dx;

	dx = (double*) malloc(sizeof(double) * N);

	fsval = 1.0;
	MAXITER = 200*N;
	niter = 0;
	eps = 1.0e-15; // Use macheps program

	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}

	retval = nel_min(funcpt,xi,N,dx,fsval,MAXITER,&niter,eps,xf);

	//printf("Iterations %d \n", niter);

	free(dx);
	return retval;
}

static int mvalue(int N) {
	int mval;

	if (N <= 10) {
		mval = N;
	} else if (N > 10 && N <= 20) {
		mval = 10;
	} else if (N > 20 && N <= 200) {
		mval = 15;
	} else if ( N > 200) {
		mval = 20;
	}

	return mval;
}

int fminunc(double (*funcpt)(double *,int),int N,double *xi,int method,double *xf) {
	int i,retval,MAXITER,niter,m;
	double fsval,eps,gtol,stol,ftol,xtol,delta;
	double *dx;

	dx = (double*) malloc(sizeof(double) * N);
		/*
		 * Method 0 - Nelder-Mead
		 * Method 1 - Newton Line Search
		 * Method 2 - Newton Trust Region - Hook Step
		 * Method 3 - Newton Trust Region - Double Dog-Leg
		 * Method 4 - Conjugate Gradient
		 * Method 5 - BFGS
		 * Method 6 - Limited Memory BFGS
		 */
	fsval = 1.0;
	MAXITER = 200*N;
	niter = 0;
	delta = -1.0; // Trust Region default
	eps = macheps(); // Use macheps program

	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}

	if (method == 0) {
		retval = nel_min(funcpt,xi,N,dx,fsval,MAXITER,&niter,eps,xf);
	} else if (method == 1) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_func(funcpt,xi,N,dx,fsval,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 2) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,xi,N,dx,fsval,delta,0,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 3) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,xi,N,dx,fsval,delta,1,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 4) {
		gtol = pow(eps,1.0/3.0);
		ftol = gtol * gtol;
		xtol = gtol * ftol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = conjgrad_min_lin(funcpt,xi,N,dx,MAXITER,&niter,eps,gtol,ftol,xtol,xf);
	} else if (method == 5) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = bfgs_min(funcpt,xi,N,dx,fsval,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 6) {
		gtol = pow(eps,1.0/3.0);
		ftol = gtol * gtol;
		xtol = gtol * ftol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		m = mvalue(N);
		retval = bfgs_l_min(funcpt,xi,N,m,dx,fsval,MAXITER,&niter,eps,gtol,ftol,xtol,xf);
	}


	printf("Iterations %d \n", niter);

	free(dx);
	return retval;
}
