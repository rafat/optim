/*
 * optimize.c
 *
 *  Created on: Mar 16, 2014
 *      Author: Rafat Hussain
 */

#include "optimize.h"

opt_object opt_init(int N) {
	opt_object obj = NULL;
	int i;
	double meps;

	obj = (opt_object) malloc (sizeof(struct opt_set) + sizeof(double)* (N-1));

	meps = macheps();
	obj->eps = meps;
	obj->xtol = meps;
	obj->gtol = pow(meps,(double)1.0/3.0);
	obj->ftol = obj->gtol * obj->gtol;
	obj->stol = obj->ftol;
	obj->N = N;
	obj->MaxIter = 200*N;
	obj->retval = 0;

	if (obj->MaxIter < 1000) {
		obj->MaxIter = 1000;
	}

	obj->Iter = 0;
	obj->Method = 0;
	strcpy(obj->MethodName,"Nelder-Mead");
	obj->objfunc = 0.0;
	for (i = 0; i < N;++i) {
		obj->xopt[i] = 0.0;
	}

	return obj;
}

void summary(opt_object obj) {
	int i;
	printf("\n Return Value : %d \n",obj->retval);
	printf("Method : %d %s \n",obj->Method,obj->MethodName);
	printf("Iterations : %d \n",obj->Iter);
	printf("Function Minimized At : [ ");
	for(i = 0; i < obj->N;++i) {
		printf("%g ",obj->xopt[i]);
	}
	printf(" ] \n");
	printf("Function Value : %g \n \n",obj->objfunc);

}

int fminsearch(double (*funcpt)(double *,int),int N,double *xi,double *xf) {
	int i,retval,MAXITER,niter;
	double fsval,eps;
	double *dx;

	dx = (double*) malloc(sizeof(double) * N);

	fsval = 1.0;
	MAXITER = 200*N;
	niter = 0;
	eps = macheps(); // Use macheps program

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

int fminunc(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,int method,double *xf) {
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
		retval = newton_min_func(funcpt,funcgrad,xi,N,dx,fsval,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 2) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,funcgrad,xi,N,dx,fsval,delta,0,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 3) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,funcgrad,xi,N,dx,fsval,delta,1,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 4) {
		gtol = pow(eps,1.0/3.0);
		ftol = gtol * gtol;
		xtol = eps;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = conjgrad_min_lin(funcpt,funcgrad,xi,N,dx,MAXITER,&niter,eps,gtol,ftol,xtol,xf);
	} else if (method == 5) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = bfgs_min(funcpt,funcgrad,xi,N,dx,fsval,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 6) {
		gtol = pow(eps,1.0/3.0);
		ftol = gtol * gtol;
		xtol = eps;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		m = mvalue(N);
		retval = bfgs_l_min(funcpt,funcgrad,xi,N,m,dx,fsval,MAXITER,&niter,eps,gtol,ftol,xtol,xf);
	} else {
		printf("Method Value should be one of 0,1,2,3,4,5 or 6. See Documentation. \n");
		exit(1);
	}


	//printf("Iterations %d \n", niter);

	free(dx);
	return retval;
}

double fminbnd(double (*funcpt)(double),double a, double b) {
	double x,t,eps;

	t = 1e-012;
	eps = macheps();

	brent_local_min(funcpt,a,b,t,eps,&x);

	return x;
}

int fminqn(double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,
		double delta,int method,double *dx,double fsval,double *xf) {
	int retval,i;
	int MAXITER,niter;
	double eps,gtol,stol;
	/*

	 * Method 1 - Newton Line Search
	 * Method 2 - Newton Trust Region - Hook Step
	 * Method 3 - Newton Trust Region - Double Dog-Leg
	 *
	 * Default Values :
	 *
	 * fsval = 1.0
	 * delta = -1.0
	 * dx = {1.0,1.0,...} - 1XN vector

	 */

	MAXITER = 200*N;
	niter = 0;

	eps = macheps(); // Use macheps program

	if (method == 1) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_func(funcpt,funcgrad,xi,N,dx,fsval,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 2) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,funcgrad,xi,N,dx,fsval,delta,0,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 3) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,funcgrad,xi,N,dx,fsval,delta,1,MAXITER,&niter,eps,gtol,stol,xf);
	} else {
		printf("Method Value should be one of 1,2 or 3. See Documentation. \n");
		exit(1);
	}


	return retval;
}

void optimize(opt_object obj,double (*funcpt)(double *,int),void(*funcgrad)(double *, int,double *),int N,double *xi,
		int method) {
	int i,m;
	double fsval,delta;
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
	delta = -1.0; // Trust Region default


	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}
	obj->Method = method;
	obj->Iter = 0;

	if (obj->N != N) {
		printf("The Object is initialized for a problem of size %d \n",obj->N);
		printf("Please Reinitialize the object for the new size %d using opt_init command \n",N);
		exit(1);
	}

	if (method == 0) {
		strcpy(obj->MethodName,"Nelder-Mead");
		obj->retval = nel_min(funcpt,xi,obj->N,dx,fsval,obj->MaxIter,&obj->Iter,obj->eps,obj->xopt);
	} else if (method == 1) {
		strcpy(obj->MethodName,"Newton Line Search");
		obj->retval = newton_min_func(funcpt,funcgrad,xi,obj->N,dx,fsval,obj->MaxIter,&obj->Iter,obj->eps,obj->gtol,
				obj->stol,obj->xopt);
	} else if (method == 2) {
		strcpy(obj->MethodName,"Newton Trust Region - Hook Step");

		obj->retval = newton_min_trust(funcpt,funcgrad,xi,obj->N,dx,fsval,delta,0,obj->MaxIter,&obj->Iter,obj->eps,
				obj->gtol,obj->stol,obj->xopt);
	} else if (method == 3) {
		strcpy(obj->MethodName,"Newton Trust Region - Double Dog-Leg");

		obj->retval = newton_min_trust(funcpt,funcgrad,xi,obj->N,dx,fsval,delta,1,obj->MaxIter,&obj->Iter,obj->eps,
				obj->gtol,obj->stol,obj->xopt);
	} else if (method == 4) {
		strcpy(obj->MethodName,"Conjugate Gradient");

		obj->retval = conjgrad_min_lin(funcpt,funcgrad,xi,obj->N,dx,obj->MaxIter,&obj->Iter,obj->eps,obj->gtol,
				obj->ftol,obj->xtol,obj->xopt);
	} else if (method == 5) {
		strcpy(obj->MethodName,"BFGS");

		obj->retval = bfgs_min(funcpt,funcgrad,xi,obj->N,dx,fsval,obj->MaxIter,&obj->Iter,obj->eps,obj->gtol,
				obj->stol,obj->xopt);
	} else if (method == 6) {
		strcpy(obj->MethodName,"Newton Line Search");
		m = mvalue(N);
		obj->retval = bfgs_l_min(funcpt,funcgrad,xi,obj->N,m,dx,fsval,obj->MaxIter,&obj->Iter,obj->eps,
				obj->gtol,obj->ftol,obj->xtol,obj->xopt);
	} else {
		strcpy(obj->MethodName,"NULL");
		printf("Method Value should be one of 0,1,2,3,4,5 or 6. See Documentation. \n");
		exit(1);
	}

	obj->objfunc = funcpt(obj->xopt,obj->N);

	free(dx);
}
