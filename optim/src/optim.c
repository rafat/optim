/*
 ============================================================================
 Name        : optim.c
 Author      : Rafat Hussain
 Version     :
 Copyright   : 
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "optimize.h"
#include "testfunctions.h"

void printGradient(void(*funcgrad) (double *,int ,double *),double *xi,int N) {

	double *jac;

	jac = (double*) malloc(sizeof(double) * N);

	if (funcgrad == NULL) {
		printf("\n No gradient information \n");
	} else {
		funcgrad(xi,N,jac);
		mdisplay(jac,1,N);
	}

	free(jac);
}

int main(void) {
	int N,j,i,retval,method;
	double *xi,*xf;
	double a,b,oup;
	opt_object optim;

	N = 5;
	a = 0.3;
	b = 1;

	xi = (double*) malloc(sizeof(double) * N);
	xf = (double*) malloc(sizeof(double) * N);

	for (i = 0; i < N;++i) {
		xi[i] = 1;
	}
	xi[0] = 0;
	//xi[1] = 1; xi[2] = 0; xi[3] = 1;
	//xi[9] = -21.4;
	method = 6;
	oup = fminbnd(humps,a,b);
	printf("OUP %g \n",oup);
	printf("%-25s%-20s%-20s \n","Method","Return Value","Objective Function");

	optim = opt_init(N);

	for(j = 0; j < 7;++j) {
		//retval = fminunc(myvalue,myvaluegrad,N,xi,j,xf);
		//printf("%-25d%-20d%-20g \n",j,retval,myvalue(xf,N));
		//mdisplay(xf,1,N);
		optimize(optim,myvalue,myvaluegrad,N,xi,j);
		summary(optim);
		for (i = 0; i < N;++i) {
			xi[i] = 1;
		}
		xi[0] = -1.2;
	}

	//retval = fminunc(brown,browngrad,N,xi,method,xf);

	//printf("Return Value %d Objective Function %g \n",retval,powell(xf,N));

	//printf("Function minimized at : ");
	//mdisplay(xf,1,N);

	//printf("FDX %g \n", (double) FDX);
	//printGradient(powellgrad,xf,N);

	free(xi);
	free(xf);

	return 0;
}
