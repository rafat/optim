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
	int N,i,retval,method;
	double *xi,*xf;

	N = 4;

	xi = (double*) malloc(sizeof(double) * N);
	xf = (double*) malloc(sizeof(double) * N);

	for (i = 0; i < N;++i) {
		xi[i] = 3;
	}
	xi[1] = -1; xi[2] = 0; xi[3] = 1;
	//xi[9] = -21.4;
	method = 5;
	retval = fminunc(myvalue,NULL,N,xi,method,xf);
	//retval = fminunc(myvalue,myvaluegrad,N,xi,method,xf);

	printf("Return Value %d Objective Function %g \n",retval,myvalue(xf,N));

	printf("Function minimized at : ");
	mdisplay(xf,1,N);

	//printf("FDX %g \n", (double) FDX);
	printGradient(myvaluegrad,xf,N);

	free(xi);
	free(xf);

	return 0;
}
