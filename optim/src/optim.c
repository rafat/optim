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
	retval = fminunc(func1,N,xi,method,xf);

	printf("Return Value %d Objective Function %g \n",retval,func1(xf,N));

	printf("Function minimized at : ");
	mdisplay(xf,1,N);

	//printf("FDX %g \n", (double) FDX);


	free(xi);
	free(xf);

	return 0;
}
