/*
 * testfunctions.c
 *
 *  Created on: Mar 16, 2014
 *      Author: HOME
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "testfunctions.h"


double myvalue
(
    double   *x,
    int       n
)
{
    double f, t ;
    int i ;
    f = 0. ;
    for (i = 0; i < n; i++)
    {
        //t = i+1 ;
        //t = sqrt (t) ;
        //f += exp (x [i]) - t*x [i] ;
    	f += x[i] * x[i] * x[i] * x[i];
    }
    return (f) ;
}

void myvaluegrad(
		double *x,
		int n,
		double *jac
)
{
	int i;
	for(i = 0; i < n; ++i) {
		jac[i] = 4 * x[i] * x[i] * x[i];
	}
}

double quartic(double *x,int N) {
	double f;
	// Powell's Quartic Function
	f = (x[0] + 10 * x[1] ) * (x[0] + 10 * x[1] ) + 5 * (x[2] - x[3]) * (x[2] - x[3])
	+ (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2])
	+ 10 * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) ;
	return f;
}

void quarticgrad(double *x,int N,double *g) {
	double t1,t2,t3,t4;
	// Gradient of Powell's Quartic Function
	t1 = (x[0] + 10 * x[1]);
	t2 = (x[2] - x[3]);
	t3 = (x[1] - 2 * x[2]);
	t4 = (x[0] - x[3]);

	g[0] = 2*t1 + 40*t4*t4*t4;
	g[1] = 20*t1 + 4*t3*t3*t3;
	g[2] = 10*t2 - 8*t3*t3*t3;
	g[3] = -10*t2 - 40*t4*t4*t4;

}

double rosenbrock(double *x,int N) {
	double f,alpha,alpha2;
	alpha = 1;
	alpha2 = alpha*alpha;
	f = 100 * (x[0]*x[0]*alpha2 - x[1]/alpha)* (x[0]*x[0]*alpha2 - x[1]/alpha) + (1 - x[0]*alpha) * (1 - x[0]*alpha);

	return f;
}

void rosenbrockgrad(double *x,int N,double *g) {
	double t,alpha,alpha2;
	alpha = 1;
	alpha2 = alpha * alpha;

	t = (x[0]*x[0]*alpha2 - x[1]/alpha);

	g[0] = 400*t*alpha2*x[0] - 2 * alpha * (1 - x[0]*alpha);
	g[1] = -200*t/alpha;
}

double func4(double *x,int N) {
	double f;
	f = pow((x[0]-2.0),4.0) + pow((x[0]-2.0),2.0) * x[1]*x[1] + (x[1] + 1.0) * (x[1] + 1.0);


	return f;
}

double func1(double *x,int N)
{
	double f;
	double pi,alpha,alpha2;
	pi = 3.14159;
	alpha = 1;
	alpha2 = alpha*alpha;

	//f = pow((x[0]-2.0),4.0) + pow((x[0]-2.0),2.0) * x[1]*x[1] + (x[1] + 1.0) * (x[1] + 1.0);
	//f = 100 * (x[0]*x[0]*alpha2 - x[1]/alpha)* (x[0]*x[0]*alpha2 - x[1]/alpha) + (1 - x[0]*alpha) * (1 - x[0]*alpha);
	//f = (x[0] + 2 * x[1] - 7) * (x[0] + 2 * x[1] - 7) + (2*x[0] + x[1] - 5) * (2*x[0] + x[1] - 5);
/*
	f = (x[0] + 10 * x[1] ) * (x[0] + 10 * x[1] ) + 5 * (x[2] - x[3]) * (x[2] - x[3])
	+ (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2])
	+ 10 * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) ;
*/
	/*
	f = (1.0 /(1.0 + (x[0] - x[1]) * (x[0] - x[1]))) + sin(pi *x[1] * x[2] / 2.0)
	+ exp(-(((x[0]+x[2])/x[1]) - 2) * (((x[0] +x[2])/x[1]) - 2));
	f = -1.0 * f;
	 */

	//f = (10000*x[0]*x[1] - 1) * (10000*x[0]*x[1] - 1) + (exp(-x[0]) + exp(-x[1]) - 1.0001) * (exp(-x[0]) + exp(-x[1]) - 1.0001);
	//f = log(f);

	//f = (x[0] - 1e06)*(x[0] - 1e06) + (x[1] - 2*1e-06)*(x[1] - 2*1e-06) +(x[0]*x[1] - 2)*(x[0]*x[1] - 2);
	f = 100 * (x[1]-x[0]*x[0]) * (x[1]-x[0]*x[0]) + ( 1.0 - x[0] ) * ( 1.0 - x[0] )  + 90 *(x[3]-x[2]*x[2])*(x[3]-x[2]*x[2])
	+ ( 1.0 - x[2])*( 1.0 - x[2]) + 10 * (x[1] + x[3] - 2)*(x[1] + x[3] - 2) + 0.1 * (x[1] - x[3]) *(x[1] - x[3]);
	return f;
}
