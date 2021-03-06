
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "vector_tool.h"
//#include "variable.h"
#include "class.h"
#include <iostream>
#include <fstream>
#include "JacobiPoly.h"
#include "lapacke.h"

using namespace std; 
 // interpolate base
//......................................
class Block cv[Nblock];
class Block_2d cv2d[Nblock];


void base(double *xlist,int n,double x,double *h){	
	
	for (int j = 0; j < n; j++) {
		h[j]=1;
		for (int i = 0; i < n; i++) {
			if (i!=j){
				double m=(x-xlist[i])/(xlist[j]-xlist[i]);
				h[j]=h[j]*m;
				}
			}
		}
		
}//base


 // first order derivative of interpolate base
//......................................
void der_base(double *xlist, int n,double x,double *hx){	

	double* h  = new double[n];
	double m;
	base((double*)xlist,n,x,(double*)h);
	
	for (int j = 0; j < n; j++) {
		
		hx[j]=0.0;

		for (int i = 0; i < n; i++) {
			if (i!=j){
				m=1.0/(x-xlist[i]);
				hx[j]=hx[j]+m;
			}
		}
		
		hx[j]=hx[j]*h[j]; 	
	}
	
	
	delete [] h;	
}//der_base



////////////////////////////////////////////////////////////////////
/////////////////////       Set up Begin       /////////////////////
////////////////////////////////////////////////////////////////////

void set_up(void){

	//double *xl = malloc(nl*sizeof(double));
  //double *xg = malloc(ng*sizeof(double));
  //double *h = malloc(ng*sizeof(double));


	double* h  = new double[ng];
	double* hx = new double[ng];
	int ini_sin=0;
	int ini_test=0;
	int ini_cos_2d=0;





	/***********************************************************************************
								BLOCK	 MESH
	***********************************************************************************/


 // Set up 1-D mesh
/* OLD set of points from (Kopriva Kolias , Chebyshev points)
	for (int i = 0; i < nl; i++) {
		xl[i]=0.5*(1-cos(M_PI*i/N));
	}

	for (int i = 0; i < nl; i++) {
		xg[i]=0.5*(1-cos(M_PI*(2*i+1)/(2*N)));
	//	printf("%f  ",xg[i]);
	}
*/

/*new set of points using Prateek's Jacobipoly library*/
	double *J1= get_JACOBI(0.0, 0.0, ng);

	xg = get_EIGS(J1, ng);
	for (int i = 0; i < ng; i++) {
		xg[i]=0.5*(1.0 + xg[i]);
	}

	double *J2 = get_JACOBI_LOBATTO(0.0, 0.0, nl);

	xl = get_EIGS(J2, nl);
	for (int i = 0; i < nl; i++) {
		xl[i]=0.5*(1.0 + xl[i]);
	}

	// 	2-D mesh
	for (int i = 0; i < ng; i++) {
		for (int j = 0; j < ng; j++) {
		//printf("%d  %f \n",i, XXg[i][j]);
			XX[i*ng+j]=xg[j];
			YY[i*ng+j]=xg[i];
		}
	}



	/***********************************************************************************
								INITIAL CONDITION
	***********************************************************************************/


	 // set_up 1-d initial condition 
	if (ini_sin==1){

		for (int nb=0;nb<Nblock;nb++){
			for (int i = 0; i < ng; i++) {
			cv[nb].ug_new_class[i] = cos(2*M_PI/Nblock*(xg[i]+nb)); 
			}
		}

	}

	if (ini_test==1){

		for (int nb=0;nb<Nblock;nb++){
			for (int i = 0; i < ng; i++) {
			cv[nb].ug_new_class[i] = 0.0;
			}
		}

	}

		// set_up 2-d initial condition 

	if (ini_cos_2d==1){

		for (int nb=0;nb<Nblock;nb++){
			for (int i = 0; i < ng; i++) {
				for (int j = 0; j < ng; j++) {
					//cv2d[nb].ug_new_class[i*ng+j] = cos(2*M_PI/Nblock*(xg[i]+nb))*cos(2*M_PI/Nblock*(xg[j]+nb));
					//cv2d[nb].ug_new_class[i*ng+j] = 2*M_PI/Nblock*(xg[j]+nb);
				}//end j
			}//end i
		}//end block

	}//end if 2d cos







	/***********************************************************************************
									 OPERATORS
	***********************************************************************************/


 // Operator g_l[nl][ng]:from g to l 
	double* temp1d_g=new double[ng];
	for (int i = 0; i < nl; i++) {
		base((double*)xg,ng,xl[i],(double*)temp1d_g);
		for (int j = 0; j < ng; j++) {
			g2l[i*ng+j]=temp1d_g[j];
		}
	}
	delete []temp1d_g;

 // Operator l_g[ng][nl]:from l --> g 
	double* temp1d_l=new double[nl];
	for (int i = 0; i < ng; i++) {
		base((double*)xl,nl,xg[i],(double*)temp1d_l);
		for (int j = 0; j < nl; j++) {
			l2g[i*nl+j]=temp1d_l[j];
		}
	}


 // Operator lxx_flux[ng][nl]:from l to g when do first derivative
	
	for (int i = 0; i < ng; i++) {
		der_base((double*)xl,nl,xg[i],(double*)temp1d_l);
		for (int j = 0; j < nl; j++) {
			lxx_flux[i*nl+j]=temp1d_l[j];
		}
	}
	

 // Operator lxx_flux[ng][nl]:from l to g when do first derivative
	
	for (int i = 0; i < ng; i++) {
		der_base((double*)xl,nl,xg[i],(double*)temp1d_l);
		for (int j = 0; j < nl; j++) {
			lyy_flux[i*nl+j]=temp1d_l[j];
		}
	}
	
	delete []temp1d_l;

 // for test
//......................................
	double* y  = new double[nl];
	double* hy  = new double[nl];
	double temp=0;
	//base((double*)xl,nl,1.0,(double*)h);
	der_base((double*)xl,nl,xg[3],(double*)hy);

	for (int i = 0; i < nl; i++) {

		y[i]=xl[i]*xl[i]*xl[i];

		temp=temp+y[i]*hy[i];
	}


	delete []y;





}//set_up







