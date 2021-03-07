
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "set_up.h"
#include <iostream>
#include <fstream>
using namespace std; 
double* P1x					= new double[nl];
double* P2x					= new double[nl];
double* P3x					= new double[nl];
double* P4x					= new double[nl];
double* P1y					= new double[nl];
double* P2y					= new double[nl];
double* P3y					= new double[nl];
double* P4y					= new double[nl];

double* der_P1x				= new double[ng];
double* der_P2x				= new double[ng];
double* der_P3x				= new double[ng];
double* der_P4x				= new double[ng];
double* der_P1y				= new double[ng];
double* der_P2y				= new double[ng];
double* der_P3y				= new double[ng];
double* der_P4y				= new double[ng];	

void get_ARC_eq(double*Px,double*Py,double *x,int n,double theta2,double theta1,double r,double cx,double cy){

	for (int j=0;j<n;j++){
		Px[j]=r*cos(theta1+(theta2-theta1)*x[j])+cx;	
		Py[j]=r*sin(theta1+(theta2-theta1)*x[j])+cy;
	}
}

void get_LINE_eq(double*P,double *x,int n,double x_end ,double x_0){

	for (int j=0;j<n;j++){
		P[j]=x_0+(x_end-x_0)*x[j];	
	}	

}


void Initial_mapping(void){
	//read_test();
	read_mesh();

	for (int nb=0;nb<Nblock;nb++){
		cv2d[nb].Block_shape(nb);
    if (my_rank==0){
	     printf("cv:%d,my theta:begin:%f,end:%f\n",nb,cv2d[nb].my_theta[0],cv2d[nb].my_theta[1] );
    }
		cv2d[nb].mapping(nb);
	}
		printf("DONE mapping\n");	
	delete [] P1x,P2x,P3x,P4x;
	delete [] P1y,P2y,P3y,P4y;	
	delete [] der_P1x,der_P2x,der_P3x,der_P4x;
	delete [] der_P1y,der_P2y,der_P3y,der_P4y;		
}

//get_para_func(void)
//get_metrix
//get J_inv
//get real mesh on g

void Block_2d::mapping(int icv)
{


//===============================================================
 //void get_para_func(void)
//===============================================================
	double temp_revx[nl];
	double temp_revy[nl];

	if (my_HOLE_fa[0]!=-2){
		get_ARC_eq(P1x,P1y,xl,nl,my_theta[1],my_theta[0],my_R,my_c[0],my_c[1]);

		for (int i=0;i<nl;i++){	
			temp_revx[i]=P1x[nl-i-1];
			temp_revy[i]=P1y[nl-i-1];			
		}
		for (int i=0;i<nl;i++){	
			P1x[i]=temp_revx[i];
			P1y[i]=temp_revy[i];			
		}
	}
	else{
		get_LINE_eq(P1x,xl,nl,x_no[1],x_no[0]);
		get_LINE_eq(P1y,xl,nl,y_no[1],y_no[0]);		
	}



	if (my_HOLE_fa[1]!=-2){
		get_ARC_eq(P2x,P2y,xl,nl,my_theta[1],my_theta[0],my_R,my_c[0],my_c[1]);
		for (int i=0;i<nl;i++){	
			temp_revx[i]=P2x[nl-i-1];
			temp_revy[i]=P2y[nl-i-1];			
		}
		for (int i=0;i<nl;i++){	
			P2x[i]=temp_revx[i];
			P2y[i]=temp_revy[i];			
		}		
	}
	else{
		get_LINE_eq(P2x,xl,nl,x_no[2],x_no[1]);
		get_LINE_eq(P2y,xl,nl,y_no[2],y_no[1]);		
	}



	if (my_HOLE_fa[2]!=-2){
		get_ARC_eq(P3x,P3y,xl,nl,my_theta[1],my_theta[0],my_R,my_c[0],my_c[1]);
	}
	else{
		get_LINE_eq(P3x,xl,nl,x_no[2],x_no[3]);
		get_LINE_eq(P3y,xl,nl,y_no[2],y_no[3]);		
	}



	if (my_HOLE_fa[3]!=-2){
		get_ARC_eq(P4x,P4y,xl,nl,my_theta[1],my_theta[0],my_R,my_c[0],my_c[1]);
	}
	else{
		get_LINE_eq(P4x,xl,nl,x_no[3],x_no[0]);
		get_LINE_eq(P4y,xl,nl,y_no[3],y_no[0]);	



	}
			

		dot(lxx_flux,ng,nl,P1x,nl,1,der_P1x) ;
		dot(lxx_flux,ng,nl,P1y,nl,1,der_P1y) ;
		dot(lxx_flux,ng,nl,P2x,nl,1,der_P2x) ;
		dot(lxx_flux,ng,nl,P2y,nl,1,der_P2y) ;
		dot(lxx_flux,ng,nl,P3x,nl,1,der_P3x) ;
		dot(lxx_flux,ng,nl,P3y,nl,1,der_P3y) ;
		dot(lxx_flux,ng,nl,P4x,nl,1,der_P4x) ;
		dot(lxx_flux,ng,nl,P4y,nl,1,der_P4y) ;													




//===============================================================
 //void get_metrix(void)
//===============================================================

 //void get_xX( void)

	for (int i=0;i<nl;i++){					//y	coord		
		for (int j=0;j<ng;j++){				//x coord
			xX[(nl-i-1)*ng+j]=(1-xl[i])*der_P1x[j]
						+xl[i]*der_P3x[j]
						-P4x[i]
						+P2x[i]
						+x_no[0]*(1-xl[i])
						-x_no[1]*(1-xl[i])
						-x_no[2]*xl[i]
						+x_no[3]*xl[i];	
		}
	}	


//void get_yX( void)

	for (int i=0;i<nl;i++){					//y	coord		
		for (int j=0;j<ng;j++){				//x coord
			yX[(nl-i-1)*ng+j]=(1-xl[i])*der_P1y[j]
						+xl[i]*der_P3y[j]
						-P4y[i]
						+P2y[i]
						+y_no[0]*(1-xl[i])
						-y_no[1]*(1-xl[i])
						-y_no[2]*xl[i]
						+y_no[3]*xl[i];	
		}
	}

//void get_xY( void)

	for (int i=0;i<ng;i++){					//y	coord		
		for (int j=0;j<nl;j++){				//x coord
			xY[(ng-i-1)*nl+j]= -P1x[j]
						+P3x[j]
						+(1-xl[j])*der_P4x[i]
						+xl[j]*der_P2x[i]
						+x_no[0]*(1-xl[j])
						+x_no[1]*xl[j]
						-x_no[2]*xl[j]
						-x_no[3]*(1-xl[j]);	

		}
	}	

//void get_yY( void)
	for (int i=0;i<ng;i++){					//y	coord		
		for (int j=0;j<nl;j++){				//x coord
			yY[(ng-i-1)*nl+j]= -P1y[j]
						+P3y[j]
						+(1-xl[j])*der_P4y[i]
						+xl[j]*der_P2y[i]
						+y_no[0]*(1-xl[j])
						+y_no[1]*xl[j]
						-y_no[2]*xl[j]
						-y_no[3]*(1-xl[j]);	
		
	
		}	
	}	

//===============================================================
 //void get_J_inv( void)
//===============================================================


	double* J_xX = new double[ng*ng];
	double* J_xY = new double[ng*ng];
	double* J_yX = new double[ng*ng];
	double* J_yY = new double[ng*ng];


	double* J_P1x					= new double[ng];
	double* J_P2x					= new double[ng];
	double* J_P3x					= new double[ng];
	double* J_P4x					= new double[ng];
	double* J_P1y					= new double[ng];
	double* J_P2y					= new double[ng];
	double* J_P3y					= new double[ng];
	double* J_P4y					= new double[ng];

	dot(l2g,ng,nl,P1x,nl,1,J_P1x) ;	
	dot(l2g,ng,nl,P1y,nl,1,J_P1y) ;	
	dot(l2g,ng,nl,P2x,nl,1,J_P2x) ;	
	dot(l2g,ng,nl,P2y,nl,1,J_P2y) ;	
	dot(l2g,ng,nl,P3x,nl,1,J_P3x) ;	
	dot(l2g,ng,nl,P3y,nl,1,J_P3y) ;	
	dot(l2g,ng,nl,P4x,nl,1,J_P4x) ;	
	dot(l2g,ng,nl,P4y,nl,1,J_P4y) ;						

	for (int i=0;i<ng;i++){					//y	coord		
		for (int j=0;j<ng;j++){				//x coord
			J_xX[(ng-i-1)*ng+j]=(1-xg[i])*der_P1x[j]
						+xg[i]*der_P3x[j]
						-J_P4x[i]
						+J_P2x[i]
						+x_no[0]*(1-xg[i])
						-x_no[1]*(1-xg[i])
						-x_no[2]*xg[i]
						+x_no[3]*xg[i];	
		}
	}	

	for (int i=0;i<ng;i++){					//y	coord		
		for (int j=0;j<ng;j++){				//x coord
			J_yX[(ng-i-1)*ng+j]=(1-xg[i])*der_P1y[j]
						+xg[i]*der_P3y[j]
						-J_P4y[i]
						+J_P2y[i]
						+y_no[0]*(1-xg[i])
						-y_no[1]*(1-xg[i])
						-y_no[2]*xg[i]
						+y_no[3]*xg[i];	
		}
	}

	for (int i=0;i<ng;i++){					//y	coord		
		for (int j=0;j<ng;j++){				//x coord
			J_xY[(ng-i-1)*ng+j]= -J_P1x[j]
						+J_P3x[j]
						+(1-xg[j])*der_P4x[i]
						+xg[j]*der_P2x[i]
						+x_no[0]*(1-xg[j])
						+x_no[1]*xg[j]
						-x_no[2]*xg[j]
						-x_no[3]*(1-xg[j]);	
		}
	}	

	for (int i=0;i<ng;i++){					//y	coord		
		for (int j=0;j<ng;j++){				//x coord
			J_yY[(ng-i-1)*ng+j]= -J_P1y[j]
						+J_P3y[j]
						+(1-xg[j])*der_P4y[i]
						+xg[j]*der_P2y[i]
						+y_no[0]*(1-xg[j])
						+y_no[1]*xg[j]
						-y_no[2]*xg[j]
						-y_no[3]*(1-xg[j]);	

		}
	}




	for (int i=0;i<ng;i++){					
		for (int j=0;j<ng;j++){
		J_inv[i*ng+j]=1.0/	(J_xX[i*ng+j]*J_yY[i*ng+j]-J_xY[i*ng+j]*J_yX[i*ng+j]);

		}
	}


//===============================================================
 //void get_mesh_cord(void)
//===============================================================

	double* xxl 				 = new double[nl*nl];//real coordinate		
	double* yyl 				 = new double[nl*nl];
	double* Ttemp_nlng			 = new double[nl*ng];
	double* temp_nlng			 = new double[nl*ng];

	for (int i=0;i<nl;i++){					//y	coord		
		for (int j=0;j<nl;j++){				//x coord
			xxl[(nl-i-1)*nl+j]=(1-xl[i])*P1x[j]
						+xl[i]*P3x[j]
						+(1-xl[j])*P4x[i]
						+xl[j]*P2x[i]
						-x_no[0]*(1-xl[j])*(1-xl[i])
						-x_no[1]*xl[j]*(1-xl[i])
						-x_no[2]*xl[j]*xl[i]
						-x_no[3]*(1-xl[j])*xl[i];	
			yyl[(nl-i-1)*nl+j]=(1-xl[i])*P1y[j]
						+xl[i]*P3y[j]
						+(1-xl[j])*P4y[i]
						+xl[j]*P2y[i]
						-y_no[0]*(1-xl[j])*(1-xl[i])
						-y_no[1]*xl[j]*(1-xl[i])
						-y_no[2]*xl[j]*xl[i]
						-y_no[3]*(1-xl[j])*xl[i];	

		}
	}	

		
	//ul_x=ug.dot(g2l.T)
	transposition(l2g,ng,nl,Ttemp_nlng) ;
	dot(xxl,nl,nl,Ttemp_nlng,nl,ng,temp_nlng) ;	
	dot(l2g,ng,nl,temp_nlng,nl,ng,xxg) ;
	dot(yyl,nl,nl,Ttemp_nlng,nl,ng,temp_nlng) ;	
	dot(l2g,ng,nl,temp_nlng,nl,ng,yyg) ;		




}






 // mapping F,G
//......................................
void Flux_first_map(double *F_on_x,double *G_on_y,double *tempx1_ngng,double *tempy1_ngng,int nb){	
	//Fx=  FX*Xx+FY*Yx
	//Fx= (FX*yY-FY*yX)/J
	//Fx= d(F_on_x*yY)dX-d(G_on_y*yX)dY
	//Fx= A-B
	//Fy=  FX*Xy+FY*Yy
	//Fy=(-FX*xY+FY*xX)/J	
	//Fy= -d(F_on_x*xY)dX+d(G_on_y*xX)dY
	//Fy=-C+D

	double* F_bar  = new double[ng*ng];
	double* G_bar  = new double[ng*ng];

	double* A  = new double[ng*nl]; //x
	double* C  = new double[ng*nl];	//x

	double* B  = new double[nl*ng];	//Y
	double* D  = new double[nl*ng];	//Y

	double* A_temp_ngnl =new double[ng*nl];	
	double* C_temp_ngnl =new double[ng*nl];		
	double* B_temp_nlng =new double[nl*ng];	
	double* D_temp_nlng =new double[nl*ng];	

	double* Ttemp_nlng  = new double[nl*ng];	
	//py:temp=lxx_flux.dot(Flux_l) temp[ng,1]:first order derivative
	// temp2=temp/x_yita
	transposition(lxx_flux,ng,nl,Ttemp_nlng) ;

	//A
	vec_product(Flux_l_x,cv2d[nb].yY,ng,nl,A_temp_ngnl) ;
	//C
	vec_product(Flux_l_x,cv2d[nb].xY,ng,nl,C_temp_ngnl)	;
	//B
	vec_product(Flux_l_y,cv2d[nb].yX,nl,ng,B_temp_nlng) ;
	//D
	vec_product(Flux_l_y,cv2d[nb].xX,ng,nl,D_temp_nlng)	;



	//scalar_profuct(D_temp_nlng,nl,ng,-1,D_temp_nlng);	

	//F_bar=A-B
	dot(A_temp_ngnl,ng,nl,Ttemp_nlng,nl,ng,tempx1_ngng) ;	
	dot(lyy_flux,ng,nl,B_temp_nlng,nl,ng,tempy1_ngng) ;
	scalar_profuct(tempy1_ngng,ng,ng,-1,tempy1_ngng);	

	vec_minus(tempx1_ngng,tempy1_ngng,ng,ng,F_bar);
	//..................................................................................		
//..................................................................................	


//..................................................................................		
//..................................................................................


	//G_bar=-C+D
	dot(C_temp_ngnl,ng,nl,Ttemp_nlng,nl,ng,tempx1_ngng) ;	
	dot(lyy_flux,ng,nl,D_temp_nlng,nl,ng,tempy1_ngng) ;
	scalar_profuct(tempy1_ngng,ng,ng,-1,tempy1_ngng);		
	vec_minus(tempy1_ngng,tempx1_ngng,ng,ng,G_bar) ;


	
	get_value(tempx1_ngng,F_bar,ng,nl) ;
	get_value(tempy1_ngng,G_bar,nl,ng) ;
		
}//Flux_map


void Flux_map(double *F_on_x, double *G_on_x,double *F_on_y, double *G_on_y,int nb){	

	double* F_bar  = new double[ng*nl];
	double* G_bar  = new double[nl*ng];


	for (int i=0;i<ng;i++){					
		for (int j=0;j<nl;j++){
		F_bar[i*nl+j]= cv2d[nb].yY[i*nl+j]*F_on_x[i*nl+j]-cv2d[nb].xY[i*nl+j]*G_on_x[i*nl+j];
		}
	}


	for (int i=0;i<nl;i++){					
		for (int j=0;j<ng;j++){
		G_bar[i*ng+j]=-cv2d[nb].yX[i*ng+j]*F_on_y[i*ng+j]+cv2d[nb].xX[i*ng+j]*G_on_y[i*ng+j];
		}
	}
	get_value(F_on_x,F_bar,ng,nl) ;
	get_value(G_on_y,G_bar,nl,ng) ;
			
}//Flux_map

void div_Jacobian(double *F, double *G,int nb){	

	for (int i=0;i<ng;i++){					
		for (int j=0;j<ng;j++){
		F[i*ng+j]= cv2d[nb].J_inv[i*ng+j]*F[i*ng+j];
		G[i*ng+j]= cv2d[nb].J_inv[i*ng+j]*G[i*ng+j];		
		}
	}
			
}//div_Jacobian



