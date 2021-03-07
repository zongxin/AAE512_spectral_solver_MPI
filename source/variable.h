#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>


const int N=5;
const int nl=N+1;
const int ng=N;

const int Nblock=112;
const int Ncv=112;
const int NT_HOLEfa=40; //FAKE!!!!
const int N_HOLEfa=8; //FAKE!!!!
const int N_WALLcv=22;
const int N_WALLfa=24;
const int N_HOLE=5;

const int NT_chunk=3;
const int N_chunk1=16;
const int N_chunk2=24;
const int N_chunk3=16;
const int N_mpifa=6;

double t=0;
double dx=1.0/N/3;
double dt=dx*dx*0.5*0.5*0.5*0.5*5*0.1*20;
double x_yita=Nblock;
double bc =5;
double wall_bc=1;
double hole_bc=5;


//one D variable
double* ug_new   = new double[ng];
double* ug 		 = new double[ng];
double* ul		 = new double[nl];

double* g2l		 = new double[nl*ng];
double* l2g		 = new double[ng*nl];

double* Flux_l   = new double[nl];
double* lxx_flux = new double[ng*nl];
double* lyy_flux = new double[ng*nl];
double* xl       = new double[nl];
double* xg       = new double[ng];


//two D variable
double* ug_new_2d    = new double[ng*ng];
double* ug_2d 		 = new double[ng*ng];
double* ul_2d_x		 = new double[ng*nl];
double* ul_2d_y		 = new double[nl*ng];

double* Flux_l_x     = new double[ng*nl];
double* Flux_l_y     = new double[nl*ng];

double* Flux_lxF 	 = new double[ng*nl];
double* Flux_lxG	 = new double[ng*nl];
double* Flux_lyF 	 = new double[nl*ng];
double* Flux_lyG	 = new double[nl*ng];

double* XX	 		 = new double[ng*ng];
double* YY	 		 = new double[ng*ng];

int ierr;
int my_size, my_rank;
int* my_chunk	;
int my_Ncv;

//void nothing(void){
//	printf("nothing", );
//}
