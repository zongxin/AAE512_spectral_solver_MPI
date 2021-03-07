#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
//#include "variable.h"
#include "in_put.h"
using namespace std;
class Block{
public:
 	double ratio;
 	double offset;
 	double* ug_new_class;// 		 = new double[ng];
  	double* ul_class;// 		 = new double[nl];
  	double* Flux_lx_class;// 		 = new double[nl];
	
	Block()
	{
		ug_new_class = new double[ng];
		ul_class = new double[nl];
		Flux_lx_class = new double[nl];
	}
	~Block(){

	}
				

};



class Block_2d{
 public:
		
		double* ug_new_class;

		//first flux on l
		double* ul_x_class;
		double* ul_y_class;

		//second flux on l
		double* Flux_lx_class;
		double* Flux_ly_class;

		double* Flux_lx_class_F;
		double* Flux_lx_class_G;
		double* Flux_ly_class_G;
		double* Flux_ly_class_F;

		//Geometry info input
		double  x_no[4];
		double  y_no[4];
		double  x_fa[4];
		double  y_fa[4];

		int	 	my_fa[4];
		int 	my_WALL_fa[4];
		int 	my_HOLE_fa[4];
		int 	my_nei[4];
		double  my_theta[2]; //begin,end
		double  my_c[2];// x, y
		double  my_R;
		
		double* xxg;//real coordinate		
		double* yyg;	

		//mapping
 		double* xX;
 		double* xY;
 		double* yX;
 		double* yY;
 		double* J_inv;



 		double* my_recv;
 		double* my_sent;
 		double* my_sent_F;
 		double* my_sent_G;
 		double* my_recv_F;
 		double* my_recv_G;
 		double* mpi_test ;
 		void Block_shape(int icv);
 		void mapping(int icv);

	Block_2d()
	{
		ug_new_class 		 = new double[ng*ng];

		//first flux on l
		ul_x_class 		     = new double[ng*nl];
		ul_y_class 		     = new double[nl*ng];

		//second flux on l
		Flux_lx_class 		 = new double[ng*nl];
		Flux_ly_class 		 = new double[nl*ng];

		Flux_lx_class_F		 = new double[ng*nl];
		Flux_lx_class_G		 = new double[ng*nl];
		Flux_ly_class_G		 = new double[nl*ng];
		Flux_ly_class_F		 = new double[nl*ng];

		//Geometry info input
		x_no[4];
		y_no[4];
		x_fa[4];
		y_fa[4];

		my_fa[4];
		my_WALL_fa[0]=-2;
		my_WALL_fa[1]=-2;
		my_WALL_fa[2]=-2;
		my_WALL_fa[3]=-2;

		my_HOLE_fa[0]=-2;
		my_HOLE_fa[1]=-2;
		my_HOLE_fa[2]=-2;
		my_HOLE_fa[3]=-2;

		my_nei[4];
		my_theta[0]=10; //begin,end
		my_theta[1]=10; //begin,end
		my_c[2];// x, y
		my_R=1;


		xxg 					 = new double[ng*ng];//real coordinate		
		yyg 					 = new double[ng*ng];	

		//mapping
 		xX = new double[nl*ng];
 		xY = new double[ng*nl];
 		yX = new double[nl*ng];
 		yY = new double[ng*nl];
 		J_inv  = new double[ng*ng];

 		my_recv= new double[ng];
 		my_sent= new double[ng];

 		my_sent_F= new double[ng];
 		my_sent_G= new double[ng];
 		my_recv_F= new double[ng];
 		my_recv_G= new double[ng];

 		mpi_test = new double[ng*nl];
 		//Block_shape();
 		//mapping();
	}

	~Block_2d(){

	}

};

void Block_2d::Block_shape(int icv)
{

	for (int i=0;i<NT_HOLEfa;i++){

		if (icv==HOLE1_cv[i]){
			my_c[0]=c_list[0];
			my_c[1]=c_list[1];
			break;
		}
		if (icv==HOLE2_cv[i]){
			my_c[0]=c_list[2];
			my_c[1]=c_list[3];
			break;
		}	
		if (icv==HOLE3_cv[i]){
			my_c[0]=c_list[4];
			my_c[1]=c_list[5];
			break;
		}	
		if (icv==HOLE4_cv[i]){

			my_c[0]=c_list[6];
			my_c[1]=c_list[7];
			break;
		}	
		if (icv==HOLE5_cv[i]){
			my_c[0]=c_list[8];
			my_c[1]=c_list[9];
			break;
		}								
	}

	for (int i=0;i<4;i++){

		x_no[i]		= nox_ocv[icv*4+i];
		y_no[i]		= noy_ocv[icv*4+i];
		x_fa[i]		= fax_ocv[icv*4+i];
		y_fa[i]		= fax_ocv[icv*4+i];
		my_fa[i]	= faocv  [icv*4+i];
		my_nei[i]	= neiocv [icv*4+i];
	}
	int ii;

	for (int j=0;j<NT_HOLEfa;j++){


		for (int i=0;i<4;i++){

			if (my_fa[i]==HOLE_fa[j]){

				//get wall number
				my_HOLE_fa[i]=my_fa[i];

				//get theta 
				double dx=x_no[i]-my_c[0];
				double dy=y_no[i]-my_c[1];	
			
				if ((dx==0)&&(dy<0)){
					my_theta[1]=M_PI*3/2;					
				}			
				else if ((dx==0)&&(dy>0)){
					my_theta[1]=M_PI*0.5;					
				}	
				else if (dx!=0){
					my_theta[1]=atan(dy/dx);					
				}	




				if (dx<0){
					my_theta[1]=my_theta[1]+M_PI;
				}

				if ((dx>0)&&(dy<=0)){
					my_theta[1]=my_theta[1]+M_PI*2;		
	
				}	
				if (icv==30)	{
        //cout<<"@@@@@inside angle"<<endl;
        //cout<<"dx:"<<dx<<"dy:"<<dy<<"theta"<<my_theta[1]<<endl;
			}
				my_theta[0]=my_theta[1]-M_PI/8.0;															
			}//end ==
		}
		
	}//end hole

	for (int j=0;j<N_WALLfa;j++){
		for (int i=0;i<4;i++){

			if (my_fa[i]==WALL_fa[j]){	
				my_WALL_fa[i]=my_fa[i];
			}
		}
	}//end wall		
	

}

