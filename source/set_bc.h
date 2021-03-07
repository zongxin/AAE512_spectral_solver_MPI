#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "MPI_utilities.h"
void Initial_bc(int nb){

	//int l=nb-1,r=nb+1,u=nb-N_col,d=nb+N_col;
	int l,r,u,d;
	int this_fa;
	int corr;
  int INTER;

	for (int i=0;i<4;i++){

		this_fa=cv2d[nb].my_fa[i];
    for (int j=0;j<N_mpifa;j++){
    if (nb==45){
        //cout<<"Nb:"<<nb<<"falist1 "<<fa1_list[j]<<endl;
        //cout<<"Nb:"<<nb<<"falist2 "<<fa2_list[j]<<endl;
     }

      if ((this_fa==fa1_list[j])||(this_fa==fa2_list[j])){
          INTER=1;
          break;
      }
      else{
          INTER=0;
      }

    }
  



		if(i==3){

			if (cv2d[nb].my_nei[i]!=-1){

				l=cv2d[nb].my_nei[i];

				for (int j=0;j<4;j++){
					if(cv2d[l].my_fa[j]==this_fa){
						corr=j;
					}
				}//end corr	

				for (int i=0;i<ng;i++){	
					if (corr==2){
						ul_2d_x[nl*i+0]=(cv2d[l].ul_y_class[0*ng+i]+ul_2d_x[nl*i+0])/2.0;								
					}		
					else if (corr==0){
						ul_2d_x[nl*i+0]=(cv2d[l].ul_y_class[(nl-1)*ng+ng-1-i]+ul_2d_x[nl*i+0])/2.0;						
					}			
					else{

            if (INTER==1){
			      	ul_2d_x[nl*i+0]=(cv2d[nb].my_recv[i]+ul_2d_x[nl*i+0])/2.0;	

            }
            else{
			      	ul_2d_x[nl*i+0]=(cv2d[l].ul_x_class[nl*i+nl-1]+ul_2d_x[nl*i+0])/2.0;		
            }
            						
					}									

				}//end 7 or not

			}//end interface
			
			else if (cv2d[nb].my_WALL_fa[i]!=-2)	{						
				for (int i=0;i<ng;i++){	
					ul_2d_x[nl*i+0]=(wall_bc+ul_2d_x[nl*i+0])*0.5;
					//ul_2d_x[nl*i+0]=(3+ul_2d_x[nl*i+0])*0.5;
				}
			}		
			else if (cv2d[nb].my_HOLE_fa[i]!=-2)	{
				for (int i=0;i<ng;i++){	
					ul_2d_x[nl*i+0]=(hole_bc+ul_2d_x[nl*i+0])*0.5;
				}
			}		
			else{ 
				printf("COULD NOT MATCH BC\n");
				printf("nb:%d,fa:%d,my_nei:%d,wall:%d,hole:%d\n",
					nb,this_fa,cv2d[nb].my_nei[i],cv2d[nb].my_WALL_fa[i],cv2d[nb].my_HOLE_fa[i]);							
			}	
													
		}// end i3


		if(i==1){
			if (cv2d[nb].my_nei[i]!=-1){

				r=cv2d[nb].my_nei[i];

				for (int j=0;j<4;j++){
					if(cv2d[r].my_fa[j]==this_fa){
						corr=j;
					}
				}	

				for (int i=0;i<ng;i++){	

					if (corr==2){
						ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_y_class[0*ng+ng-1-i]+ul_2d_x[nl*i+nl-1])/2.0;	
						//zongxin here
						//ul_2d_x[nl*i+nl-1]=i*10;

					}
					else if (corr==0){
						ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_y_class[(nl-1)*ng+i]+ul_2d_x[nl*i+nl-1])/2.0;

					}
					else{

            if (INTER==1){
			      	ul_2d_x[nl*i+nl-1]=(cv2d[nb].my_recv[i]+ul_2d_x[nl*i+nl-1])/2.0;	

            }
            else{
			      	ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_x_class[nl*i+0]+ul_2d_x[nl*i+nl-1])/2.0;		
            }
						
					}	

				}



			}//end interface
			
			else if (cv2d[nb].my_WALL_fa[i]!=-2)	{
				for (int i=0;i<ng;i++){	
					ul_2d_x[nl*i+nl-1]=(ul_2d_x[nl*i+nl-1]+wall_bc)*0.5;
					//ul_2d_x[nl*i+nl-1]=(ul_2d_x[nl*i+nl-1]+3)*0.5;	
				}
			}		
			else if (cv2d[nb].my_HOLE_fa[i]!=-2)	{
				for (int i=0;i<ng;i++){	
					ul_2d_x[nl*i+nl-1]=(ul_2d_x[nl*i+nl-1]+hole_bc)*0.5;								
				}	
			}										
			else{
			printf("COULD NOT MATCH BC\n");
			printf("nb:%d,fa:%d,my_nei:%d,wall:%d,hole:%d\n",
				nb,this_fa,cv2d[nb].my_nei[i],cv2d[nb].my_WALL_fa[i],cv2d[nb].my_HOLE_fa[i]);						
			}		
												
		}//end i1

		if(i==2){

			if (cv2d[nb].my_nei[i]!=-1){

				u=cv2d[nb].my_nei[i];

				for (int j=0;j<4;j++){
					if(cv2d[u].my_fa[j]==this_fa){
						corr=j;
					}
				}	

				for (int j=0;j<ng;j++){	
					if (corr==1){
						ul_2d_y[0*ng+j]=(cv2d[u].ul_x_class[nl*(ng-1-j)+nl-1]+ul_2d_y[0*ng+j])/2.0;	
						//zongxin here
						//ul_2d_y[0*ng+j]=(ng-1-i)*10;								
					}
					else if(corr==3){
						ul_2d_y[0*ng+j]=(cv2d[u].ul_x_class[nl*j+0]+ul_2d_y[0*ng+j])/2.0;
					}		
					else{
						ul_2d_y[0*ng+j]=(cv2d[u].ul_y_class[(nl-1)*ng+j]+ul_2d_y[0*ng+j])/2.0;		//									
					}						

				}

			}//end interface
		
			else if (cv2d[nb].my_WALL_fa[i]!=-2)	{
				for (int j=0;j<ng;j++){								
					ul_2d_y[0*ng+j]=(wall_bc+ul_2d_y[0*ng+j])*0.5;	

					//ul_2d_y[0*ng+j]=(2+ul_2d_y[0*ng+j])*0.5;	
				}
			}		
			else if (cv2d[nb].my_HOLE_fa[i]!=-2)	{
				for (int j=0;j<ng;j++){	
					ul_2d_y[0*ng+j]=(hole_bc+ul_2d_y[0*ng+j])*0.5;		
				}
			}										
			else{
			printf("COULD NOT MATCH BC\n");
			printf("nb:%d,fa:%d,my_nei:%d,wall:%d,hole:%d\n",
				nb,this_fa,cv2d[nb].my_nei[i],cv2d[nb].my_WALL_fa[i],cv2d[nb].my_HOLE_fa[i]);						
			}	

		}//end i2

		if(i==0){
			if (cv2d[nb].my_nei[i]!=-1){

				d=cv2d[nb].my_nei[i];

				for (int j=0;j<4;j++){
					if(cv2d[d].my_fa[j]==this_fa){
						corr=j;
					}
				}	

				for (int j=0;j<ng;j++){	
					if (corr==1){
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_x_class[nl*j+nl-1]+ul_2d_y[(nl-1)*ng+j])/2.0;
					}	

					else if (corr==3){
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_x_class[nl*(ng-1-j)+0]+ul_2d_y[(nl-1)*ng+j])/2.0;	
					}			
					else{
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_y_class[0*ng+j]+ul_2d_y[(nl-1)*ng+j])/2.0;	
					}
				}


			}// end interface
	
			else if (cv2d[nb].my_WALL_fa[i]!=-2)	{
				for (int j=0;j<ng;j++){	
					ul_2d_y[(nl-1)*ng+j]=(ul_2d_y[(nl-1)*ng+j]+wall_bc)*0.5;	
					//ul_2d_y[(nl-1)*ng+j]=(ul_2d_y[(nl-1)*ng+j]+0)*0.5;	
				}
			}		
			else if (cv2d[nb].my_HOLE_fa[i]!=-2)	{
				for (int j=0;j<ng;j++){	
					ul_2d_y[(nl-1)*ng+j]=(ul_2d_y[(nl-1)*ng+j]+hole_bc)*0.5;	
				}	
			}										
			else{
			printf("COULD NOT MATCH BC\n");
			printf("nb:%d,fa:%d,my_nei:%d,wall:%d,hole:%d\n",
				nb,this_fa,cv2d[nb].my_nei[i],cv2d[nb].my_WALL_fa[i],cv2d[nb].my_HOLE_fa[i]);						
			}	
	
		}//end i0	

	}	
}












/*
                          if (nb==37){
                              cout<<"Nb:"<<nb<<" this fa:"<<this_fa<<"inter:"<<INTER<<endl;
                           }
*/








void Gradient_bc(int nb){
	int l,r,u,d;
	int this_fa;
	int corr;
  int INTER;
	for (int i=0;i<4;i++){

		this_fa=cv2d[nb].my_fa[i];
    for (int j=0;j<N_mpifa;j++){
      if ((this_fa==fa1_list[j])||(this_fa==fa2_list[j])){
          INTER=1;
          break;
      }
      else{
          INTER=0;
      }
    }
		if(i==3){
			if (cv2d[nb].my_nei[i]!=-1){

				l=cv2d[nb].my_nei[i];

				for (int j=0;j<4;j++){
					if(cv2d[l].my_fa[j]==this_fa){
						corr=j;
					}
				}	

								
				for (int i=0;i<ng;i++){	
					if (corr==2){
						Flux_lxF[nl*i+0]=(cv2d[l].Flux_ly_class_F[0*ng+i]+Flux_lxF[nl*i+0])/2.0;
						Flux_lxG[nl*i+0]=(cv2d[l].Flux_ly_class_G[0*ng+i]+Flux_lxG[nl*i+0])/2.0;								
					}		
					else if (corr==0){
						Flux_lxF[nl*i+0]=(cv2d[l].Flux_ly_class_F[(nl-1)*ng+ng-1-i]+Flux_lxF[nl*i+0])/2.0;
						Flux_lxG[nl*i+0]=(cv2d[l].Flux_ly_class_G[(nl-1)*ng+ng-1-i]+Flux_lxG[nl*i+0])/2.0;							
					}			
					else{

            if (INTER==1){	
						  Flux_lxF[nl*i+0]=(cv2d[nb].my_recv_F[i]+Flux_lxF[nl*i+0])/2.0;			//
						  Flux_lxG[nl*i+0]=(cv2d[nb].my_recv_G[i]+Flux_lxG[nl*i+0])/2.0;		
            }
            else{
						  Flux_lxF[nl*i+0]=(cv2d[l].Flux_lx_class_F[nl*i+nl-1]+Flux_lxF[nl*i+0])/2.0;			//
						  Flux_lxG[nl*i+0]=(cv2d[l].Flux_lx_class_G[nl*i+nl-1]+Flux_lxG[nl*i+0])/2.0;		
            }

								
					}						

				}
			}								
		}//end i3

		if(i==1){
			if (cv2d[nb].my_nei[i]!=-1){

				r=cv2d[nb].my_nei[i];

				for (int j=0;j<4;j++){
					if(cv2d[r].my_fa[j]==this_fa){
						corr=j;
					}
				}	

				for (int i=0;i<ng;i++){	

					if (corr==2){
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_ly_class_F[0*ng+ng-1-i]+Flux_lxF[nl*i+nl-1])/2.0;
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_ly_class_G[0*ng+ng-1-i]+Flux_lxG[nl*i+nl-1])/2.0;								
					}
					else if (corr==0){
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_ly_class_F[(nl-1)*ng+i]+Flux_lxF[nl*i+nl-1])/2.0;
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_ly_class_G[(nl-1)*ng+i]+Flux_lxG[nl*i+nl-1])/2.0;								
					}
					else{

            if (INTER==1){	
						  Flux_lxF[nl*i+nl-1]=(cv2d[nb].my_recv_F[i]+Flux_lxF[nl*i+nl-1])/2.0;			//
						  Flux_lxG[nl*i+nl-1]=(cv2d[nb].my_recv_G[i]+Flux_lxG[nl*i+nl-1])/2.0;		
            }
            else{
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_lx_class_F[nl*i+0]+Flux_lxF[nl*i+nl-1])/2.0;		//
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_lx_class_G[nl*i+0]+Flux_lxG[nl*i+nl-1])/2.0;	
            }

					}	

				}
			}																						
		}//end i1

		if(i==2){
			if (cv2d[nb].my_nei[i]!=-1){

				u=cv2d[nb].my_nei[i];

				for (int j=0;j<4;j++){
					if(cv2d[u].my_fa[j]==this_fa){
						corr=j;
					}
				}	


				for (int j=0;j<ng;j++){	
					if (corr==1){
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_lx_class_F[nl*(ng-1-j)+nl-1]+Flux_lyF[0*ng+j])/2.0;
						Flux_lyG[0*ng+j]=(cv2d[u].Flux_lx_class_G[nl*(ng-1-j)+nl-1]+Flux_lyG[0*ng+j])/2.0;								

					}
					else if(corr==3){
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_lx_class_F[nl*j+0]+Flux_lyF[0*ng+j])/2.0;
						Flux_lyG[0*ng+j]=(cv2d[u].Flux_lx_class_G[nl*j+0]+Flux_lyG[0*ng+j])/2.0;								
					}		
					else{
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_ly_class_F[(nl-1)*ng+j]+Flux_lyF[0*ng+j])/2.0;	
						Flux_lyG[0*ng+j]=(cv2d[u].Flux_ly_class_G[(nl-1)*ng+j]+Flux_lyG[0*ng+j])/2.0;		//									
					}						

				}
			}									
		}// end i2

		if(i==0){
			if (cv2d[nb].my_nei[i]!=-1){

				d=cv2d[nb].my_nei[i];

				for (int j=0;j<4;j++){
					if(cv2d[d].my_fa[j]==this_fa){
						corr=j;
					}
				}			

				for (int j=0;j<ng;j++){	

					if (corr==1){
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_F[nl*j+nl-1]+Flux_lyF[(nl-1)*ng+j])/2.0;
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_G[nl*j+nl-1]+Flux_lyG[(nl-1)*ng+j])/2.0;								
					}	

					else if (corr==3){
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_F[nl*(ng-1-j)+0]+Flux_lyF[(nl-1)*ng+j])/2.0;	
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_G[nl*(ng-1-j)+0]+Flux_lyG[(nl-1)*ng+j])/2.0;									
					}	
					else{
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_ly_class_F[0*ng+j]+Flux_lyF[(nl-1)*ng+j])/2.0;	
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_ly_class_G[0*ng+j]+Flux_lyG[(nl-1)*ng+j])/2.0;	//									
					}									

				}
			}//end interface										
		}//end i2		
	}	
}
