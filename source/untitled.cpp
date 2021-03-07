

//
//   u_t + u u_x = u_xx
//
// To compile:
//    g++ -Wall -o test untitled.cpp
//
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include<sys/time.h>
#include "set_bc.h"
#include "out_put.h"


using namespace std;

int one_d_wave = 0  ;
int one_d_heat = 0  ;
int two_d_heat = 1  ;

int OUT        = 0  ;
int TEST       = 0  ;
int oo_test = 10;
int oo_file = 1000;

int main(int argc, char *argv[]) {  

  ierr = MPI_Init(&argc, &argv);	
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &my_size);

  struct timeval time1;
  gettimeofday(&time1,NULL);


	t=0;
	set_up();
	Initial_mapping();
	output_ini();
  Initial_MPI();
  double* test   = new double[ng*nl];
	for (int i=0;i<ng;i++){	
  	for (int j=0;j<nl;j++){	
    test[nl*i+j]=i*0.1+j+4; 
  }
}

for (int n=1;n<20000;n++){

		if (one_d_wave == 1){

			double* k   = new double[ng];
			double* temp   = new double[ng];
			double* temp2   = new double[ng];			

			printf("dx=%f",(xl[1]-xl[0])/5 );
			printf("dt=%f",dt/((xl[1]-xl[0])/5)) ;
			printf("\n" );

			// ug=ug_new;	
			// //dot_product((double*)A,m,n,(double*)B,x,y,(double*)C) 
			// dot((double*)g2l,nl,ng,(double*)ug,ng,1,(double*)ul) ;	
			// ul[0]=ul4[nl-1];
			// //ul[nl-1]=0;
			// Flux_l=ul;	
			// //py:temp=lxx_flux.dot(Flux_l)
			// dot((double*)lxx_flux,ng,nl,(double*)Flux_l,nl,1,(double*)temp) ;
			// //py:k=-temp*c/x_xi x_xi=0.2
			// scalar_profuct((double*)temp,ng,1,-5.0,(double*)k);
			// //py:ug_new=k*dt+ug;
			// scalar_profuct((double*)k, ng,1,dt,(double*)temp2);
			// vec_sum((double*)temp2,(double*)ug,ng,1,(double*)ug_new) ;	
	
	 		// ug1=ug1_new;	
			// //dot_product((double*)A,m,n,(double*)B,x,y,(double*)C) 
			// dot((double*)g2l,nl,ng,(double*)ug1,ng,1,(double*)ul1) ;	
			// ul1[0]=ul[nl-1];
			// Flux_l=ul1;				
			// //py:temp=lxx_flux.dot(Flux_l)
			// dot((double*)lxx_flux,ng,nl,(double*)Flux_l,nl,1,(double*)temp) ;
			// //py:k=-temp*c/x_xi x_xi=0.2
			// scalar_profuct((double*)temp,ng,1,-5.0,(double*)k);
			// //py:ug_new=k*dt+ug;
			// scalar_profuct((double*)k, ng,1,dt,(double*)temp2);
			// vec_sum((double*)temp2,(double*)ug1,ng,1,(double*)ug1_new) ;	//


 	 		// ug2=ug2_new;	
			// //dot_product((double*)A,m,n,(double*)B,x,y,(double*)C) 
			// dot((double*)g2l,nl,ng,(double*)ug2,ng,1,(double*)ul2) ;	
			// ul2[0]=ul1[nl-1];
			// //ul[nl-1]=0;
			// Flux_l=ul2;	
			// //py:temp=lxx_flux.dot(Flux_l)
			// dot((double*)lxx_flux,ng,nl,(double*)Flux_l,nl,1,(double*)temp) ;
			// //py:k=-temp*c/x_xi x_xi=0.2
			// scalar_profuct((double*)temp,ng,1,-5.0,(double*)k);
			// //py:ug_new=k*dt+ug;
			// scalar_profuct((double*)k, ng,1,dt,(double*)temp2);
			// vec_sum((double*)temp2,(double*)ug2,ng,1,(double*)ug2_new) ;	//			

	 		// ug3=ug3_new;	
			// //dot_product((double*)A,m,n,(double*)B,x,y,(double*)C) 
			// dot((double*)g2l,nl,ng,(double*)ug3,ng,1,(double*)ul3) ;	
			// ul3[0]=ul2[nl-1];
			// //ul[nl-1]=0;
			// Flux_l=ul3;	
			// //py:temp=lxx_flux.dot(Flux_l)
			// dot((double*)lxx_flux,ng,nl,(double*)Flux_l,nl,1,(double*)temp) ;
			// //py:k=-temp*c/x_xi x_xi=0.2
			// scalar_profuct((double*)temp,ng,1,-5.0,(double*)k);
			// //py:ug_new=k*dt+ug;
			// scalar_profuct((double*)k, ng,1,dt,(double*)temp2);
			// vec_sum((double*)temp2,(double*)ug3,ng,1,(double*)ug3_new) ;	//

	 		// ug4=ug4_new;	
			// //dot_product((double*)A,m,n,(double*)B,x,y,(double*)C) 
			// dot((double*)g2l,nl,ng,(double*)ug4,ng,1,(double*)ul4) ;	
			// ul4[0]=ul3[nl-1];
			// //ul[nl-1]=0;
			// Flux_l=ul4;	
			// //py:temp=lxx_flux.dot(Flux_l)
			// dot((double*)lxx_flux,ng,nl,(double*)Flux_l,nl,1,(double*)temp) ;
			// //py:k=-temp*c/x_xi x_xi=0.2
			// scalar_profuct((double*)temp,ng,1,-5.0,(double*)k);
			// //py:ug_new=k*dt+ug;
			// scalar_profuct((double*)k, ng,1,dt,(double*)temp2);
			// vec_sum((double*)temp2,(double*)ug4,ng,1,(double*)ug4_new) ;	

												
			printf("\n" );
			for (int i=0;i<ng;i++){
				for (int j=0;j<1;j++){
					printf("%f,",ug[i]);
				}
			}	
		}//end one D wave	
		

		if (one_d_heat == 1){

			double* k  		  = new double[ng];
			double* temp   	  = new double[ng];
			double* Flux_lx   = new double[nl];
			double* temp2     = new double[ng];	
			double* temp3     = new double[ng];
			//double* temp333     = new double[nl];		


			// Get ul from ug
			//......................................

			for (int nnb=0;nnb<my_Ncv;nnb++){
        int nb;
        nb=my_chunk[nnb];
				//ug=cv[nb].ug_new_class;	
				get_value(ug,cv[nb].ug_new_class,ng,1);				
				//dot_product((double*)A,m,n,(double*)B,x,y,(double*)C) 
				dot(g2l,nl,ng,ug,ng,1,ul) ;	
				//cv[nb].ul_class=ul;
				get_value(cv[nb].ul_class,ul,nl,1);

			}	//end loop on block--->ug 2 ul 



			// Set BC on ul,get flux:"Flux_l" on ul
			// Get first flux:"Flux_l=ul" on l mesh
			// Get gradient of flux on g mesh
			// g-->l,get flux value "Flux_lx"on l mesh
			//......................................
			for (int nnb=0;nnb<my_Ncv;nnb++){
        int nb;
        nb=my_chunk[nnb];
				//ul=cv[nb].ul_class;
				get_value(ul,cv[nb].ul_class,nl,1);

				if (nb==0){
					ul[0]=(5.0+ul[0])/2.0;
					//ul[0] = (ul[0] + cv[Nblock-1].ul_class[nl-1])/2.0;
					ul[nl-1] = (cv[nb+1].ul_class[0]+ul[nl-1])/2.0;
				}

				if ((nb>0)&&(nb<(Nblock-1))){
					ul[0]    = (cv[nb-1].ul_class[nl-1]+ul[0])/2.0;
					ul[nl-1] = (cv[nb+1].ul_class[0]+ul[nl-1])/2.0;
				}
				
				if (nb==(Nblock-1)){
					ul[0]    = (cv[nb-1].ul_class[nl-1]+ul[0])/2.0;
					ul[nl-1]=(5.0+ul[nl-1])/2.0;
					//ul[nl-1] = (ul[nl-1] + cv[0].ul_class[0])/2.0;
				}
					
				//Flux_l=ul;	
				get_value(Flux_l,ul,nl,1);
				//py:temp=lxx_flux.dot(Flux_l) temp[ng,1]:first order derivative
				// temp2=temp/x_yita
				dot(lxx_flux,ng,nl,Flux_l,nl,1,temp) ;	
				scalar_profuct(temp,ng,1,1.0/x_yita,temp2);
				// Flux_lx[nl,1]:simply interpolate first derivative on l mesh				
				dot(g2l,nl,ng,temp2,ng,1,Flux_lx) ;	

				get_value(cv[nb].Flux_lx_class,Flux_lx,nl,1);

			}	//end loop on block--->get first derivative

	
			// Make sure the continuity of second flux on interface
			// Get gradient of second flux on g mesh
			// update ug_new=k*dt+ug_old
			//......................................
			for (int nnb=0;nnb<my_Ncv;nnb++){
        int nb;
        nb=my_chunk[nnb];
				get_value(Flux_lx,cv[nb].Flux_lx_class,nl,1);		
				get_value(ug,cv[nb].ug_new_class,ng,1);				

				//interface flux
				if (nb==0){
					Flux_lx[nl-1] = (cv[nb+1].Flux_lx_class[0]+Flux_lx[nl-1])/2.0;
				}

				if ((nb>0)&&(nb<(Nblock-1))){
					Flux_lx[0]    = (cv[nb-1].Flux_lx_class[nl-1]+Flux_lx[0])/2.0;
					Flux_lx[nl-1] = (cv[nb+1].Flux_lx_class[0]+Flux_lx[nl-1])/2.0;
				}
				
				if (nb==(Nblock-1)){
					Flux_lx[0]    = (cv[nb-1].Flux_lx_class[nl-1]+Flux_lx[0])/2.0;
				}


				//temp2[ng,1]:second order derivative
				dot(lxx_flux,ng,nl,Flux_lx,nl,1,temp2) ;
				//py:k=temp2/(x_yita)**2
				scalar_profuct(temp2,ng,1,1.0/x_yita,k);

				//py:ug_new=k*dt+ug;				
				scalar_profuct(k, ng,1,dt,temp3);
				vec_sum(temp3,ug,ng,1,temp2) ;

				get_value(cv[nb].ug_new_class,temp2,ng,1);

								
			}////end loop on block--->get second derivative

					printf("\n");	
					printf("~~~~~~~~\n");		

		}//end one D heat	
		if (two_d_heat==1){

			double* k_2d  		  = new double[ng*ng];

			double* Ttemp_ngnl     = new double[ng*nl];
			double* Ttemp_nlng     = new double[nl*ng];

			double* tempx1_ngng     = new double[ng*ng];
			double* tempx2_ngng     = new double[ng*ng];			
			double* tempy1_ngng     = new double[ng*ng];
			double* tempy2_ngng     = new double[ng*ng];			
			double* tempt1_ngng     = new double[ng*ng];
			double* tempt2_ngng     = new double[ng*ng];


			// initial metrix for mapping
			// Get ul from ug
			//......................................

			for (int nnb=0;nnb<my_Ncv;nnb++){
        int nb;
        nb=my_chunk[nnb];

			//......................................
				//ug=cv[nb].ug_new_class;	
				get_value(ug_2d,cv2d[nb].ug_new_class,ng,ng);				
				//ul_x=ug.dot(g2l.T)
				transposition(g2l,nl,ng,Ttemp_ngnl) ;
				dot(ug_2d,ng,ng,Ttemp_ngnl,ng,nl,cv2d[nb].ul_x_class) ;	
				//ul_y=(g_2.dot(ug))
				dot(g2l,nl,ng,ug_2d,ng,ng,cv2d[nb].ul_y_class) ;	
        prapere_test(nb,1);
        Sent(nb,cv2d[nb].ul_x_class,1);
        //Sent(nb,cv2d[nb].mpi_test,1);




			}	//end loop on block--->ug 2 ul 



			// Set BC on ul,get flux:"Flux_l" on ul
			// Get first flux:"Flux_l=ul" on l mesh
			// Get gradient of flux on g mesh
			// g-->l,get flux value "Flux_lx"on l mesh
			//......................................
			for (int nnb=0;nnb<my_Ncv;nnb++){
        int nb;
        nb=my_chunk[nnb];
				//ul_x_class= new double[ng*nl];
        Recv(nb,cv2d[nb].my_recv,1);
/*
        if ((nb==13)||(nb==46)){
          cout<<"my nb:"<<nb<<" my recv:"<<endl;
          for (int i=0;i<ng;i++) {
           cout<<cv2d[nb].my_recv[i]<<" ";
          }
          cout<<endl;
          cout<<endl;
        }

*/

				get_value(ul_2d_x,cv2d[nb].ul_x_class,ng,nl);
				get_value(ul_2d_y,cv2d[nb].ul_y_class,nl,ng);

				Initial_bc(nb);

				//Flux_l=ul;	
				get_value(Flux_l_x,ul_2d_x,ng,nl);
				get_value(Flux_l_y,ul_2d_y,nl,ng);
		
			/////////////////////////////////////////		
			/////////////////////////////////////////				
				//mapping(Flux_l_x, Flux_l_y, nb);

				Flux_first_map(Flux_l_x,Flux_l_y,tempx1_ngng,tempy1_ngng,nb);	

						if(n%oo_test==100000000000000){
						/////////////////////////////////////////
						//                TEST
						/////////////////////////////////////////
            if (my_rank==1){
						  if (nb==14){				
							  printf("\n");	
                cout<<"!!!!!!!!!   "<<nb<<endl;

							  for (int i=0;i<nl;i++){					
								  for (int j=0;j<ng;j++){
								  printf(",%f ",cv2d[nb].xX[i*ng+j]);		
								  }
								  printf("\n");
							  }
							  printf("\n");	
						  }

						  if (nb==28){				
							  printf("\n");	
                cout<<"!!!!!!!!!   "<<nb<<endl;

							  for (int i=0;i<nl;i++){					
								  for (int j=0;j<ng;j++){
								  printf(",%f ",cv2d[nb].xX[i*ng+j]);		
								  }
								  printf("\n");
							  }
							  printf("\n");	
						  }


						}
						/////////////////////////////////////////		
						/////////////////////////////////////////
				}	
				div_Jacobian(tempx1_ngng, tempy1_ngng, nb);
				//test_der(tempx1_ngng,tempy1_ngng);

				//scalar_profuct(temp,ng,1,1.0/x_yita,temp2);			
				//g-->l
				transposition(g2l,nl,ng,Ttemp_ngnl) ;
				dot(tempx1_ngng,ng,ng,Ttemp_ngnl,ng,nl,cv2d[nb].Flux_lx_class_F) ;	
				dot(tempy1_ngng,ng,ng,Ttemp_ngnl,ng,nl,cv2d[nb].Flux_lx_class_G) ;	

				dot(g2l,nl,ng,tempx1_ngng,ng,ng,cv2d[nb].Flux_ly_class_F) ;	
				dot(g2l,nl,ng,tempy1_ngng,ng,ng,cv2d[nb].Flux_ly_class_G) ;	

        Sent(nb,cv2d[nb].Flux_lx_class_F,2);
        Sent(nb,cv2d[nb].Flux_lx_class_G,3);







			}	//end loop on block--->get first derivative


			// Make sure the continuity of second flux on interface
			// Get gradient of second flux on g mesh
			// update ug_new=k*dt+ug_old
			//......................................
			for (int nnb=0;nnb<my_Ncv;nnb++){
        int nb;
        nb=my_chunk[nnb];
        //cout<<"my_rank:"<<my_rank<<"my cv:"<<nb<<" "<<my_Ncv<<endl;


				int l,r,u,d;
				int this_fa;
				int corr;
				//int l=nb-1,r=nb+1,u=nb-N_col,d=nb+N_col;
        Recv(nb,cv2d[nb].my_recv_F,2);
        Recv(nb,cv2d[nb].my_recv_G,3);

				get_value(Flux_lxF,cv2d[nb].Flux_lx_class_F,ng,nl);		//^_^
				get_value(Flux_lxG,cv2d[nb].Flux_lx_class_G,ng,nl);	
				get_value(Flux_lyF,cv2d[nb].Flux_ly_class_F,nl,ng);	
				get_value(Flux_lyG,cv2d[nb].Flux_ly_class_G,nl,ng);		//^_^

				get_value(ug_2d,cv2d[nb].ug_new_class,ng,ng);	

				Gradient_bc(nb);

				/////////////////////////////////////////		
				/////////////////////////////////////////	
				Flux_map(Flux_lxF,Flux_lxG,Flux_lyF,Flux_lyG, nb);

				transposition(lxx_flux,ng,nl,Ttemp_nlng) ;
				dot(Flux_lxF,ng,nl,Ttemp_nlng,nl,ng,tempx1_ngng) ;	
				dot(lyy_flux,ng,nl,Flux_lyG,nl,ng,tempy2_ngng) ;

				scalar_profuct(tempy2_ngng,ng,ng,-1,tempy1_ngng);	

				div_Jacobian(tempx1_ngng, tempy1_ngng, nb);

				//py:ug_new=k*dt+ug;	
				vec_sum(tempx1_ngng,tempy1_ngng,ng,ng,k_2d) ;			
				scalar_profuct(k_2d, ng,ng,dt,tempt1_ngng);
				vec_sum(tempt1_ngng,ug_2d,ng,ng,tempt2_ngng) ;

				get_value(cv2d[nb].ug_new_class,tempt2_ngng,ng,ng);

				if(n%oo_test==0){
						/////////////////////////////////////////
						//                TEST
						/////////////////////////////////////////
            if (my_rank==1){
						  if (nb==26){				
							  printf("\n");	
							  printf("iteration:%d\n",n);
							  printf("time:%d\n",t);
							  printf("aaCFD dt recomod%f\n",0.5*(xl[1]-xl[0])*(xl[1]-xl[0]) );
							  printf("dt=%f\n",dt );	

							  for (int i=0;i<ng;i++){					
								  for (int j=0;j<ng;j++){
								  printf(",%f ",tempt2_ngng[i*ng+j]);		
								  }
								  printf("\n");
							  }
							  printf("\n");	
						  }
						}
						/////////////////////////////////////////		
						/////////////////////////////////////////
				}	//end test

			
			}////end loop on block--->get second derivative


		}//end two D heat
	

			if (n==1){
        Sent_mesh();
    MPI_Barrier(MPI_COMM_WORLD);
        Recv_mesh();
        if (my_rank==0){
				out_put_mesh();	
        }	
				printf("mesh file saved!\n" );	
			}	





			if (n%oo_file==0){

        Sent_solution();
    MPI_Barrier(MPI_COMM_WORLD);
        Recv_solution();
        if (my_rank==0){
			    out_put_solution(n);	
        }

				printf("solution %d saved!\n",n );	
			
			}





		t=t+dt;

		}//time loop
    if (my_rank==0){
      struct timeval time2;
      gettimeofday(&time2,NULL);
      printf("wall clock used: %f seconds\n",
         (double)(time2.tv_sec-time1.tv_sec)+(double)(time2.tv_usec-time1.tv_usec)*1.0e-6);
    }

    ierr = MPI_Finalize();	
    return 0;  
	}


