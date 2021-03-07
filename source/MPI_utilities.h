#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "transformation.h"
#include<mpi.h>



void Initial_MPI(void){

  if (my_rank==0){
    my_Ncv      =N_chunk1;
    my_chunk		 = new int[my_Ncv ];
	  for (int i=0;i<my_Ncv;i++){
      my_chunk[i]=chunk1[i];
    }
  }


  if (my_rank==1){
    my_Ncv      =N_chunk2;
    my_chunk		 = new int[my_Ncv ];
	  for (int i=0;i<my_Ncv;i++){
      my_chunk[i]=chunk2[i];
    }
  }

  if (my_rank==2){
    my_Ncv      =N_chunk3;
    my_chunk		 = new int[my_Ncv ];
	  for (int i=0;i<my_Ncv;i++){
      my_chunk[i]=chunk3[i];
    } 
  }
  if (my_rank==1){
  printf("my rank:%d\n",my_rank );
	for (int n=0;n<my_Ncv;n++){
  printf("%d ",my_chunk[n]);
  }
  printf("\n" );}
}

void Sent(int icv,double *A,int tag){
//A=cv2d[icv].ul_2d_x[nl*i+0]
  int r,l;
  r=cv2d[icv].my_nei[1];
  l=cv2d[icv].my_nei[3];

  if (my_rank==0){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk1_rcv[i]){
				for (int i=0;i<ng;i++){	
				  cv2d[icv].my_sent[i]=A[nl*i+nl-1];
        }
        MPI_Send(cv2d[icv].my_sent, ng, MPI_DOUBLE, 1, (1000)*tag+r, MPI_COMM_WORLD);
      }
    }



  }



  if (my_rank==1){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk2_rcv[i]){
				for (int i=0;i<ng;i++){	
				  cv2d[icv].my_sent[i]=A[nl*i+nl-1];
/*
                    if ((icv==26)||(icv==47)){
                      cout<<"my nb:"<<icv<<" my send:"<<endl;
                       cout<<cv2d[icv].my_sent[i]<<" ";
                      cout<<endl;
                    }
*/
        }
        MPI_Send(cv2d[icv].my_sent, ng, MPI_DOUBLE, 2, (1000)*tag+r, MPI_COMM_WORLD);
      }

      if (icv==chunk2_lcv[i]){
				for (int i=0;i<ng;i++){	
				  cv2d[icv].my_sent[i]=A[nl*i+0];
/*
                    if ((icv==26)||(icv==47)){
                      cout<<"my nb:"<<icv<<" my send:"<<endl;
                       cout<<cv2d[icv].my_sent[i]<<" ";
                      cout<<endl;
                    }

*/
        }
        MPI_Send(cv2d[icv].my_sent, ng, MPI_DOUBLE, 0, (1000)*tag+l, MPI_COMM_WORLD);
      }

    }

  }


  if (my_rank==2){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk3_lcv[i]){
				for (int i=0;i<ng;i++){	
				  cv2d[icv].my_sent[i]=A[nl*i+0];
        }
        MPI_Send(cv2d[icv].my_sent, ng, MPI_DOUBLE, 1, (1000)*tag+l, MPI_COMM_WORLD);
      }
    }
  }

}



void prapere_test(int icv,int tag){
//A=cv2d[icv].ul_2d_x[nl*i+0]
  int r,l;
  r=cv2d[icv].my_nei[1];
  l=cv2d[icv].my_nei[3];

  if (my_rank==0){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk1_rcv[i]){
	              for (int k=0;k<ng;k++){	
                	for (int j=0;j<nl;j++){	
                  cv2d[icv].mpi_test[nl*k+j]=r+0.1*k+0.01*j; 
                }
              }
      }
    }



  }



  if (my_rank==1){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk2_rcv[i]){
	              for (int k=0;k<ng;k++){	
                	for (int j=0;j<nl;j++){	
                  cv2d[icv].mpi_test[nl*k+j]=r+0.1*k+0.01*j; 
                }
              }
      }

      if (icv==chunk2_lcv[i]){
	              for (int k=0;k<ng;k++){	
                	for (int j=0;j<nl;j++){	
                  cv2d[icv].mpi_test[nl*k+j]=l+0.1*k+0.01*j; 
                }
              }
      }

    }

  }


  if (my_rank==2){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk3_lcv[i]){
	              for (int k=0;k<ng;k++){	
                	for (int j=0;j<nl;j++){	
                  cv2d[icv].mpi_test[nl*k+j]=l+0.1*k+0.01*j; 
                }
              }
      }
    }
  }

}



void Recv(int icv,double *b,int tag){
//b=cv2d[icv].my_recv
  int r,l;
  r=cv2d[icv].my_nei[1];
  l=cv2d[icv].my_nei[3];

  if (my_rank==0){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk1_rcv[i]){
        MPI_Recv(b, ng, MPI_DOUBLE, 1, 1000*tag+icv, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    }
  }


  if (my_rank==1){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk2_rcv[i]){
        MPI_Recv(b, ng, MPI_DOUBLE, 2, 1000*tag+icv, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }

      if (icv==chunk2_lcv[i]){

        MPI_Recv(b, ng, MPI_DOUBLE, 0, 1000*tag+icv, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }

    }

  }


  if (my_rank==2){

	  for (int i=0;i<N_mpifa;i++){
      if (icv==chunk3_lcv[i]){

        MPI_Recv(b, ng, MPI_DOUBLE, 1, 1000*tag+icv, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    }
  }

}


void Sent_solution(void){
  if (my_rank!=0){

			for (int nnb=0;nnb<my_Ncv;nnb++){
        int nb;
        nb=my_chunk[nnb];
        MPI_Send(cv2d[nb].ug_new_class, ng*ng, MPI_DOUBLE, 0, nb, MPI_COMM_WORLD);
      }
  }
}

void Recv_solution(void){
  if (my_rank==0){
			for (int nnb=0;nnb<N_chunk2;nnb++){
        int nb;
        nb=chunk2[nnb];
        MPI_Recv(cv2d[nb].ug_new_class, ng*ng, MPI_DOUBLE, 1, nb, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
  }

  if (my_rank==0){
			for (int nnb=0;nnb<N_chunk3;nnb++){
        int nb;
        nb=chunk3[nnb];
        MPI_Recv(cv2d[nb].ug_new_class, ng*ng, MPI_DOUBLE, 2, nb, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
  }
}



void Sent_mesh(void){
  if (my_rank!=0){

			for (int nnb=0;nnb<my_Ncv;nnb++){
        int nb;
        nb=my_chunk[nnb];
        MPI_Send(cv2d[nb].xxg, ng*ng, MPI_DOUBLE, 0, 5000+nb, MPI_COMM_WORLD);
        MPI_Send(cv2d[nb].yyg, ng*ng, MPI_DOUBLE, 0, 6000+nb, MPI_COMM_WORLD);
      }

  }
}


void Recv_mesh(void){
  if (my_rank==0){
			for (int nnb=0;nnb<N_chunk2;nnb++){
        int nb;
        nb=chunk2[nnb];
        MPI_Recv(cv2d[nb].xxg, ng*ng, MPI_DOUBLE, 1, 5000+nb, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(cv2d[nb].yyg, ng*ng, MPI_DOUBLE, 1, 6000+nb, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
  }

  if (my_rank==0){
			for (int nnb=0;nnb<N_chunk3;nnb++){
        int nb;
        nb=chunk3[nnb];
        MPI_Recv(cv2d[nb].xxg, ng*ng, MPI_DOUBLE, 2, 5000+nb, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(cv2d[nb].yyg, ng*ng, MPI_DOUBLE, 2, 6000+nb, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
  }
}



