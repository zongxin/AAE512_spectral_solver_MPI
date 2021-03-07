//#include<iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;


//void test(double *p,const int m,const int n,double *pshape,const int nx2,const int ny2)  

void dot(double *A,const int m,const int n,double *B,const int x,const int y,double *C)  {  

  if (n!=x){
    printf("n!=x" );

  }
  if ((n==0)||(m==0)||(x==0)||(y==0)){
    printf("dimention number should be zero!" );

  }

  if((m==1)&&(y==1)){
    C[0]=0;
    for(int i=0;i<n;i++)  {
      C[0]=C[0]+A[i]*B[i];
      }    
    }//varified

  else if((m==1)&&(y!=1)){

    for(int j=0;j<y;j++)  {
      double temp=0;
      for (int i=0;i<n;i++){
        temp=temp+A[i]*B[i*y+j];
      }
      C[j]=temp;
      }    
    }//varified

  else if((m!=1)&&(y==1)){

    for(int i=0;i<m;i++)  {
      double temp=0;
      for (int j=0;j<n;j++){
        temp=temp+A[i*n+j]*B[j];
      }
      C[i]=temp;
      }    
    }//varified

  else if((m!=1)&&(n==1)){
    for(int i=0;i<m;i++)  {
      for ( int k=0;k<y;k++){
        C[i*y+k]=A[i]*B[k];
       }
      }    
    } //varified

  else if((m!=1)&&(y!=1)){
    int k;
    double temp;
    for(int i=0;i<m;i++)  {


      for ( k=0;k<y;k++){
        temp=0.0;

        for (int j=0;j<n;j++){

        temp=temp+A[i*n+j]*B[j*y+k];

         }
 
        C[i*y+k]=temp;

       }

      }  

    } //verified!
}  

void scalar_profuct(double *A,const int m,const int n,const double k,double *C)  {  
  int D;
  if (m==1){ D=n;}
  if (n==1){ D=m;}

  if((m==1)||(n==1)){

    for(int i=0;i<D;i++)  {
      C[i]=A[i]*k;
      }    
    }//varified

  else {

    for(int i=0;i<m;i++)  {
      for (int j=0;j<n;j++){
        C[i*n+j]=A[i*n+j]*k;
        }
      }    
    }//varified

}  

void vec_sum(double *A,double *B,const int m,const int n,double *C)  {  
  int D;
  if (m==1){ D=n;}
  if (n==1){ D=m;}

  if((m==1)||(n==1)){

    for(int i=0;i<D;i++)  {
      C[i]=A[i]+B[i];
      }    

    }//varified

  else {

    for(int i=0;i<m;i++)  {
      for (int j=0;j<n;j++){
        C[i*n+j]=A[i*n+j]+B[i*n+j];
        }
      }    
    }//varified
}  

void vec_minus(double *A,double *B,const int m,const int n,double *C)  {  
  int D;
  if (m==1){ D=n;}
  if (n==1){ D=m;}

  if((m==1)||(n==1)){

    for(int i=0;i<D;i++)  {
      C[i]=A[i]-B[i];
      }    

    }//varified

  else {

    for(int i=0;i<m;i++)  {
      for (int j=0;j<n;j++){
        C[i*n+j]=A[i*n+j]-B[i*n+j];
        }
      }    
    }//varified
}  

void vec_product(double *A,double *B,const int m,const int n,double *C)  {  
  int D;
  if (m==1){ D=n;}
  if (n==1){ D=m;}

  if((m==1)||(n==1)){

    for(int i=0;i<D;i++)  {
      C[i]=A[i]*B[i];
      }    

    }//varified

  else {

    for(int i=0;i<m;i++)  {
      for (int j=0;j<n;j++){
        C[i*n+j]=A[i*n+j]*B[i*n+j];
        }
      }    
    }//varified
}  


void get_value(double *A,double *B,const int m,const int n)  {  
  int D;
  if (m==1){ D=n;}
  if (n==1){ D=m;}

  if((m==1)||(n==1)){

    for(int i=0;i<D;i++)  {
      A[i]=B[i];
      }    

    }//varified

  else {

    for(int i=0;i<m;i++)  {
      for (int j=0;j<n;j++){
        A[i*n+j]=B[i*n+j];
        }
      }    
    }//varified
}  

void transposition(double *A,const int m,const int n,double *T)  {  
    //A[m,n] B[n,m]
    for(int i=0;i<m;i++)  {
      for (int j=0;j<n;j++){
        T[j*m+i]=A[i*n+j];     
        }
      }    
}  

