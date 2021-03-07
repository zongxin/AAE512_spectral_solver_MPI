#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "variable.h"
using namespace std;  


int* noocv		 = new int[Ncv*4];
int* faocv		 = new int[Ncv*4];
int* neiocv		 = new int[Ncv*4];

double* c_list		 = new double[2*N_HOLE];
double* nox_ocv		 = new double[Ncv*4];
double* noy_ocv		 = new double[Ncv*4];
double* fax_ocv		 = new double[Ncv*4];
double* fay_ocv		 = new double[Ncv*4];

int* HOLE_cv		 = new int[NT_HOLEfa];
int* HOLE_fa		 = new int[NT_HOLEfa];

int* HOLE1_cv		 = new int[N_HOLEfa];
int* HOLE1_fa		 = new int[N_HOLEfa];
int* HOLE2_cv		 = new int[N_HOLEfa];
int* HOLE2_fa		 = new int[N_HOLEfa];
int* HOLE3_cv		 = new int[N_HOLEfa];
int* HOLE3_fa		 = new int[N_HOLEfa];
int* HOLE4_cv		 = new int[N_HOLEfa];
int* HOLE4_fa		 = new int[N_HOLEfa];
int* HOLE5_cv		 = new int[N_HOLEfa];
int* HOLE5_fa		 = new int[N_HOLEfa];

int* WALL_fa		 = new int[N_WALLfa];
int* WALL_cv		 = new int[N_WALLcv]; 

int* chunk1 		 = new int[N_chunk1];
int* chunk2 		 = new int[N_chunk2];
int* chunk3 		 = new int[N_chunk3];
int *fa1_list        = new int[N_mpifa] ;
int *fa2_list        = new int[N_mpifa] ;
int *chunk1_rcv		 = new int[N_mpifa] ;
int *chunk2_lcv		 = new int[N_mpifa] ;
int *chunk2_rcv		 = new int[N_mpifa] ;
int *chunk3_lcv		 = new int[N_mpifa] ;

void read_mesh_info(void){

	std::ifstream fin("./mesh_folder/boundary.txt", std::ios::in);
	char line[1024]={0};
	double x ;
	std::string name ;


	while(fin.getline(line, sizeof(line)))
	{

		std::stringstream word(line);


		word >> name;

		if(name=="HOLE_cv"){
			for(int i=0;i<NT_HOLEfa;i++){
				word>>HOLE_cv[i];
			}
		}

		if(name=="HOLE_fa"){
			for(int i=0;i<NT_HOLEfa;i++){
				word>>HOLE_fa[i];
			}
		}
		/*****************/

		if(name=="HOLE1_cv"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE1_cv[i];
			}
		}

		if(name=="HOLE1_fa"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE1_fa[i];
			}
		}
		/*****************/

		if(name=="HOLE2_cv"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE2_cv[i];
			}
		}
		if(name=="HOLE2_fa"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE2_fa[i];
			}
		}
		/*****************/
		if(name=="HOLE3_cv"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE3_cv[i];
			}
		}
		if(name=="HOLE3_fa"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE3_fa[i];
			}
		}
		/*****************/
		if(name=="HOLE4_cv"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE4_cv[i];
			}
		}
		if(name=="HOLE4_fa"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE4_fa[i];
			}
		}
		/*****************/
		if(name=="HOLE5_cv"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE5_cv[i];
			}
		}
		if(name=="HOLE5_fa"){
			for(int i=0;i<N_HOLEfa;i++){
				word>>HOLE5_fa[i];
			}
		}						
		//..........................................


		if(name=="WALL_cv"){
			for(int i=0;i<N_WALLcv;i++){
				word>>WALL_cv[i];
			}
		}


		if(name=="WALL_fa"){
			for(int i=0;i<N_WALLfa;i++){
				word>>WALL_fa[i];
			}
		}

						
		if(name=="C_list"){
			for(int i=0;i<2*N_HOLE;i++){
				word>>c_list[i];
			}
		}
		//..........................................
		//..............   MPI STUFF   .............
		//..........................................	

		if(name=="chunk1"){
			for(int i=0;i<N_chunk1;i++){
				word>>chunk1[i];
			}
		}			
		if(name=="chunk2"){
			for(int i=0;i<N_chunk2;i++){
				word>>chunk2[i];
			}
		}	

		if(name=="chunk3"){
			for(int i=0;i<N_chunk3;i++){
				word>>chunk3[i];
			}
		}	

		if(name=="fa1_list"){
			for(int i=0;i<N_mpifa;i++){
				word>>fa1_list[i];
			}
		}		

		if(name=="fa2_list"){
			for(int i=0;i<N_mpifa;i++){
				word>>fa2_list[i];
			}
		}	

		if(name=="chunk1_rcv"){
			for(int i=0;i<N_mpifa;i++){
				word>>chunk1_rcv[i];
			}
		}	

		if(name=="chunk2_lcv"){
			for(int i=0;i<N_mpifa;i++){
				word>>chunk2_lcv[i];
			}
		}

		if(name=="chunk2_rcv"){
			for(int i=0;i<N_mpifa;i++){
				word>>chunk2_rcv[i];
			}
		}				

		if(name=="chunk3_lcv"){
			for(int i=0;i<N_mpifa;i++){
				word>>chunk3_lcv[i];
			}
		}			
	}
					
	fin.clear();
	fin.close();

}




void read_faocv(void){

	std::ifstream fin("./mesh_folder/faocv.txt", std::ios::in);
	char line[1024]={0};
	int n_line=0;
	while(fin.getline(line, sizeof(line)))
	{	
		std::stringstream word(line);

		word >> faocv[n_line*4+0];
		word >> faocv[n_line*4+1];	
		word >> faocv[n_line*4+2];
		word >> faocv[n_line*4+3];	


		n_line++;	
	}

	fin.clear();
	fin.close();
}



void read_neiocv(void){

	std::ifstream fin("./mesh_folder/neiocv.txt", std::ios::in);
	char line[1024]={0};
	int n_line=0;

	while(fin.getline(line, sizeof(line)))
	{
		std::stringstream word(line);

		word >> neiocv[n_line*4+0];
		word >> neiocv[n_line*4+1];	
		word >> neiocv[n_line*4+2];
		word >> neiocv[n_line*4+3];	
		
		n_line++;	
	}

	fin.clear();
	fin.close();

}

void read_noocv(void){
	printf("reading boundary info ...\n");

	std::ifstream fin("./mesh_folder/noocv.txt", std::ios::in);
	char line[1024]={0};
	int n_line=0;

	while(fin.getline(line, sizeof(line)))
	{
		std::stringstream word(line);

		word >> noocv[n_line*4+0];
		word >> noocv[n_line*4+1];	
		word >> noocv[n_line*4+2];
		word >> noocv[n_line*4+3];	
		
		n_line++;	
	}

	fin.clear();
	fin.close();

}

void read_noy_ocv(void){

	std::ifstream fin("./mesh_folder/noy_ocv.txt", std::ios::in);
	char line[1024]={0};
	int n_line=0;

	while(fin.getline(line, sizeof(line)))
	{
		std::stringstream word(line);

		word >> noy_ocv[n_line*4+0];
		word >> noy_ocv[n_line*4+1];	
		word >> noy_ocv[n_line*4+2];
		word >> noy_ocv[n_line*4+3];
		
		n_line++;	
	}

	fin.clear();
	fin.close();

}

void read_nox_ocv(void){

	std::ifstream fin("./mesh_folder/nox_ocv.txt", std::ios::in);
	char line[1024]={0};
	int n_line=0;

	while(fin.getline(line, sizeof(line)))
	{
		std::stringstream word(line);

		word >> nox_ocv[n_line*4+0];
		word >> nox_ocv[n_line*4+1];	
		word >> nox_ocv[n_line*4+2];
		word >> nox_ocv[n_line*4+3];	
		n_line++;	
	}

	fin.clear();
	fin.close();


}

void read_mesh(void){

	read_mesh_info();

	read_noocv();
	read_faocv();
	read_neiocv();
	read_nox_ocv();
	read_noy_ocv();

}



/*
void read_test(void)
{
std::ifstream fin("./mesh_folder/faocv.txt", std::ios::in);
char line[1024]={0};
std::string x = "";
std::string y = "";
std::string z = "";
while(fin.getline(line, sizeof(line)))
{
std::stringstream word(line);
word >> x;
word >> y;
word >> z;
std::cout << "x: " << x << std::endl;
std::cout << "y: " << y << std::endl;
std::cout << "z: " << z << std::endl;
}
fin.clear();
fin.close();
}
*/
