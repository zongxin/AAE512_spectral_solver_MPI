#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;  



void out_put_solution(int n ){

	int one_d = 0;
	int two_d = 1;
	fstream fout;

    string num;
    stringstream ss;
    ss << n;
    ss >> num;//或者 res = ss.str();
    string filename="./temp/solution_";
    string tail=".dat";
    string name=filename+num+tail;
	fout.open(name.c_str(),ios::out);
	if (one_d==1){
		for (int nb=0;nb<Nblock;nb++){
			for (int i=0;i<ng;i++){
				//printf("%f \n", );
				fout << xg[i]+nb<<" "<<cv2d[nb].ug_new_class[i]<<"\n";
			}	
		}
	} //1d

	if (two_d==1){
		for (int nb=0;nb<Nblock;nb++){
			for (int i=0;i<ng;i++){
				for (int j=0;j<ng;j++){
					fout <<cv2d[nb].ug_new_class[i*ng+j]<<" ";
				}
			fout<<"\n";	
			}	
		}
	} //1d
	fout.close();   
}//end solution





void out_put_mesh(void ){
	fstream fout_x;
	fstream fout_y;

    string namex="./temp/mesh_x.dat";
    string namey="./temp/mesh_y.dat";

	fout_x.open(namex.c_str(),ios::out);
	for (int nb=0;nb<Nblock;nb++){		
		for (int i=0;i<ng;i++){
			for (int j=0;j<ng;j++){
				fout_x <<cv2d[nb].xxg[i*ng+j]<<" ";
			}
		fout_x<<"\n";	
		}	
	}
	fout_x.close();   

	fout_y.open(namey.c_str(),ios::out);	
	for (int nb=0;nb<Nblock;nb++){

		for (int i=0;i<ng;i++){
			for (int j=0;j<ng;j++){
				fout_y <<cv2d[nb].yyg[i*ng+j]<<" ";
			}
		fout_y<<"\n";	
		}	
	}
	fout_y.close();   
}//end mesh


void output_ini(void){

	int two_d = 1;
	fstream fout;

    string name="./temp/Initial.dat";
	fout.open(name.c_str(),ios::out);

	if (two_d==1){
		for (int nb=0;nb<Nblock;nb++){
			for (int i=0;i<ng;i++){
				for (int j=0;j<ng;j++){
					fout <<cv2d[nb].ug_new_class[i*ng+j]<<" ";
				}
			fout<<"\n";	
			}	
		}
	} //1d
	fout.close();   
}//end solution



void test_der(double* Fx,double* Fy){

	fstream fout;
    string name1="dx.dat";
    string name2="dy.dat" ; 
	fout.open(name1.c_str(),ios::out);
	for (int nb=0;nb<Nblock;nb++){
		for (int i=0;i<ng;i++){
			for (int j=0;j<ng;j++){
				fout <<Fx[i*ng+j]<<" ";
			}
		fout<<"\n";	
		}	
	}

	fout.close();   

	fout.open(name2.c_str(),ios::out);
	for (int nb=0;nb<Nblock;nb++){
		for (int i=0;i<ng;i++){
			for (int j=0;j<ng;j++){
				fout <<Fy[i*ng+j]<<" ";
			}
		fout<<"\n";	
		}	
	}

	fout.close();   
}//end solution


void test_l_mesh(double* Fx,double* Fy){

	fstream fout;
    string name1="dx.dat";
    string name2="dy.dat" ; 
	fout.open(name1.c_str(),ios::out);
	for (int nb=0;nb<Nblock;nb++){
		for (int i=0;i<ng;i++){
			for (int j=0;j<nl;j++){
				fout <<Fx[i*nl+j]<<" ";
			}
		fout<<"\n";	
		}	
	}

	fout.close();   

	fout.open(name2.c_str(),ios::out);
	for (int nb=0;nb<Nblock;nb++){
		for (int i=0;i<nl;i++){
			for (int j=0;j<ng;j++){
				fout <<Fy[i*ng+j]<<" ";
			}
		fout<<"\n";	
		}	
	}

	fout.close();   
}//end solution
