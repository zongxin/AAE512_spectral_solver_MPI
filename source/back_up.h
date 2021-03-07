				if (Nblock>1){
				if (nb==0){
					for (int i=0;i<ng;i++){					
						ul_2d_x[nl*i+0]=(bc+ul_2d_x[nl*i+0])*0.5;	
						ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_x_class[nl*i+0]+ul_2d_x[nl*i+nl-1])/2.0;					
					}
					for (int j=0;j<ng;j++){					
						ul_2d_y[0*ng+j]=(bc+ul_2d_y[0*ng+j])*0.5;	
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_y_class[0*ng+j]+ul_2d_y[(nl-1)*ng+j])/2.0;						
					}

				}

				else if (nb==N_col-1){					
					for (int i=0;i<ng;i++){	
						ul_2d_x[nl*i+0]=(cv2d[l].ul_x_class[nl*i+nl-1]+ul_2d_x[nl*i+0])/2.0;				
						ul_2d_x[nl*i+nl-1]=(bc+ul_2d_x[nl*i+nl-1])*0.5;	
					}
					for (int j=0;j<ng;j++){					
						ul_2d_y[0*ng+j]=(bc+ul_2d_y[0*ng+j])*0.5;	
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_y_class[0*ng+j]+ul_2d_y[(nl-1)*ng+j])/2.0;
					}
				}
		
				


				else if (nb==(N_row-1)*N_col){
					for (int i=0;i<ng;i++){					
						ul_2d_x[nl*i+0]=(ul_2d_x[nl*i+0]+bc)*0.5;	
						ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_x_class[nl*i+0]+ul_2d_x[nl*i+nl-1])/2.0;	
					}
					for (int j=0;j<ng;j++){		
						ul_2d_y[0*ng+j]=(cv2d[u].ul_y_class[(nl-1)*ng+j]+ul_2d_y[0*ng+j])/2.0;				
						ul_2d_y[(nl-1)*ng+j]=(ul_2d_y[(nl-1)*ng+j]+bc)*0.5;	
					}
				}

				else if (nb==Nblock-1){

					for (int i=0;i<ng;i++){	
						ul_2d_x[nl*i+0]=(cv2d[l].ul_x_class[nl*i+nl-1]+ul_2d_x[nl*i+0])/2.0;				
						ul_2d_x[nl*i+nl-1]=(ul_2d_x[nl*i+nl-1]+bc)*0.5;	
					}
					for (int j=0;j<ng;j++){					
						ul_2d_y[0*ng+j]=(cv2d[u].ul_y_class[(nl-1)*ng+j]+ul_2d_y[0*ng+j])/2.0;
						ul_2d_y[(nl-1)*ng+j]=(ul_2d_y[(nl-1)*ng+j]+bc)*0.5;	
					}
				}


				else if (nb_row==0){

					for (int i=0;i<ng;i++){					
						ul_2d_x[nl*i+0]=(cv2d[l].ul_x_class[nl*i+nl-1]+ul_2d_x[nl*i+0])/2.0;	
						ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_x_class[nl*i+0]+ul_2d_x[nl*i+nl-1])/2.0;		
					}					
					for (int j=0;j<ng;j++){					
						ul_2d_y[0*ng+j]=(ul_2d_y[0*ng+j]+bc)*0.5;	
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_y_class[0*ng+j]+ul_2d_y[(nl-1)*ng+j])/2.0;	
					}				

				}
	
				else if (nb_row==N_row-1){
					for (int i=0;i<ng;i++){					
						ul_2d_x[nl*i+0]=(cv2d[l].ul_x_class[nl*i+nl-1]+ul_2d_x[nl*i+0])/2.0;	
						ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_x_class[nl*i+0]+ul_2d_x[nl*i+nl-1])/2.0;		
					}
					for (int j=0;j<ng;j++){					
						ul_2d_y[0*ng+j]=(cv2d[u].ul_y_class[(nl-1)*ng+j]+ul_2d_y[0*ng+j])/2.0;	
						ul_2d_y[(nl-1)*ng+j]=(ul_2d_y[(nl-1)*ng+j]+bc)*0.5;	
					}

				}

				else if (nb_col==0){

					for (int i=0;i<ng;i++){					
						ul_2d_x[nl*i+0]=(ul_2d_x[nl*i+0]+bc)*0.5;	
						ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_x_class[nl*i+0]+ul_2d_x[nl*i+nl-1])/2.0;	
					}
					for (int j=0;j<ng;j++){	
						ul_2d_y[0*ng+j]=(cv2d[u].ul_y_class[(nl-1)*ng+j]+ul_2d_y[0*ng+j])/2.0;	
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_y_class[0*ng+j]+ul_2d_y[(nl-1)*ng+j])/2.0;	
					}					
				}

				else if (nb_col==N_col-1){

					for (int i=0;i<ng;i++){		
						ul_2d_x[nl*i+0]=(cv2d[l].ul_x_class[nl*i+nl-1]+ul_2d_x[nl*i+0])/2.0;			
						ul_2d_x[nl*i+nl-1]=(ul_2d_x[nl*i+nl-1]+bc)*0.5;	
					}
					for (int j=0;j<ng;j++){	
						ul_2d_y[0*ng+j]=(cv2d[u].ul_y_class[(nl-1)*ng+j]+ul_2d_y[0*ng+j])/2.0;	
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_y_class[0*ng+j]+ul_2d_y[(nl-1)*ng+j])/2.0;	
					}					
				}


				else {
					for (int i=0;i<ng;i++){					
						ul_2d_x[nl*i+0]=(cv2d[l].ul_x_class[nl*i+nl-1]+ul_2d_x[nl*i+0])/2.0;	
						ul_2d_x[nl*i+nl-1]=(cv2d[r].ul_x_class[nl*i+0]+ul_2d_x[nl*i+nl-1])/2.0;		
					}
					for (int j=0;j<ng;j++){	
						ul_2d_y[0*ng+j]=(cv2d[u].ul_y_class[(nl-1)*ng+j]+ul_2d_y[0*ng+j])/2.0;	
						ul_2d_y[(nl-1)*ng+j]=(cv2d[d].ul_y_class[0*ng+j]+ul_2d_y[(nl-1)*ng+j])/2.0;	
					}
				}
				}//Nblock>1



				if (Nblock==1){
					for (int i=0;i<ng;i++){					
						ul_2d_x[nl*i+0]=(ul_2d_x[nl*i+0]+bc)*0.5;	
						ul_2d_x[nl*i+nl-1]=(ul_2d_x[nl*i+nl-1]+bc)*0.5;	
					}
					for (int j=0;j<ng;j++){		
						ul_2d_y[0*ng+j]=(ul_2d_y[0*ng+j]+bc)*0.5;					
						ul_2d_y[(nl-1)*ng+j]=(ul_2d_y[(nl-1)*ng+j]+bc)*0.5;	
					}				
				}















				//if (false){
				if (Nblock>1){			
				if (nb==0){
					for (int i=0;i<ng;i++){					
						//..Flux_lx[nl*i+0]=(bc+Flux_lx[nl*i+0])*0.5;	
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_lx_class_F[nl*i+0]+Flux_lxF[nl*i+nl-1])/2.0;	//				
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_lx_class_G[nl*i+0]+Flux_lxG[nl*i+nl-1])/2.0;	
					}
					for (int j=0;j<ng;j++){					
						//..Flux_l_y[0*ng+j]=(bc+Flux_l_y[0*ng+j])*0.5;	
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_ly_class_F[0*ng+j]+Flux_lyF[(nl-1)*ng+j])/2.0;	
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_ly_class_G[0*ng+j]+Flux_lyG[(nl-1)*ng+j])/2.0;	//				
					}

				}


				else if (nb==N_col-1){					
					for (int i=0;i<ng;i++){	
						Flux_lxF[nl*i+0]=(cv2d[l].Flux_lx_class_F[nl*i+nl-1]+Flux_lxF[nl*i+0])/2.0;		//		
						Flux_lxG[nl*i+0]=(cv2d[l].Flux_lx_class_G[nl*i+nl-1]+Flux_lxG[nl*i+0])/2.0;	
						//..Flux_lx[nl*i+nl-1]=(bc+Flux_lx[nl*i+nl-1])*0.5;	
					}
					for (int j=0;j<ng;j++){					
						//..Flux_l_y[0*ng+j]=(bc+Flux_l_y[0*ng+j])*0.5;	
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_ly_class_F[0*ng+j]+Flux_lyF[(nl-1)*ng+j])/2.0;
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_ly_class_G[0*ng+j]+Flux_lyG[(nl-1)*ng+j])/2.0;	//
					}
				}
		
				

				else if (nb==(N_row-1)*N_col){
					for (int i=0;i<ng;i++){					
						//..Flux_lx[nl*i+0]=(Flux_lx[nl*i+0]+bc)*0.5;	
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_lx_class_F[nl*i+0]+Flux_lxF[nl*i+nl-1])/2.0;	//
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_lx_class_G[nl*i+0]+Flux_lxG[nl*i+nl-1])/2.0;	
					}
					for (int j=0;j<ng;j++){		
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_ly_class_F[(nl-1)*ng+j]+Flux_lyF[0*ng+j])/2.0;	
						Flux_lyG[0*ng+j]=(cv2d[u].Flux_ly_class_G[(nl-1)*ng+j]+Flux_lyG[0*ng+j])/2.0;		//		
						//..Flux_l_y[(nl-1)*ng+j]=(Flux_l_y[(nl-1)*ng+j]+bc)*0.5;	
					}
				} 

				else if (nb==Nblock-1){

					for (int i=0;i<ng;i++){	
						Flux_lxF[nl*i+0]=(cv2d[l].Flux_lx_class_F[nl*i+nl-1]+Flux_lxF[nl*i+0])/2.0;		//
						Flux_lxG[nl*i+0]=(cv2d[l].Flux_lx_class_G[nl*i+nl-1]+Flux_lxG[nl*i+0])/2.0;			
						//..Flux_lx[nl*i+nl-1]=(Flux_lx[nl*i+nl-1]+bc)*0.5;	
					}
					for (int j=0;j<ng;j++){					
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_ly_class_F[(nl-1)*ng+j]+Flux_lyF[0*ng+j])/2.0;	
						Flux_lyG[0*ng+j]=(cv2d[u].Flux_ly_class_G[(nl-1)*ng+j]+Flux_lyG[0*ng+j])/2.0;	//
						//..Flux_l_y[(nl-1)*ng+j]=(Flux_l_y[(nl-1)*ng+j]+bc)*0.5;	
					}
				}

				else if (nb_row==0){

					for (int i=0;i<ng;i++){					
						Flux_lxF[nl*i+0]=(cv2d[l].Flux_lx_class_F[nl*i+nl-1]+Flux_lxF[nl*i+0])/2.0;		//
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_lx_class_F[nl*i+0]+Flux_lxF[nl*i+nl-1])/2.0;	//	

						Flux_lxG[nl*i+0]=(cv2d[l].Flux_lx_class_G[nl*i+nl-1]+Flux_lxG[nl*i+0])/2.0;	
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_lx_class_G[nl*i+0]+Flux_lxG[nl*i+nl-1])/2.0;		
					}					
					for (int j=0;j<ng;j++){					
						//..Flux_l_y[0*ng+j]=(Flux_l_y[0*ng+j]+bc)*0.5;	
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_ly_class_F[0*ng+j]+Flux_lyF[(nl-1)*ng+j])/2.0;							
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_ly_class_G[0*ng+j]+Flux_lyG[(nl-1)*ng+j])/2.0;	//

					}				

				}

				else if (nb_row==N_row-1){
					for (int i=0;i<ng;i++){					
						Flux_lxF[nl*i+0]=(cv2d[l].Flux_lx_class_F[nl*i+nl-1]+Flux_lxF[nl*i+0])/2.0;		//
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_lx_class_F[nl*i+0]+Flux_lxF[nl*i+nl-1])/2.0;	//	

						Flux_lxG[nl*i+0]=(cv2d[l].Flux_lx_class_G[nl*i+nl-1]+Flux_lxG[nl*i+0])/2.0;	
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_lx_class_G[nl*i+0]+Flux_lxG[nl*i+nl-1])/2.0;	
					}

					for (int j=0;j<ng;j++){					
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_ly_class_F[(nl-1)*ng+j]+Flux_lyF[0*ng+j])/2.0;	
						Flux_lyG[0*ng+j]=(cv2d[u].Flux_ly_class_G[(nl-1)*ng+j]+Flux_lyG[0*ng+j])/2.0;	//
						//..Flux_l_y[(nl-1)*ng+j]=(Flux_l_y[(nl-1)*ng+j]+bc)*0.5;	
					}

				}
	
				else if (nb_col==0){


					for (int i=0;i<ng;i++){					
						//..Flux_lx[nl*i+0]=(Flux_lx[nl*i+0]+bc)*0.5;	
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_lx_class_F[nl*i+0]+Flux_lxF[nl*i+nl-1])/2.0;	//
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_lx_class_G[nl*i+0]+Flux_lxG[nl*i+nl-1])/2.0;	
					}
					for (int j=0;j<ng;j++){	
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_ly_class_F[(nl-1)*ng+j]+Flux_lyF[0*ng+j])/2.0;	
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_F[0*ng+j]+Flux_lyF[(nl-1)*ng+j])/2.0;	

						Flux_lyG[0*ng+j]=(cv2d[u].Flux_ly_class_G[(nl-1)*ng+j]+Flux_lyG[0*ng+j])/2.0;		//
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_G[0*ng+j]+Flux_lyG[(nl-1)*ng+j])/2.0;	//
					}					
				}

				else if (nb_col==N_col-1){

					for (int i=0;i<ng;i++){		
						Flux_lxF[nl*i+0]=(cv2d[l].Flux_lx_class_F[nl*i+nl-1]+Flux_lxF[nl*i+0])/2.0;		//	
						Flux_lxG[nl*i+0]=(cv2d[l].Flux_lx_class_G[nl*i+nl-1]+Flux_lxG[nl*i+0])/2.0;								
						//..Flux_lx[nl*i+nl-1]=(Flux_lx[nl*i+nl-1]+bc)*0.5;	
					}
					for (int j=0;j<ng;j++){	
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_ly_class_F[(nl-1)*ng+j]+Flux_lyF[0*ng+j])/2.0;	
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_F[0*ng+j]+Flux_lyF[(nl-1)*ng+j])/2.0;	

						Flux_lyG[0*ng+j]=(cv2d[u].Flux_ly_class_G[(nl-1)*ng+j]+Flux_lyG[0*ng+j])/2.0;			//
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_G[0*ng+j]+Flux_lyG[(nl-1)*ng+j])/2.0;		//
					}					
				}

				else {
					for (int i=0;i<ng;i++){					
						Flux_lxF[nl*i+0]=(cv2d[l].Flux_lx_class_F[nl*i+nl-1]+Flux_lxF[nl*i+0])/2.0;			//
						Flux_lxF[nl*i+nl-1]=(cv2d[r].Flux_lx_class_F[nl*i+0]+Flux_lxF[nl*i+nl-1])/2.0;		//

						Flux_lxG[nl*i+0]=(cv2d[l].Flux_lx_class_G[nl*i+nl-1]+Flux_lxG[nl*i+0])/2.0;		
						Flux_lxG[nl*i+nl-1]=(cv2d[r].Flux_lx_class_G[nl*i+0]+Flux_lxG[nl*i+nl-1])/2.0;		
					}
					for (int j=0;j<ng;j++){	
						Flux_lyF[0*ng+j]=(cv2d[u].Flux_ly_class_F[(nl-1)*ng+j]+Flux_lyF[0*ng+j])/2.0;	
						Flux_lyF[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_F[0*ng+j]+Flux_lyF[(nl-1)*ng+j])/2.0;	

						Flux_lyG[0*ng+j]=(cv2d[u].Flux_ly_class_G[(nl-1)*ng+j]+Flux_lyG[0*ng+j])/2.0;		//
						Flux_lyG[(nl-1)*ng+j]=(cv2d[d].Flux_lx_class_G[0*ng+j]+Flux_lyG[(nl-1)*ng+j])/2.0;	//
					}
				}

				}//Nblock>1

				