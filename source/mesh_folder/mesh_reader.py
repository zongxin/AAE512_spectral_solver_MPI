import os
import sys
import numpy as np
import scipy as sp 
import pylab as plt
from pdb import set_trace
from matplotlib import rc as matplotlibrc
import umesh_reader,copy
#import math_tool as ma
#import uns_plot as uns
import matplotlib as mpl
from matplotlib import cm

#the following enables LaTeX typesetting, which will cause the plotting to take forever..  
matplotlibrc('text.latex', preamble='\usepackage{color}')
matplotlibrc('text',usetex=True)
matplotlibrc('font', family='serif')

Problem1 = True
HAVE_HOLE= True
MUL 	 = True

LR_FLIP=0
UD_FLIP=0





fig_width = 10
fig_height = 10
textFontSize   = 10
gcafontSize    = 26
lineWidth      = 2  
icemcfd_project_folder = 'mesh_folder/'

####################################################################################
####################################  Problem1  ####################################
####################################################################################
if Problem1:
	#filename = "nine.msh"
	#filename = 'single11.msh'
	#filename = 'fluent.msh'
	#filename = 'one.msh'	
	#filename = 'two.msh'
	#filename = 'one_skew.msh'		
	#filename = 'quad.msh'
	#filename =	'7_shape.msh'
	#filename = 'upper.msh'	
	#filename = 'eye.msh'
	#filename  = 'final.msh'	
	filename  = 'final2.msh'	
	c_list=np.zeros((5,2))
	c_list[0,:]=np.array([10,6 ])
	c_list[1,:]=np.array([6 ,10])	
	c_list[2,:]=np.array([2 ,6])	
	c_list[3,:]=np.array([6 ,2])	
	c_list[4,:]=np.array([6 ,6])		

	Phi_dict={}
	kk=0
	if kk ==0:
		mshfile_fullpath = filename
		part_names, xy_no, xy_fa, xy_cv, noofa, cvofa, faono, faocv, partofa = \
						umesh_reader.read_unstructured_grid(mshfile_fullpath,node_reordering=True)
		if UD_FLIP==1:
			filename = 'down.msh'
			xy_no[:,1]=- xy_no[:,1]
			xy_fa[:,1]=- xy_fa[:,1]
			xy_cv[:,1]=- xy_cv[:,1]						
		if True:
			############################################################
			########################    Group   ########################
			############################################################
			B_fa = {}			

			i = 0
			nno = xy_no.shape[0]
			ncv = xy_cv.shape[0]
			nfa = xy_fa.shape[0]
			phi  = np.zeros(nno)			
	
			for ifa, part in enumerate(partofa):
			  if part != 'SOLID':
			    B_fa[i] = ifa
			    i = i+1

			B_fa=list(B_fa.values())			
			############################################################
			######################    face list  #######################
			############################################################

			WALL_fa={}
			for ifa, part in enumerate(partofa):
			  if part == 'WALL':
			    WALL_fa[i] = ifa
			    i = i+1
			WALL_fa=np.array(list(WALL_fa.values())	)
			WALL_nfa=len(WALL_fa)



			HOLE1_fa={}
			for ifa, part in enumerate(partofa):
			  if part == 'HOLE1':
			    HOLE1_fa[i] = ifa
			    i = i+1
			HOLE1_fa=np.array(list(HOLE1_fa.values())	)
			HOLE1_nfa=len(HOLE1_fa)

			HOLE2_fa={}
			for ifa, part in enumerate(partofa):
			  if part == 'HOLE2':
			    HOLE2_fa[i] = ifa
			    i = i+1
			HOLE2_fa=np.array(list(HOLE2_fa.values())	)
			HOLE2_nfa=len(HOLE2_fa)

			HOLE3_fa={}
			for ifa, part in enumerate(partofa):
			  if part == 'HOLE3':
			    HOLE3_fa[i] = ifa
			    i = i+1
			HOLE3_fa=np.array(list(HOLE3_fa.values())	)
			HOLE3_nfa=len(HOLE3_fa)

			HOLE4_fa={}
			for ifa, part in enumerate(partofa):
			  if part == 'HOLE4':
			    HOLE4_fa[i] = ifa
			    i = i+1
			HOLE4_fa=np.array(list(HOLE4_fa.values())	)
			HOLE4_nfa=len(HOLE4_fa)

			HOLE5_fa={}
			for ifa, part in enumerate(partofa):
			  if part == 'HOLE5':
			    HOLE5_fa[i] = ifa
			    i = i+1
			HOLE5_fa=np.array(list(HOLE5_fa.values())	)
			HOLE5_nfa=len(HOLE5_fa)

			HOLE_fa={}
			for ifa, part in enumerate(partofa):
			  if part == 'HOLE':
			    HOLE_fa[i] = ifa
			    i = i+1
			HOLE_fa=np.array(list(HOLE_fa.values())	)
			HOLE_nfa=len(HOLE_fa)

			if MUL:
				HOLE_fa=np.zeros((5,HOLE1_nfa))
				HOLE_fa[0,:]=HOLE1_fa
				HOLE_fa[1,:]=HOLE2_fa
				HOLE_fa[2,:]=HOLE3_fa
				HOLE_fa[3,:]=HOLE4_fa
				HOLE_fa[4,:]=HOLE5_fa												
				HOLE_fa=HOLE_fa.flatten()
				HOLE_nfa=len(HOLE_fa)
			#set_trace()	


			############################################################
			###################    hole face sort  #####################
			############################################################
			if (MUL==False):
				c=np.zeros(2)
				x=np.zeros(HOLE_nfa)
				y=np.zeros(HOLE_nfa)
				for i,fa in enumerate(HOLE_fa)	:	

					x[i]=xy_fa[fa,0]
					y[i]=xy_fa[fa,1]

				# create a new list of corners which includes angles
				cornersWithAngles = []
				ii=-1

				for xx, yy in zip(x,y):
					ii=ii+1
					dx = xx - c[0]
					dy = yy - c[1]
				 	an = np.arctan(dy/dx)
					if (dx < 0):
						an = an +np.pi 
					if ((dx<0)&(dy<0)):
						an=an+2*np.pi
					
					cornersWithAngles.append((HOLE_fa[ii], an))
				# sort it using the angles
				cornersWithAngles.sort(key = lambda tup: tup[1])
				temp = np.array(cornersWithAngles)
				for i,fa in enumerate(HOLE_fa)	:
					HOLE_fa[i]=temp[i][0]


			#................................................
			c=c_list[0,:]
			x=np.zeros(HOLE1_nfa)
			y=np.zeros(HOLE1_nfa)
			for i,fa in enumerate(HOLE1_fa)	:	

				x[i]=xy_fa[fa,0]
				y[i]=xy_fa[fa,1]

			# create a new list of corners which includes angles
			cornersWithAngles = []
			ii=-1

			for xx, yy in zip(x,y):
				ii=ii+1
				dx = xx - c[0]
				dy = yy - c[1]
			 	an = np.arctan(dy/dx)
				if (dx < 0):
					an = an +np.pi 
				if ((dx<0)&(dy<0)):
					an=an+2*np.pi
				
				cornersWithAngles.append((HOLE1_fa[ii], an))
			# sort it using the angles
			cornersWithAngles.sort(key = lambda tup: tup[1])
			temp = np.array(cornersWithAngles)
			for i,fa in enumerate(HOLE1_fa)	:
				HOLE1_fa[i]=temp[i][0]


			#................................................
			c=c_list[1,:]
			x=np.zeros(HOLE2_nfa)
			y=np.zeros(HOLE2_nfa)
			for i,fa in enumerate(HOLE2_fa)	:	

				x[i]=xy_fa[fa,0]
				y[i]=xy_fa[fa,1]

			# create a new list of corners which includes angles
			cornersWithAngles = []
			ii=-1

			for xx, yy in zip(x,y):
				ii=ii+1
				dx = xx - c[0]
				dy = yy - c[1]
			 	an = np.arctan(dy/dx)
				if (dx < 0):
					an = an +np.pi 
				if ((dx<0)&(dy<0)):
					an=an+2*np.pi
				
				cornersWithAngles.append((HOLE2_fa[ii], an))
			# sort it using the angles
			cornersWithAngles.sort(key = lambda tup: tup[1])
			temp = np.array(cornersWithAngles)
			for i,fa in enumerate(HOLE2_fa)	:
				HOLE2_fa[i]=temp[i][0]

			#................................................
			c=c_list[2,:]
			x=np.zeros(HOLE3_nfa)
			y=np.zeros(HOLE3_nfa)
			for i,fa in enumerate(HOLE3_fa)	:	

				x[i]=xy_fa[fa,0]
				y[i]=xy_fa[fa,1]

			# create a new list of corners which includes angles
			cornersWithAngles = []
			ii=-1

			for xx, yy in zip(x,y):
				ii=ii+1
				dx = xx - c[0]
				dy = yy - c[1]
			 	an = np.arctan(dy/dx)
				if (dx < 0):
					an = an +np.pi 
				if ((dx<0)&(dy<0)):
					an=an+2*np.pi
				
				cornersWithAngles.append((HOLE3_fa[ii], an))
			# sort it using the angles
			cornersWithAngles.sort(key = lambda tup: tup[1])
			temp = np.array(cornersWithAngles)
			for i,fa in enumerate(HOLE3_fa)	:
				HOLE3_fa[i]=temp[i][0]

			#................................................
			c=c_list[3,:]
			x=np.zeros(HOLE4_nfa)
			y=np.zeros(HOLE4_nfa)
			for i,fa in enumerate(HOLE4_fa)	:	

				x[i]=xy_fa[fa,0]
				y[i]=xy_fa[fa,1]

			# create a new list of corners which includes angles
			cornersWithAngles = []
			ii=-1

			for xx, yy in zip(x,y):
				ii=ii+1
				dx = xx - c[0]
				dy = yy - c[1]
			 	an = np.arctan(dy/dx)
				if (dx < 0):
					an = an +np.pi 
				if ((dx<0)&(dy<0)):
					an=an+2*np.pi
				
				cornersWithAngles.append((HOLE4_fa[ii], an))
			# sort it using the angles
			cornersWithAngles.sort(key = lambda tup: tup[1])
			temp = np.array(cornersWithAngles)
			for i,fa in enumerate(HOLE4_fa)	:
				HOLE4_fa[i]=temp[i][0]
			#................................................
			c=c_list[4,:]
			x=np.zeros(HOLE5_nfa)
			y=np.zeros(HOLE5_nfa)
			for i,fa in enumerate(HOLE5_fa)	:	

				x[i]=xy_fa[fa,0]
				y[i]=xy_fa[fa,1]

			# create a new list of corners which includes angles
			cornersWithAngles = []
			ii=-1

			for xx, yy in zip(x,y):
				ii=ii+1
				dx = xx - c[0]
				dy = yy - c[1]
			 	an = np.arctan(dy/dx)
				if (dx < 0):
					an = an +np.pi 
				if ((dx<0)&(dy<0)):
					an=an+2*np.pi
				
				cornersWithAngles.append((HOLE5_fa[ii], an))
			# sort it using the angles
			cornersWithAngles.sort(key = lambda tup: tup[1])
			temp = np.array(cornersWithAngles)
			for i,fa in enumerate(HOLE5_fa)	:
				HOLE5_fa[i]=temp[i][0]

			############################################################
			###################    WALL CV.        #####################
			############################################################
			B_cv  = {}
			B_no  = {}
			B_noo = {}
			B_cv0 = {}
			B_cvv = {}
			for i, ifa in enumerate(B_fa):
				B_noo[i] = noofa[ifa]
				B_no = list(set(B_no).union(set(B_noo[i])))
				B_cvv[i] = cvofa[ifa]
				B_cv0 = list(set(B_cv0).union(set(B_cvv[i])))			

			B_cv = list(B_cv0) # makes a copy
			# B_cvnew= B_cv # does NOT make a copy, it create a new variable that
			#                 points to the same list.. 
			
			B_cv.remove(-1)			
			B_cv=np.array(B_cv)


			#................................................
			WALL_cv  = {}
			WALL_no  = {}
			WALL_noo = {}
			WALL_cv0 = {}
			WALL_cvv = {}
			for i, ifa in enumerate(WALL_fa):
				WALL_noo[i] = noofa[ifa]
				WALL_no = list(set(WALL_no).union(set(WALL_noo[i])))
				WALL_cvv[i] = cvofa[ifa]
				WALL_cv0 = list(set(WALL_cv0).union(set(WALL_cvv[i])))			

			WALL_cv = list(WALL_cv0) # makes a copy
			# WALL_cvnew= WALL_cv # does NOT make a copy, it create a new variaWALLle that
			#                 points to the same list.. 
			
			WALL_cv.remove(-1)			
			WALL_cv=np.array(WALL_cv)


			WALL_ncv	=len(WALL_cv)


			# ................................................
			# HOLE_cv  = {}
			# HOLE_no  = {}
			# HOLE_noo = {}
			# HOLE_cv0 = {}
			# HOLE_cvv = {}
			# for i, ifa in enumerate(HOLE_fa):
			# 	HOLE_noo[i] = noofa[ifa]
			# 	HOLE_no = list(set(HOLE_no).union(set(HOLE_noo[i])))
			# 	HOLE_cvv[i] = cvofa[ifa]
			# 	HOLE_cv0 = list(set(HOLE_cv0).union(set(HOLE_cvv[i])))			

			# HOLE_cv = list(HOLE_cv0) # makes a copy
			# # HOLE_cvnew= HOLE_cv # does NOT make a copy, it create a new variaHOLEle that
			# #                 points to the same list.. 
			# if (HAVE_HOLE):
			# 	HOLE_cv.remove(-1)	
			# HOLE_cv=np.array(HOLE_cv)		
			# HOLE_ncv=len(HOLE_cv)


			HOLE1_cv=np.zeros(HOLE1_nfa)
			for i,fa in enumerate(HOLE1_fa):
				HOLE1_cv[i]=np.int(cvofa[fa,0]+cvofa[fa,1]+1)

			HOLE2_cv=np.zeros(HOLE2_nfa)
			for i,fa in enumerate(HOLE2_fa):
				HOLE2_cv[i]=np.int(cvofa[fa,0]+cvofa[fa,1]+1)	

			HOLE3_cv=np.zeros(HOLE3_nfa)
			for i,fa in enumerate(HOLE3_fa):
				HOLE3_cv[i]=np.int(cvofa[fa,0]+cvofa[fa,1]+1)

			HOLE4_cv=np.zeros(HOLE4_nfa)
			for i,fa in enumerate(HOLE4_fa):
				HOLE4_cv[i]=np.int(cvofa[fa,0]+cvofa[fa,1]+1)

			HOLE5_cv=np.zeros(HOLE5_nfa)
			for i,fa in enumerate(HOLE5_fa):
				HOLE5_cv[i]=np.int(cvofa[fa,0]+cvofa[fa,1]+1)



			if MUL:
				HOLE_cv=np.zeros((5,HOLE1_nfa))
				HOLE_cv[0,:]=HOLE1_cv
				HOLE_cv[1,:]=HOLE2_cv
				HOLE_cv[2,:]=HOLE3_cv
				HOLE_cv[3,:]=HOLE4_cv
				HOLE_cv[4,:]=HOLE5_cv												
				HOLE_cv=HOLE_cv.flatten()
				HOLE_ncv=len(HOLE_cv)		
			else:		
				HOLE_cv=np.zeros(HOLE_nfa)
				for i,fa in enumerate(HOLE_fa):
					HOLE_cv[i]=np.int(cvofa[fa,0]+cvofa[fa,1]+1)										
			#................................................

	

			#Obtain noocv 
			noocv = {}
			for icv, face in enumerate(faocv):
				temp  = list(face)
				temp2 = {}
				for i, iface in enumerate(temp):
					temp2 = list(set(temp2).union(set(noofa[iface])))
				noocv[icv] = temp2
			noocv=np.array(list(noocv.values()))



			for i in range(ncv)	:	
				nolist=noocv[i,:]	
				xy=np.zeros((4,2))
				for ino,no in enumerate(nolist):
					xy[ino,:]=xy_no[no,:]
				x=xy[:,0]
				y=xy[:,1]

				cx = sum(x)/len(x)
				cy = sum(y)/len(y)
				# create a new list of corners which includes angles
				cornersWithAngles = []
				ii=-1

				for xx, yy in zip(x,y):
					ii=ii+1
					dx = xx - cx
					dy = yy - cy
				 	an = np.arctan(dy/dx)
					if dx < 0:
						an = an +np.pi 
					cornersWithAngles.append((nolist[ii], an))
				# sort it using the angles
				cornersWithAngles.sort(key = lambda tup: tup[1])
				temp = np.array(cornersWithAngles)

				noocv[i,0]=temp[3][0]
				noocv[i,1]=temp[0][0]
				noocv[i,2]=temp[1][0]
				noocv[i,3]=temp[2][0]				









			faocv=np.array(faocv)

			for i in range(ncv)	:
				falist=faocv[i,:]	
				temp=np.zeros(4)
				for ifa,fa in enumerate(falist):
					no1=noofa[fa][0]
					no2=noofa[fa][1]

					if (((noocv[i,0]==no1)&(noocv[i,1]==no2))|((noocv[i,0]==no2)&(noocv[i,1]==no1))):
						temp[0]=fa

					elif (((noocv[i,1]==no1)&(noocv[i,2]==no2))|((noocv[i,1]==no2)&(noocv[i,2]==no1))):
						temp[1]=fa

					elif (((noocv[i,2]==no1)&(noocv[i,3]==no2))|((noocv[i,2]==no2)&(noocv[i,3]==no1))):
						temp[2]=fa

					elif (((noocv[i,3]==no1)&(noocv[i,0]==no2))|((noocv[i,3]==no2)&(noocv[i,0]==no1))):
						temp[3]=fa

					else:
						print "wrong in sort order for faocv"
				faocv[i,:]=temp



			neiocv=np.zeros((ncv,4))
			for i in range(ncv)	:
				falist=faocv[i,:]	

				for ifa,fa in enumerate(falist):
					temp=cvofa[fa,:]
					if temp[0]!=i:
						neiocv[i,ifa]=temp[0]
					else:
						neiocv[i,ifa]=temp[1]







			nox_ocv=np.zeros(noocv.shape)
			noy_ocv=np.zeros(noocv.shape)

			for i in range(ncv):
				nolist=noocv[i,:]
				for j in range(4):
					nox_ocv[i,j]=xy_no[nolist[j]][0]
					noy_ocv[i,j]=xy_no[nolist[j]][1]

			fax_ocv=np.zeros(noocv.shape)
			fay_ocv=np.zeros(noocv.shape)

			for i in range(ncv):
				nolist=noocv[i,:]
				for j in range(4):
					fax_ocv[i,j]=xy_fa[nolist[j]][0]
					fay_ocv[i,j]=xy_fa[nolist[j]][1]


			#Obtain cv of inner notes: cvono_in
			In_no = range(nno)
			for i, b_no in enumerate(B_no):
				In_no.remove(b_no)			

			In_fa = range(nfa)
			for i, b_fa in enumerate(B_fa):
				In_fa.remove(b_fa)			





			MPI_size=3
			fa1=(np.max(xy_no[:,0])-np.min(xy_no[:,0]))/3.0+np.min(xy_no[:,0])
			fa2=fa1*2.0

			chunk1 = {}
			chunk2 = {}
			chunk3 = {}
			for i in range(ncv):
				if xy_cv[i,0]<fa1:
					chunk1[str(i)]=i
				elif (xy_cv[i,0]>fa1)&(xy_cv[i,0]<fa2):
					chunk2[str(i)]=i
				else:					
					chunk3[str(i)]=i
			chunk1=chunk1.values()
			chunk2=chunk2.values()
			chunk3=chunk3.values()		

			fa1_list={}
			fa2_list={}			
			for i in range(nfa):
				if (np.abs(xy_fa[i,0]-fa1)< 0.00001):
					fa1_list[str(i)]=i
				elif (np.abs(xy_fa[i,0]-fa2)< 0.00001):
					fa2_list[str(i)]=i;
			fa1_list=fa1_list.values()
			fa2_list=fa2_list.values()

			chunk1_rcv=np.zeros(len(fa1_list))
			chunk2_lcv=np.zeros(len(fa1_list))			
			chunk2_rcv=np.zeros(len(fa2_list))
			chunk3_lcv=np.zeros(len(fa2_list))

			for i,fa in enumerate(fa1_list):
				if (cvofa[fa][0] in chunk1)&(cvofa[fa][1] in chunk2):
					chunk1_rcv[i]=cvofa[fa][0]
					chunk2_lcv[i]=cvofa[fa][1]
				elif (cvofa[fa][0] in chunk2)&(cvofa[fa][1] in chunk1):
					chunk1_rcv[i]=cvofa[fa][1]
					chunk2_lcv[i]=cvofa[fa][0]

			for i,fa in enumerate(fa2_list):
				if (cvofa[fa][0] in chunk2)&(cvofa[fa][1] in chunk3):
					chunk2_rcv[i]=cvofa[fa][0]
					chunk3_lcv[i]=cvofa[fa][1]
				elif (cvofa[fa][0] in chunk3)&(cvofa[fa][1] in chunk2):
					chunk2_rcv[i]=cvofa[fa][1]
					chunk3_lcv[i]=cvofa[fa][0]
			set_trace()
#			xy_no=xy_no*1.0/3
#			xy_fa=xy_fa*1.0/3
#			nox_ocv=nox_ocv*1.0/3
#			noy_ocv=noy_ocv*1.0/3			
# 			fax_ocv=nox_ocv*1.0/3
#			fay_ocv=noy_ocv*1.0/3	
			############################################################
			########################    Out put  #######################
			############################################################			
			f=open("xyono.txt",'w')
			for i in range(nno):
				f.write(str(xy_no[i,0])+" "+str(xy_no[i,1])+"\n")
			f.close()	

			f=open("xyofa.txt",'w')
			for i in range(nno):
				f.write(str(xy_fa[i,0])+" "+str(xy_fa[i,1])+"\n")
			f.close()	

			f=open("faocv.txt",'w')
			for i in range(ncv):
				f.write(str(faocv[i][0])+" "+str(faocv[i][1])+" "+str(faocv[i][2])+" "+str(faocv[i][3])+"\n")
			f.close()	

			f=open("noocv.txt",'w')
			for i in range(ncv):
				f.write(str(noocv[i][0])+" "+str(noocv[i][1])+" "+str(noocv[i][2])+" "+str(noocv[i][3])+"\n")
			f.close()	

			f=open("neiocv.txt",'w')
			for i in range(ncv):
				f.write(str(np.int(neiocv[i][0]))+" "+str(np.int(neiocv[i][1]))+" "+str(np.int(neiocv[i][2]))+" "+str(np.int(neiocv[i][3]))+"\n")
			f.close()	

			f=open("nox_ocv.txt",'w')
			for i in range(ncv):
				f.write(str(nox_ocv[i][0])+" "+str(nox_ocv[i][1])+" "+str(nox_ocv[i][2])+" "+str(nox_ocv[i][3])+"\n")
			f.close()	

			f=open("noy_ocv.txt",'w')
			for i in range(ncv):
				f.write(str(noy_ocv[i][0])+" "+str(noy_ocv[i][1])+" "+str(noy_ocv[i][2])+" "+str(noy_ocv[i][3])+"\n")
			f.close()	

			f=open("fax_ocv.txt",'w')
			for i in range(ncv):
				f.write(str(fax_ocv[i][0])+" "+str(fax_ocv[i][1])+" "+str(fax_ocv[i][2])+" "+str(fax_ocv[i][3])+"\n")
			f.close()	

			f=open("fay_ocv.txt",'w')
			for i in range(ncv):
				f.write(str(fay_ocv[i][0])+" "+str(fay_ocv[i][1])+" "+str(fay_ocv[i][2])+" "+str(fay_ocv[i][3])+"\n")
			f.close()



			f=open("boundary.txt",'w')
			###########################
			f.write("HOLE_cv")
			if (HAVE_HOLE):			
				for i in range(HOLE_ncv):
					f.write(" "+str(np.int(HOLE_cv[i])))
			else:
				f.write(" "+str(np.int(-1)))	
			f.write("\n")

			f.write("HOLE1_cv")
			if (HAVE_HOLE):			
				for i in range(HOLE1_nfa):
					f.write(" "+str(np.int(HOLE1_cv[i])))
			else:
				f.write(" "+str(np.int(-1)))	
			f.write("\n")

			f.write("HOLE3_cv")
			if (HAVE_HOLE):			
				for i in range(HOLE3_nfa):
					f.write(" "+str(np.int(HOLE3_cv[i])))
			else:
				f.write(" "+str(np.int(-1)))	
			f.write("\n")		

			f.write("HOLE4_cv")
			if (HAVE_HOLE):			
				for i in range(HOLE4_nfa):
					f.write(" "+str(np.int(HOLE4_cv[i])))
			else:
				f.write(" "+str(np.int(-1)))	
			f.write("\n")		

			f.write("HOLE5_cv")
			if (HAVE_HOLE):			
				for i in range(HOLE5_nfa):
					f.write(" "+str(np.int(HOLE5_cv[i])))
			else:
				f.write(" "+str(np.int(-1)))	
			f.write("\n")		
			f.write("HOLE2_cv")
			if (HAVE_HOLE):			
				for i in range(HOLE2_nfa):
					f.write(" "+str(np.int(HOLE2_cv[i])))
			else:
				f.write(" "+str(np.int(-1)))	
			f.write("\n")														
			#.........................................
			HOLE_fa=HOLE_fa.flatten()
			f.write("HOLE_fa")
			if (HAVE_HOLE):
				for i in range(len(HOLE_fa)):
					f.write(" "+str(np.int(HOLE_fa[i])))
			else:
				f.write(" "+str(np.int(-1)))
			f.write("\n")

 			f.write("HOLE1_fa")
			if (HAVE_HOLE):
				for i in range(HOLE1_nfa):
					f.write(" "+str(HOLE1_fa[i]))
			else:
				f.write(" "+str(np.int(-1)))
			f.write("\n")

 			f.write("HOLE2_fa")
			if (HAVE_HOLE):
				for i in range(HOLE2_nfa):
					f.write(" "+str(HOLE2_fa[i]))
			else:
				f.write(" "+str(np.int(-1)))
			f.write("\n")		

 			f.write("HOLE3_fa")
			if (HAVE_HOLE):
				for i in range(HOLE3_nfa):
					f.write(" "+str(HOLE3_fa[i]))
			else:
				f.write(" "+str(np.int(-1)))
			f.write("\n")	

 			f.write("HOLE4_fa")
			if (HAVE_HOLE):
				for i in range(HOLE4_nfa):
					f.write(" "+str(HOLE4_fa[i]))
			else:
				f.write(" "+str(np.int(-1)))
			f.write("\n")	

 			f.write("HOLE5_fa")
			if (HAVE_HOLE):
				for i in range(HOLE5_nfa):
					f.write(" "+str(HOLE5_fa[i]))
			else:
				f.write(" "+str(np.int(-1)))
			f.write("\n")														
			#.........................................		
			f.write("WALL_cv")
			for i in range(WALL_ncv):
				f.write(" "+str(np.int(WALL_cv[i])))
			f.write("\n")
			#.........................................
			f.write("WALL_fa")
			for i in range(WALL_nfa):
				f.write(" "+str(np.int(WALL_fa[i])))
			f.write("\n")


			f.write("C_list")
			cclist=c_list.flatten()
			for i in range(len(cclist)):
				f.write(" "+str((cclist[i])))
			f.write("\n")

			#.........................................
			#.........................................		
			f.write("chunk1")
			for i in range(len(chunk1)):
				f.write(" "+str((chunk1[i])))
			f.write("\n")

			f.write("chunk2")
			for i in range(len(chunk2)):
				f.write(" "+str((chunk2[i])))
			f.write("\n")

			f.write("chunk3")
			for i in range(len(chunk3)):
				f.write(" "+str((chunk3[i])))
			f.write("\n")

			f.write("fa1_list")
			for i in range(len(fa1_list)):
				f.write(" "+str((fa1_list[i])))
			f.write("\n")

			f.write("fa2_list")
			for i in range(len(fa2_list)):
				f.write(" "+str((fa2_list[i])))
			f.write("\n")

			f.write("chunk1_rcv")
			for i in range(len(chunk1_rcv)):
				f.write(" "+str(np.int((chunk1_rcv[i]))))
			f.write("\n")

			f.write("chunk2_lcv")
			for i in range(len(chunk2_lcv)):
				f.write(" "+str(np.int((chunk2_lcv[i]))))
			f.write("\n")

			f.write("chunk2_rcv")
			for i in range(len(chunk2_rcv)):
				f.write(" "+str(np.int((chunk2_rcv[i]))))
			f.write("\n")

			f.write("chunk3_lcv")
			for i in range(len(chunk3_lcv)):
				f.write(" "+str(np.int((chunk3_lcv[i]))))
			f.write("\n")

			f.close()
			print "Out put finished!"
			print "ncv:",ncv,"nfa:",nfa,"nno:",nno	
			print "HOLE_ncv:",len(HOLE_fa)	
			print "HOLE1_ncv:",HOLE1_nfa	
			print "N_chunk1:",len(chunk1)		
			print "N_chunk2:",len(chunk2)
			print "N_chunk3:",len(chunk3)	
			print "N_MPI_fa:",len(fa1_list)																
			print "WALL_ncv:",WALL_ncv	
			print "WALL_nfa:",WALL_nfa														

