import os
import sys
import numpy as np
import scipy as sp 
import scipy.sparse as scysparse
import scipy.sparse.linalg as splinalg
from scipy.io import savemat,loadmat
import pylab as plt
from pdb import set_trace
from scipy.interpolate import griddata
from matplotlib import rc as matplotlibrc
import umesh_reader,copy
import math_tool as ma
import Operator as Op
import bivariate_fit as fit
import time
import uns_plot as uns
from matplotlib.collections import PatchCollection
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import spectrum as sp
#the following enables LaTeX typesetting, which will cause the plotting to take forever..  
matplotlibrc('text.latex', preamble='\usepackage{color}')
matplotlibrc('text',usetex=True)
matplotlibrc('font', family='serif')

Problem1 = True
Problem2 = False

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

	P1_1		= False
	PATCH 	= False# need P1_1
	Config 	= True	
	P1_2 		= False  # need config
	P1_3 		= True   # need config

	filename = ['ccc.msh','test.msh','aaa.msh','ddd.msh','eee.msh','p1_C_tri.msh','p1_Cmix.msh','MeshC.msh']
	NCV=np.zeros(3)
	C_RMS=np.zeros(3)
	Phi_dict={}
	kk=0
	if kk ==0:
		mshfile_fullpath = icemcfd_project_folder + filename[kk]
		part_names, xy_no, xy_fa, xy_cv, noofa, cvofa, faono, faocv, partofa = \
						umesh_reader.read_unstructured_grid(mshfile_fullpath,node_reordering=True)
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
			  if part != 'FLUID':
			    B_fa[i] = ifa
			    i = i+1

			B_fa=list(B_fa.values())			

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

			#Obtain noocv 
			noocv = {}
			for icv, face in enumerate(faocv):
				temp  = list(face)
				temp2 = {}
				for i, iface in enumerate(temp):
					temp2 = list(set(temp2).union(set(noofa[iface])))
				noocv[icv] = temp2
			noocv=list(noocv.values())			

			#Obtain cv of inner notes: cvono_in
			In_no = range(nno)
			for i, b_no in enumerate(B_no):
				In_no.remove(b_no)			

			cvono_in = {}
			for i, i_no in enumerate(In_no):
				temp  = list(faono[i_no])
				temp2 = {}
				for i, iface in enumerate(temp):
					temp2 = list(set(temp2).union(set(cvofa[iface])))
				cvono_in[i_no] = temp2
			cvono_in=list(cvono_in.values())	

			In_fa = range(nfa)
			for i, b_fa in enumerate(B_fa):
				In_fa.remove(b_fa)			

			nno_in=len(In_no)	
			nno_B=len(B_no)

			#Obtain ncx and ncy with the same order of In_no
			#Obtain the surounding note numbers around one note,
			# store in the dictionary Sur_no
			ncx=np.zeros(len(In_no))
			ncy=np.zeros(len(In_no))
			Sur_no ={}	
			lengh=np.zeros(len(In_no))

			for i, ino in enumerate(In_no):			

			  Nc = xy_no[ino]
			  ncx[i] = Nc[0]
			  ncy[i] = Nc[1]			

			  sur_no = []
			  temp   =list([])
			  lengh[i]=len(faono[ino])
			  for inoo, ifa in enumerate(faono[ino]):			

			    temp = list(noofa[ifa])
			    temp.remove(ino)
			    sur_no.append(temp[0])			

			  Sur_no[str(ino)]=sur_no
			  INO= np.argmax(lengh)
			    
			#Obtain the coordinate of surounding note group
			#The coordinates of center point is stored at the last line of xy_sur			

			xy_sur={}
			nfaoin_no=np.zeros(len(In_no))			

			for nc,nsurlist in Sur_no.items():			

			  temp=np.zeros((len(nsurlist)+1,3))
			  temp[len(nsurlist),0]  =int(nc)
			  temp[-1,0]  =int(nc)
			  temp[-1,1:3]=xy_no[int(nc)]
			  for i, nsur in enumerate(nsurlist)  :			

			    temp[i,1:3] = xy_no[nsur]
			    temp[i,0]= nsur			

			  xy_sur[str(nc)]=temp

			#instore all nodes!!!
			#Obtain the surounding note numbers around one note,
			# store in the dictionary SSur_no
			llist=range(nno)
			#set_trace()
			nncx=xy_no[:,0]
			nncy=xy_no[:,1]
			SSur_no ={}	
			lengh=np.zeros(nno)

			for i, ino in enumerate(llist):			
	

			  sur_no = []
			  temp   =list([])
			  lengh[i]=len(faono[ino])
			  for inoo, ifa in enumerate(faono[ino]):			

			    temp = list(noofa[ifa])
			    temp.remove(ino)
			    sur_no.append(temp[0])			

			  SSur_no[str(ino)]=sur_no
			  INO= np.argmax(lengh)
			
			    
			#Obtain the coordinate of surounding note group
			#The coordinates of center point is stored at the last line of xy_sur			

			XXYY_sur={}
			nfaoin_no=np.zeros(len(In_no))			

			for nc,nsurlist in SSur_no.items():			

			  temp=np.zeros((len(nsurlist)+1,3))
			  temp[len(nsurlist),0]  =int(nc)
			  temp[-1,0]  =int(nc)
			  temp[-1,1:3]=xy_no[int(nc)]
			  for i, nsur in enumerate(nsurlist)  :			

			    temp[i,1:3] = xy_no[nsur]
			    temp[i,0]= nsur			

			  XXYY_sur[str(nc)]=temp
			#set_trace()
			############################################################
			####################    Prob   1.1      ####################
			############################################################

		if P1_1 :
			area = np.zeros(ncv)
			c_abs = np.zeros(ncv)
			k     = np.zeros(ncv)
			CFL   = np.zeros(ncv)

			for icv,falist in enumerate(faocv):			
				xx = {}
				yy = {}
				for i, note in enumerate(noocv[icv]):
					xx[i]=xy_no[note,0]
					yy[i]=xy_no[note,1]

				#xx = list(xx.values())
				#yy = list(yy.values())
				xx = np.array(xx.values())
				yy = np.array(yy.values())
				area[icv] = ma.PolyArea(xx,yy)				
				#set_trace()
				#cx = yy*np.exp(-xx)
				#cy = 0.5*yy**2*np.exp(-xx)
				cx=xx/(yy+2)
				cy=-np.log(yy+2)
				#cx=((xx+1.0)/(yy+1.0))**0.5
				#cy=-((yy+1.0)/(xx+1.0))**0.5
				#cx=xx
				#cy=-yy			
				c_abs[icv] = (np.average(cx)**2+np.average(cy)**2)**0.5
				k[icv]=c_abs[icv]/(area[icv]**0.5)

			t1=0.85/np.max(k)
			a=0.85/4*np.min(area)/t1	

			CFL1=4*a*t1/area
			CFL2=k*t1

			for i in range(ncv):
				if CFL1[i]>=CFL2[i]:
					CFL[i]=CFL1[i]
				else:
					CFL[i]=CFL2[i]
			print 't1:',t1,'a:',a,ncv
			set_trace()


			if PATCH:	

				patches=uns.contour_plot(xy_no,noocv)
						
				fig = plt.figure(0,figsize=(fig_width*0.8,fig_height*0.8))
				ax = fig.add_subplot(111,aspect='equal')		

				cmap = mpl.cm.viridis 
				pcoll = PatchCollection(patches, cmap=cmap, alpha=1, lw=0.5)
				pcoll.set_array(np.array(CFL))
				gci=ax.add_collection(pcoll)                        
				ax.set_xlim([np.min(xy_no[:,0]), np.max(xy_no[:,0])])
				ax.set_ylim([np.min(xy_no[:,1]), np.max(xy_no[:,1])])			
	 
				#pcoll.set_clim(-3, 3)
				cbar = plt.colorbar(pcoll)  
				#cbar.set_label(r'$Divergence$',fontsize=1.5*gcafontSize)  
				#cbar.set_ticks(np.linspace(-3,3,3))  
				#cbar.set_ticklabels( ('-3', '0', '3'))  
				cbar.ax.set_ylabel(r'$CFL$', fontsize=1.5*gcafontSize)
				cl = plt.getp(cbar.ax, 'ymajorticklabels')
				plt.setp(cl, fontsize=gcafontSize) 		

				ax.set_ylabel(r'$y$',fontsize=1.5*gcafontSize)
				ax.set_ylabel(r'$y$',fontsize=1.5*gcafontSize)
				plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
				plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)		

				fig.tight_layout()
				fig_name = 'p1_CFL.png'
				figure_path = '../report/figures/'
				fig_fullpath = figure_path + fig_name
				#plt.show()
				plt.savefig(fig_fullpath)
				plt.close()
				print fig_name+' saved!'

				#quiver for c

				x = xy_no[:,0]
				y = xy_no[:,1]

				cx=x/(y+2)
				cy=-np.log(y+2)
	
				X,Y = np.meshgrid(x,y)	

				fig = plt.figure(0,figsize=(fig_width,fig_height))
				ax = fig.add_subplot(111)
				ax.quiver(x,y,cx,cy,scale=20,scale_units='xy',color='r')

				ax.set_xlabel(r'$x$',fontsize=1.5*gcafontSize)
				ax.set_ylabel(r'$y$',fontsize=1.5*gcafontSize)
				plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
				plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
				fig.tight_layout()
				fig_name = 'Quiver.pdf'
				figure_path = '../report/figures/'
				fig_fullpath = figure_path + fig_name
				plt.savefig(fig_fullpath)
				plt.close()
				print fig_name+' saved!'						
			############################################################
			####################    configer      ####################
			############################################################

		if Config:

			low  = 300
			high = 500

			phi_in_new = np.zeros(nno_in)
			phi_in     = np.zeros(nno_in)
			phi_bc = np.zeros(nno_B)

			#for i, ifa in enumerate(B_fa):
			#	temp = list(noofa[ifa])
			#	no0=B_no.index(temp[0])
			#	no1=B_no.index(temp[1])
			#	if partofa[ifa]=='HOT':
			#		phi_bc[no0] = high
			#		phi_bc[no1] = high
			#	else:
			#		phi_bc[no0] = low
			#		phi_bc[no1] = low	
			for i, ifa in enumerate(B_fa):
				temp = list(noofa[ifa])
				no0=B_no.index(temp[0])
				no1=B_no.index(temp[1])
				nn0=temp[0]
				nn1=temp[1]
				#set_trace()
				if partofa[ifa]=='HOT':
					phi_bc[no0] = high
					phi_bc[no1] = high
				else:
					phi_bc[no0] = low
					phi_bc[no1] = low
			xx, Lap_biv_q = Op.Lap_biv_no2no(xy_sur,In_no,B_no,ncx,ncy)


			Gx_no2no,Gy_no2no,Gxq_no2no,Gyq_no2no = Op.Grad_no2no(xy_sur,In_no,B_no)
			#set_trace()
			#Lap_biv=Gx_no2no.dot(Gx_no2no)+Gx_no2no.dot(Gx_no2no)
			Lap_biv,Lap_biv_q = Op.Lap_grad(XXYY_sur,In_no,B_no)


























			############################################################
			####################    Prob   1.2      ####################
			############################################################
		if P1_2 :
			alist=[0.0707900178726,4.87497303266,0.513579563266,0.917301951983,0.359007482382,0.00143140329917,0.000933879990509,0.0129376104109]
			tlist=[0.0413108008701,0.0108974961798,0.012414280846422449,0.0015174696060270686,0.00139552027395,0.0111341786129,0.00859085481688,0.0153519491764]
			a  = alist[kk]
			dt = alist[kk]*15
			cx = 0
			cy = 0
			area = np.zeros(ncv)
			c_abs = np.zeros(ncv)
			k     = np.zeros(ncv)
			CFL   = np.zeros(ncv)			
			# BC treatment
			q1_bc = -a*Lap_biv_q.dot(phi_bc.T)
			q2_bc = Gxq_no2no.dot(phi_bc.T)*cx+Gyq_no2no.dot(phi_bc.T)*cy
			q_bc  = -(q1_bc+q2_bc)

			A = -a*Lap_biv+Gx_no2no*cx+Gy_no2no*cy
			#cx=y*e^(-x)
			#cy=(x,y)=0.5*y^2*e^(-x)	
			#cx = xy_no[:,1]*np.exp(-xy_no[:,0])
			#cy = 0.5*xy_no[:,1]**2*np.exp(-xy_no[:,0])	

			ii=0
			for icv,falist in enumerate(faocv):			
				xx = {}
				yy = {}
				for i, note in enumerate(noocv[icv]):
					xx[i]=xy_no[note,0]
					yy[i]=xy_no[note,1]

				#xx = list(xx.values())
				#yy = list(yy.values())
				xx = np.array(xx.values())
				yy = np.array(yy.values())
				area[icv] = ma.PolyArea(xx,yy)				
				#set_trace()
				#cx = yy*np.exp(-xx)
				#cy = 0.5*yy**2*np.exp(-xx)
				cx=xx/(yy+2)
				cy=-np.log(yy+2)
				c_abs[icv] = (np.average(cx)**2+np.average(cy)**2)**0.5
				k[icv]=c_abs[icv]/(area[icv]**0.5)

			t1=0.85/np.max(k)

			CFL1=4*a*dt/area
			CFL2=k*dt

			for i in range(ncv):
				if CFL1[i]>=CFL2[i]:
					CFL[i]=CFL1[i]
				else:
					CFL[i]=CFL2[i]
			print 't1:',dt
			print 'CFLmax:',np.max(CFL)
			set_trace()

			if False: # Second order
				while (np.sum(np.abs(phi_in_new))==0 or \
				       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.01):#	
					ii=ii+1
					phi_in=phi_in_new*1.0
					#set_trace()
					k1=q_bc-A.dot(phi_in)
					phi_temp_2=k1*dt+phi_in
					k2=q_bc-A.dot(phi_temp_2)
					m = phi_in+dt*(k1+k2)*0.5
					phi_in_new = np.array(m)
					#set_trace()
					print 'iteration:',ii
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_v:',np.max(phi_in_new)


			if False: # First order
				while (np.sum(np.abs(phi_in_new))==0 or \
				       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.01):#	
					ii=ii+1
					phi_in=phi_in_new*1.0
					#set_trace()
					m = q_bc*dt+(np.eye(len(In_no))-dt*A).dot(phi_in)
					phi_in_new = np.array(m)
					phi_in_new = phi_in_new[0,:]
					print 'iteration:',ii
					#print "phi_new i:s:",phi_in_new
					#print "phi -old is: ", phi_in
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_v:',np.max(phi_in_new)

			if False: # Forth order
				while (np.sum(np.abs(phi_in_new))==0 or \
				       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.001):#	
					ii=ii+1
					phi_in=phi_in_new*1.0
					#set_trace()
					k1=q_bc-A.dot(phi_in)

					phi_temp_2=k1*dt*0.5+phi_in
					k2=q_bc-A.dot(phi_temp_2)


					phi_temp_3=k2*dt*0.5+phi_in
					k3=q_bc-A.dot(phi_temp_3)

					phi_temp_4=k3*dt+phi_in
					k4=q_bc-A.dot(phi_temp_4)

					m = phi_in+dt*(k1+k2*2.0+k3*2.0+k4)/6.0
					phi_in_new = np.array(m)
					#set_trace()
					print 'iteration:',ii
					#print "phi_new i:s:",phi_in_new
					#print "phi -old is: ", phi_in
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_v:',np.max(phi_in_new)
					print 'min_v:',np.min(phi_in_new)
					aa= np.argmax(phi_in_new)
					no=In_no[aa]
					print xy_no[no,:]



			time1=0
			if False: # Second order
				while (np.sum(np.abs(phi_in_new))==0 or \
				       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.01):#	
					ii=ii+1
					phi_in=phi_in_new*1.0
					#set_trace()
					k1=q_bc-A.dot(phi_in)
					phi_temp_2=k1*dt+phi_in
					k2=q_bc-A.dot(phi_temp_2)
					m = phi_in+dt*(k1+k2)*0.5
					phi_in_new = np.array(m)
					time1=time1+dt
					#set_trace()
					print 'iteration:',ii
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_v:',np.max(phi_in_new)

			time2=0
			if True: # CRANK-NICHOLSON
				while (np.sum(np.abs(phi_in_new))==0 or \
				       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.001):#	
					ii=ii+1
					time2=time2+dt
					phi_in=phi_in_new*1.0
					RHS=(np.eye(A.shape[0])-0.5*A*dt).dot(phi_in)+dt*q_bc
					LHS=np.eye(A.shape[0])+0.5*A*dt
					#set_trace()
					phi_in_new=splinalg.spsolve(LHS,RHS.T, permc_spec=None, use_umfpack=True)

					#m = phi_in+dt*(k1+k2*2.0+k3*2.0+k4)/6.0
					#phi_in_new = np.array(m)
					#set_trace()
					print 'iteration:',ii
					#print "phi_new i:s:",phi_in_new
					#print "phi -old is: ", phi_in
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_v:',np.max(phi_in_new)
					print 'min_v:',np.min(phi_in_new)
					aa= np.argmax(phi_in_new)
					no=In_no[aa]
					print xy_no[no,:]		
					##########################################
					#########################################
					##########################################
					P_phi=np.zeros(nno)	
					for i in range(nno):
						#if (ii==2):
						#		set_trace()							
						if (i in In_no):
							k=In_no.index(i)
							P_phi[i]=phi_in[k]
						else:
							k=B_no.index(i)
							P_phi[i]=phi_bc[k]	
					print '**************'
					print np.max(P_phi)
					print np.min(P_phi)
					#################################################################################
					#################################################################################
					#################################################################################
					if (ii==1 or ii%15==0):
					#if False:
						#set_trace()
						x = xy_no[:,0]
						y = xy_no[:,1]#
						phi = P_phi					
						n=50
						xg = np.linspace(x.min(),x.max(),n)
						yg = np.linspace(y.min(),y.max(),n)
						X,Y = np.meshgrid(xg,yg)
						# interpolate Z values on defined grid
						Z = griddata(np.vstack((x.flatten(),y.flatten())).T, phi.flatten().T, \
								        (X,Y), method='cubic').reshape(X.shape)
						# mask nan values, so they will not appear on plot
						Zm = np.ma.masked_where(np.isnan(Z),Z)	

						center1 = np.array([3.0, 7.0])
						radius1 = 1
						DIST = ((X - center1[0])**2.0 + (Y - center1[1])**2.0)**0.5
						Z[np.where(DIST < radius1)] = np.nan		
						center2 = np.array([7.0, 7.0])				

						radius1 = 1
						DIST = ((X - center2[0])**2.0 + (Y - center2[1])**2.0)**0.5
						Z[np.where(DIST < radius1)] = np.nan						

						Z[np.where((X>2) & (X <8) & (Y>2) & (Y<3)) ]= np.nan		
						Zm = np.ma.masked_where(np.isnan(Z),Z)	

						fig = plt.figure(0,figsize=(fig_width,fig_height))
						ax = fig.add_subplot(111, aspect='equal')		
						c=ax.pcolormesh(X,Y,Zm,shading='gouraud')	
						for inos_of_fa in noofa:
							ax.plot(xy_no[inos_of_fa,0], xy_no[inos_of_fa,1], 'k-', linewidth = 0.2*lineWidth)
						cbar= fig.colorbar(c)
						cbar.ax.tick_params(labelsize=gcafontSize)							

						c.set_clim(0, 8)
				 
						cbar.ax.set_ylabel(r'$/phi$', fontsize=1.2*gcafontSize)
						cl = plt.getp(cbar.ax, 'ymajorticklabels')
						plt.setp(cl, fontsize=gcafontSize) 				

						ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
						ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
						#plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
						#plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)	
						k='time= '+'%.2f' % (time2-dt)+'s'
						ax.set_title(k, fontsize=1.2*gcafontSize)
						ax.set_xticks([])
						ax.set_yticks([])
						fig.tight_layout()
						fig_name = '4diff'+str(ii)+'.png'
						figure_path = '../report/figures/test/'
						fig_fullpath = figure_path + fig_name
						plt.savefig(fig_fullpath)
						plt.close()
						print fig_name+' saved!'	






			############################################################
			####################    Prob   1.3      ####################
			############################################################
		if P1_3 :
			alist=[0.0727433976688,4.87497303266,0.513579563266,0.917301951983,0.359007482382,0.00143140329917,0.000933879990509,0.0129376104109]
			tlist=[0.0402014811742,0.0108974961798,0.012414280846422449,0.0015174696060270686,0.00139552027395,0.0111341786129,0.00859085481688,0.0153519491764]
			a  = 0
			dt = alist[kk]*0.9

			#cx=ncx/(ncy+2)
			#cy=-np.log(ncy+2)			

			cx=ncx/(ncy+2)
			cy=-np.log(ncy+2)
			area = np.zeros(ncv)
			c_abs = np.zeros(ncv)
			k     = np.zeros(ncv)
			CFL   = np.zeros(ncv)		

			for i, ino in enumerate(In_no):
				phi_in[i]=0

			phi_in_new=phi_in

			#GX,GY,GXq,GYq=Op.Grad_no2no_upwind(xy_sur,In_no,B_no,cx,cy)

			# BC treatment
			#q1_bc = -a*Lap_biv_q.dot(phi_bc.T)
			#q2_bc = GXq.dot(phi_bc.T)*1+GYq.dot(phi_bc.T)*1
			#q_bc  = -q2_bc

			#A = GX*cx+GY*cy

			ii=0
			for icv,falist in enumerate(faocv):			
				xx = {}
				yy = {}
				for i, note in enumerate(noocv[icv]):
					xx[i]=xy_no[note,0]
					yy[i]=xy_no[note,1]

				#xx = list(xx.values())
				#yy = list(yy.values())
				xx = np.array(xx.values())
				yy = np.array(yy.values())
				area[icv] = ma.PolyArea(xx,yy)				
				#set_trace()
				#cx = yy*np.exp(-xx)
				#cy = 0.5*yy**2*np.exp(-xx)
				#cx=xx/(yy+2)
				#cy=-np.log(yy+2)
				cx=xx/(yy+2)
				cy=-np.log(yy+2)		
				c_abs[icv] = (np.average(cx)**2+np.average(cy)**2)**0.5
				k[icv]=c_abs[icv]/(area[icv]**0.5)

			t1=0.85/np.max(k)

			CFL1=4*a*dt/area
			CFL2=k*dt

			for i in range(ncv):
				if CFL1[i]>=CFL2[i]:
					CFL[i]=CFL1[i]
				else:
					CFL[i]=CFL2[i]
			print 't1:',dt
			print 'CFLmax:',np.max(CFL)


		
			if False: # Impliciy
				while (ii<1000000):#	
					ii=ii+1
					phi_in=phi_in_new*1.0

					#set_trace()
					phi_in_new=splinalg.spsolve(A*dt+np.eye(len(phi_in)),phi_in+q_bc, permc_spec=None, use_umfpack=True)
					a,b=np.linalg.eig(A*dt+np.eye(len(phi_in)))
					print np.max(a)
					set_trace()

					print 'iteration:',ii
					#print "phi_new i:s:",phi_in_new
					#print "phi -old is: ", phi_in
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_v:',np.max(phi_in_new)
					print 'min_v:',np.min(phi_in_new)
					aa= np.argmax(phi_in_new)
					no=In_no[aa]
					print xy_no[no,:]		

			time1=0
			if False: # Second order

				while (np.sum(np.abs(phi_in_new))==0 or \
				       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.01):#	
					ii=ii+1
					phi_in=phi_in_new*1.0
					#set_trace()
					k1=q_bc-A.dot(phi_in.T)
					k1=np.array(k1)
					k1=k1[0,:]
					phi_temp_2=k1*dt+phi_in
					#set_trace()
					k2=q_bc-A.dot(phi_temp_2)
					k2=np.array(k2)[0,:]
					m = phi_in+dt*(k1+k2)*0.5
					phi_in_new = np.array(m)
					time1=time1+dt
					#set_trace()
					print 'iteration:',ii
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_v:',np.max(phi_in_new)

					print 'iteration:',ii
			if False: # CRANK-NICHOLSON
				while (np.sum(np.abs(phi_in_new))==0 or \
				       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.001):#	
					ii=ii+1
					phi_in=phi_in_new*1.0
					RHS=(np.eye(A.shape[0])-0.5*A*dt).dot(phi_in)+dt*q_bc
					LHS=np.eye(A.shape[0])+0.5*A*dt
					#set_trace()
					phi_in_new=splinalg.spsolve(LHS,RHS.T, permc_spec=None, use_umfpack=True)

					#m = phi_in+dt*(k1+k2*2.0+k3*2.0+k4)/6.0
					#phi_in_new = np.array(m)
					#set_trace()
					print 'iteration:',ii
					#print "phi_new i:s:",phi_in_new
					#print "phi -old is: ", phi_in
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_v:',np.max(phi_in_new)
					print 'min_v:',np.min(phi_in_new)
					aa= np.argmax(phi_in_new)
					no=In_no[aa]
					print xy_no[no,:]	
			
			if True: # LAGR
				cx=ncx/(ncy+2)
				cy=-np.log(ncy+2)			
				newx=ncx-cx*dt
				newy=ncy-cy*dt
				X,Y = np.meshgrid(newx,newy)
				x = xy_no[:,0]
				y = xy_no[:,1]#			
				set_trace()
				ttime=0
				while (ii<1000000):
					ttime=ttime+dt
					ii=ii+1
					phi_in=phi_in_new*1.0

					P_phi=np.zeros(nno)	
					for i in range(nno):	
						if (i in In_no):
							k=In_no.index(i)
							P_phi[i]=phi_in[k]
						else:
							k=B_no.index(i)
							P_phi[i]=phi_bc[k]
					#set_trace()
					Z = griddata(np.vstack((x.flatten(),y.flatten())).T, P_phi.flatten().T, \
		         		     (newx,newy), method='linear').reshape(newx.shape)		
					phi_in_new=Z
					#set_trace()
					#AA,bb=Op.Semi_Lagr_upwind(xy_sur,In_no,B_no,cx,cy,dt)

					#phi_in_new=AA.dot(phi_in)+bb.dot(phi_bc)
					#set_trace()

					print 'iteration:',ii
					#print "phi_new i:s:",phi_in_new
					#print "phi -old is: ", phi_in
					print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
					print 'max_phi:',np.max(phi_in_new)
					print 'min_phi:',np.min(phi_in_new)
				##########################################
					#########################################
					##########################################
					P_phi=np.zeros(nno)	
					for i in range(nno):																	
						if (i in In_no):
							k=In_no.index(i)
							P_phi[i]=phi_in[k]
						else:
							k=B_no.index(i)
							P_phi[i]=phi_bc[k]	
					print '**************'
					print np.max(P_phi)
					print np.min(P_phi)
					#################################################################################
					#################################################################################
					#################################################################################
					#if (ii==1 or ii%15==0):
					if False:
						#set_trace()
						x = xy_no[:,0]
						y = xy_no[:,1]#
						phi = P_phi					
						n=50
						xg = np.linspace(x.min(),x.max(),n)
						yg = np.linspace(y.min(),y.max(),n)
						X,Y = np.meshgrid(xg,yg)
						# interpolate Z values on defined grid
						Z = griddata(np.vstack((x.flatten(),y.flatten())).T, phi.flatten().T, \
								        (X,Y), method='cubic').reshape(X.shape)
						# mask nan values, so they will not appear on plot
						Zm = np.ma.masked_where(np.isnan(Z),Z)	

						center1 = np.array([3.0, 7.0])
						radius1 = 1
						DIST = ((X - center1[0])**2.0 + (Y - center1[1])**2.0)**0.5
						Z[np.where(DIST < radius1)] = np.nan		
						center2 = np.array([7.0, 7.0])				

						radius1 = 1
						DIST = ((X - center2[0])**2.0 + (Y - center2[1])**2.0)**0.5
						Z[np.where(DIST < radius1)] = np.nan						

						Z[np.where((X>2) & (X <8) & (Y>2) & (Y<3)) ]= np.nan		
						Zm = np.ma.masked_where(np.isnan(Z),Z)	

						fig = plt.figure(0,figsize=(fig_width,fig_height))
						ax = fig.add_subplot(111, aspect='equal')		
						c=ax.pcolormesh(X,Y,Zm,shading='gouraud')	
						for inos_of_fa in noofa:
							ax.plot(xy_no[inos_of_fa,0], xy_no[inos_of_fa,1], 'k-', linewidth = 0.2*lineWidth)
						cbar= fig.colorbar(c)
						cbar.ax.tick_params(labelsize=gcafontSize)							

						c.set_clim(-50, 500)
				 
						cbar.ax.set_ylabel(r'$/phi$', fontsize=1.2*gcafontSize)
						cl = plt.getp(cbar.ax, 'ymajorticklabels')
						plt.setp(cl, fontsize=gcafontSize) 				

						ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
						ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
						#plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
						#plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)	
						k='time= '+str(round(ttime,3))+'s'
						ax.set_title(k, fontsize=1.2*gcafontSize)
						ax.set_xticks([])
						ax.set_yticks([])
						fig.tight_layout()
						fig_name = '3conv'+str(ii)+'.png'
						figure_path = '../report/figures/'
						fig_fullpath = figure_path + fig_name
						plt.savefig(fig_fullpath)
						plt.close()
						print fig_name+' saved!'				




####################################################################################
####################################  Problem2  ####################################
####################################################################################
if Problem2:
	#basic func
	if True:
		# base function for any single value
		def base(xg,x):
			h=np.ones((len(xg)))
			for j in range(len(xg)):
				for i in range(len(xg)):
					if (i!=j):
						m=(x-xg[i])/(xg[j]-xg[i])
						h[j]=h[j]*m
			return h
		# first order base func need orign y_value		

		def der_base(xl,x):
			k={}
			h=np.ones((len(xl)))		

			temp=np.zeros(len(xl))		

			for j in range(len(xl)):
				b=np.ones((len(xl)-1))		

				xx=np.ones((len(xl)-1))*x
				for i in range(len(xl)):
					k=list(copy.deepcopy(xl))
					bb=1
					k.remove(xl[j])
					b=xl[j]-k
					t=xx-k
					tt=np.ones((len(xl)-1))			
					for z,value in enumerate(t):
						bb=bb*b[z]#over
						for zz in range(len(t)):
							if (zz!=z):
								tt[z]=tt[z]*t[zz]
				temp[j]=np.sum(tt)
				h[j]=temp[j]/bb
				#set_trace()
			return h		

		def der_base2_0(xl,x):
			k={}
			#h=np.ones((len(xl)))		

			#temp=np.zeros(len(xl))
			j=0
			if j==0:
				b=np.ones((len(xl)-1))		

				xx=np.ones((len(xl)-1))*x
				for i in range(len(xl)):
					k=list(copy.deepcopy(xl))
					bb=1
					k.remove(xl[j])
					b=xl[j]-k
					t=xx-k
					tt=np.ones((len(xl)-1))			
					for z,value in enumerate(t):
						bb=bb*b[z]#over
						for zz in range(len(t)):
							if (zz!=z):
								tt[z]=tt[z]*t[zz]
				temp=np.sum(tt)
				h=temp/bb
				#set_trace()
			return h		

		def x_g2l(xg,xl):
			I=np.zeros((len(xl),len(xg)))
			for i,ixl in enumerate(xl):
				k=xl[i]
				I[i,:]=sp.base(xg,k)
			return I							

	mesh=False
	diff=False
	contour=False
	conv=False
	validation=False

	Nj=15
	N =Nj-1
	base=np.arange(Nj)*1.0
	base2=np.arange(Nj-1)*1.0
	xl=0.5*(1-np.cos(np.pi*base/N))
	xg=0.5*(1-np.cos(np.pi*(2*base2+1)/(2*N)))
	ug=np.zeros((N,N))
	ul=np.zeros((N+1,N+1))

	GX,GY=np.meshgrid(xg,xg)
	LX,LY=np.meshgrid(xl,xl)

	l_flux=np.zeros((len(ug),len(ul)))
	for i in range(len(ug)):
		l_flux[i]=der_base(xl,xg[i])	
	
	g_flux=np.zeros((len(ug),len(ul)))
	for i in range(len(ug)):
		g_flux[i]=der_base(xl,xg[i])

	lxx_flux=np.zeros((len(ug),len(ul)))
	for i in range(len(ug)):
		lxx_flux[i]=der_base(xl,xg[i])	

	hxx_flux=np.zeros((len(ug),len(ug)))
	for i in range(len(ug)):
		hxx_flux[i]=der_base(xg,xg[i])	


	xdelt=x_g2l(xg,xl)
	

	if mesh:
		fig = plt.figure(0,figsize=(fig_width*0.8,fig_height*0.8))
		ax = fig.add_subplot(111)	

		for i in range(GX.shape[0]):
			ax.axvline(GY[i,0],color='k',linestyle='-')
			ax.axhline(GX[0,i],color='k',linestyle='-')
		for i in range(LX.shape[0]):
			ax.axvline(LY[i,0],color='r',linestyle='--')
			ax.axhline(LX[0,i],color='r',linestyle='--')
		ax.set_xlim([0.0, 1.0])
		ax.set_ylim([0.0, 1.0])
		plt.xlabel(r'$x$',fontsize=1.5*gcafontSize)
		plt.ylabel(r'$y$',fontsize=1.5*gcafontSize)
		fig_name ='mesh.png'
		figure_path = '../report/figures/'
		fig_fullpath = figure_path + fig_name
		plt.savefig(fig_fullpath)
		plt.close()
		print fig_name+'saved!'
	####################################
	####################################
	####################################
	#if validation:
	if True:
		ugx=np.zeros((len(ug),len(ug)))
		ugx=np.sin(GX*2*np.pi)*np.sin(GY*2*np.pi)*0.1
		#set_trace()
		phi_in=ugx
		phi_in_new=ugx	

		ii=0
		dx = xg[1:] - xg[:-1]
		dx_min = np.min(dx)
		dt = 0.1*dx_min*0.5
		#dt = 0.1*dx_min**2
		#while (ii==0 or \
		#       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.01):#	
		time = 0.0
		print dt
		while (ii<10000):
			ii=ii+1
			phi_in=phi_in_new*1.0	

			ulx=phi_in.dot(xdelt.T)
			uly=xdelt.dot(phi_in)	

			#ulx[:,0] = np.sin(5*time)
			#uly[0,:] = np.sin(5000*time)*30.0
			#ulx[:,0] = np.sin(5000*time)*30.0
			ulx[:,0] = ulx[:,-1]
			uly[0,:] = uly[-1,:]
			
			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)
			Fx=F.dot(hxx_flux.T)
			Qy=hxx_flux.dot(Q)			

			m = phi_in-dt*(F+Q)
			phi_in_new = np.array(m)	
		
			time += dt
			print 'iteration:',ii
			print np.max(phi_in_new)
			print np.min(phi_in_new)
			#l_x
			LY_x,LX_x=np.meshgrid(xl,xg)
			LY_y,LX_y=np.meshgrid(xg,xl)

			fig = plt.figure(0,figsize=(fig_width*0.8,fig_height*0.6))
			ax = fig.add_subplot(111)
			ax = Axes3D(fig)	
			ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
			ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
			ax.set_zlabel(r'$\phi$',fontsize=1.2*gcafontSize)
			ax.plot_surface(GX, GY,phi_in_new, rstride=1, cstride=1, cmap=cm.viridis)
			fig_name ='validation0.png'
			figure_path = '../report/figures/'
			fig_fullpath = figure_path + fig_name
			plt.savefig(fig_fullpath)
			plt.close()

			fig = plt.figure(0,figsize=(fig_width*0.8,fig_height*0.6))
			ax = fig.add_subplot(111)
			ax = Axes3D(fig)	
			ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
			ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
			ax.set_zlabel(r'$\phi_x$',fontsize=1.2*gcafontSize)
			ax.plot_surface(GX, GY,F, rstride=1, cstride=1, cmap=cm.viridis)
			fig_name ='validation1.png'
			figure_path = '../report/figures/'
			fig_fullpath = figure_path + fig_name
			plt.savefig(fig_fullpath)
			plt.close()

			fig = plt.figure(0,figsize=(fig_width*0.8,fig_height*0.6))
			ax = fig.add_subplot(111)
			ax = Axes3D(fig)	
			ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
			ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
			ax.set_zlabel(r'$\phi_y$',fontsize=1.2*gcafontSize)
			ax.plot_surface(GX, GY,Q, rstride=1, cstride=1, cmap=cm.viridis)
			fig_name ='validation2.png'
			figure_path = '../report/figures/'
			fig_fullpath = figure_path + fig_name
			plt.savefig(fig_fullpath)
			plt.close()

			fig = plt.figure(0,figsize=(fig_width*0.8,fig_height*0.6))
			ax = fig.add_subplot(111)
			ax = Axes3D(fig)	
			ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
			ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
			ax.set_zlabel(r'$\phi_{xx}$',fontsize=1.2*gcafontSize)
			ax.plot_surface(GX, GY,Fx, rstride=1, cstride=1, cmap=cm.viridis)
			fig_name ='validation3.png'
			figure_path = '../report/figures/'
			fig_fullpath = figure_path + fig_name
			plt.savefig(fig_fullpath)
			plt.close()		

			fig = plt.figure(0,figsize=(fig_width*0.8,fig_height*0.6))
			ax = fig.add_subplot(111)
			ax = Axes3D(fig)	
			ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
			ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
			ax.set_zlabel(r'$\phi_{yy}$',fontsize=1.2*gcafontSize)
			ax.plot_surface(GX, GY,Qy, rstride=1, cstride=1, cmap=cm.viridis)
			fig_name ='validation4.png'
			figure_path = '../report/figures/'
			fig_fullpath = figure_path + fig_name
			plt.savefig(fig_fullpath)
			plt.close()

			set_trace()


			ax.plot_surface(GX, GY,Fx, rstride=1, cstride=1, cmap=cm.viridis)
			#plt.ylim([-2.0, 2.0])
			fig_name ='diffu'+str(ii)+'_2D.png'
			figure_path = '../report/figures/test/'
			fig_fullpath = figure_path + fig_name
			plt.savefig(fig_fullpath)
			plt.close()
			print fig_name+'saved!'
			if False:
				fig = plt.figure()
				ax = Axes3D(fig)		
				#set_trace()
				print '000'	
				ax.plot_surface(GX, GY,phi_in, rstride=1, cstride=1, cmap=cm.viridis)
				ax.plot_surface(GX, np.zeros(GY.shape),np.zeros(GX.shape)+phi_in[:,-5], rstride=1, cstride=1, color='k')

				
				#plt.ylim([-2.0, 2.0])
				fig_name ='diffu'+str(ii)+'_2D.png'
				figure_path = '../report/figures/'
				fig_fullpath = figure_path + fig_name
				#if ii ==20:
					#set_trace()
					#plt.show()
				plt.savefig(fig_fullpath)
				plt.close()
				
				print fig_name+'saved!'			
				set_trace()		




	####################################
	####################################
	####################################

	if conv:
		#ugx=np.ones((len(ug),len(ug)))
		ugx=np.sin(GX*2*np.pi)*np.sin(GY*2*np.pi)
		#set_trace()
		phi_in=ugx
		phi_in_new=ugx	

		ii=0
		dx = xg[1:] - xg[:-1]
		dx_min = np.min(dx)
		
		LY_x,LX_x=np.meshgrid(xl,xg)
		LY_y,LX_y=np.meshgrid(xg,xl)
		cx=GX**2*GY*0.001
		cy=GY**2*GX*0.001			
		#cx=np.ones(GX.shape)
		#cy=np.ones(GX.shape)
		C=(cx**2+cy**2)**0.5
		dt = 0.8/np.max(C)*dx_min*3.7
		print 'max_c:',np.max(C)
		print 'cfl=',np.max(C)/dx_min*dt
		print 'dt:',dt
		#dt = 0.1*dx_min**2
		#while (ii==0 or \
		#       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.01):#	
		time = 0.0
		set_trace()
		while (ii<10000):
			ii=ii+1
			phi_in=phi_in_new*1.0	

			ulx=phi_in.dot(xdelt.T)
			uly=xdelt.dot(phi_in)	

			#ulx[:,0] = np.sin(5*time)
			#uly[0,:] = np.sin(5000*time)*30.0
			#ulx[:,0] = np.sin(5000*time)*30.0

			#ulx[:,0] = np.cos(LY_x[:,0])
			#uly[0,:] = np.cos(LX_y[0,:])
			ulx[:,0] = ulx[:,-1]
			uly[0,:] = uly[-1,:]
			#ulx[:,0] = 3*np.ones(ulx[:,0].shape)
			#uly[0,:] = 3*np.ones(uly[0,:].shape)
			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)
			#Fx=F.dot(hxx_flux.T)
			#Qy=hxx_flux.dot(Q)			
			#k1=phi_in-dt*(F*cx+Q*cy)
			k1=-(F*cx+Q*cy)

			phi_temp_2=k1*dt*0.5+phi_in			
			ulx=phi_temp_2.dot(xdelt.T)
			uly=xdelt.dot(phi_temp_2)			
			#ulx[:,0] = np.cos(LY_x[:,0])*np.cos(LX_x[:,0])
			#uly[0,:] = np.cos(LY_y[0,:])*np.cos(LX_y[0,:])			
			ulx[:,0] = ulx[:,-1]
			uly[0,:] = uly[-1,:]		
			#ulx[:,0] = 3*np.ones(ulx[:,0].shape)
			#uly[0,:] = 3*np.ones(uly[0,:].shape)				
			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)			
			k2=-(F*cx+Q*cy)

			phi_temp_3=k2*dt*0.5+phi_in			
			ulx=phi_temp_3.dot(xdelt.T)
			uly=xdelt.dot(phi_temp_3)			
			#ulx[:,0] = np.cos(LY_x[:,0])*np.cos(LX_x[:,0])
			#uly[0,:] = np.cos(LY_y[0,:])*np.cos(LX_y[0,:])	
			ulx[:,0] = ulx[:,-1]
			uly[0,:] = uly[-1,:]		
			#ulx[:,0] = 3*np.ones(ulx[:,0].shape)
			#uly[0,:] = 3*np.ones(uly[0,:].shape)						
			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)			
			k3=-(F*cx+Q*cy)			

			phi_temp_4=k3*dt+phi_in
			ulx=phi_temp_4.dot(xdelt.T)
			uly=xdelt.dot(phi_temp_4)			
			#ulx[:,0] = np.cos(LY_x[:,0])*np.cos(LX_x[:,0])
			#uly[0,:] = np.cos(LY_y[0,:])*np.cos(LX_y[0,:])			
			ulx[:,0] = ulx[:,-1]
			uly[0,:] = uly[-1,:]
			#ulx[:,0] = 3*np.ones(ulx[:,0].shape)
			#uly[0,:] = 3*np.ones(uly[0,:].shape)			
			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)				
			k4=-(F*cx+Q*cy)	

			m = phi_in+dt*(k1+k2*2.0+k3*2.0+k4)/6.0
			phi_in_new = np.array(m)

			#set_trace()
			print 'iteration:',ii
			#print "phi_new i:s:",phi_in_new
			#print "phi -old is: ", phi_in
			print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
			print 'max_v:',np.max(phi_in_new)
			print 'min_v:',np.min(phi_in_new)
			aa= np.argmax(phi_in_new)
			time += dt


			if False:	
			#if (ii==1) or (ii%1000==0):
				fig = plt.figure(0,figsize=(fig_width,fig_height))
				ax = fig.add_subplot(111, aspect='equal')		
				c=ax.pcolormesh(GX,GY,phi_in_new,shading='gouraud')			

				cbar= fig.colorbar(c)
				cbar.ax.tick_params(labelsize=gcafontSize)									

				c.set_clim(0, 1)
						 
				cbar.ax.set_ylabel(r'$/phi$', fontsize=1.2*gcafontSize)
				cl = plt.getp(cbar.ax, 'ymajorticklabels')
				plt.setp(cl, fontsize=gcafontSize) 						

				ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
				ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
				#plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
				#plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)	
				k='time= '+str(round(time-dt,3))+'s'
				ax.set_title(k, fontsize=1.2*gcafontSize)
				ax.set_xticks([])
				ax.set_yticks([])
				fig.tight_layout()
				fig_name = 'sp_conv'+str(ii)+'.png'
				figure_path = '../report/figures/'
				fig_fullpath = figure_path + fig_name
				plt.savefig(fig_fullpath)
				plt.close()
				print fig_name+' saved!'						
	####################################
	####################################
	####################################
	if diff:
		ugx=np.zeros((len(ug),len(ug)))
		#ugx=np.sin(GX*2*np.pi)*np.sin(GY*2*np.pi)*0.1
		#set_trace()
		phi_in=ugx
		phi_in_new=ugx	

		ii=0
		dx = xg[1:] - xg[:-1]
		dx_min = np.min(dx)
		LY_x,LX_x=np.meshgrid(xl,xg)
		LY_y,LX_y=np.meshgrid(xg,xl)
		cx=GX**2*GY*0.001
		cy=GY**2*GX*0.001			
		C=(cx**2+cy**2)**0.5
		dtt = 0.8/np.max(C)*dx_min
		dt=0.5*dtt*3.5
		vis=0.8/4*dx_min**2/dtt
		print 'max_c:',np.max(C)
		print 'cfl=',4*vis*dt/dx_min**2
		print 'dt:',dt
		print 'vis:',vis
		#dt = 0.1*dx_min**2
		#while (ii==0 or \
		#       np.sum(np.abs(phi_in - phi_in_new)**2)**0.5>0.01):#	
		time = 0.0
		set_trace()
		while (ii<10000):
			ii=ii+1
			phi_in=phi_in_new*1.0	
			ulx=phi_in.dot(xdelt.T)
			uly=xdelt.dot(phi_in)			
			ulx[:,0] = 3*np.ones(ulx[:,0].shape)
			uly[0,:] = 5*np.ones(uly[0,:].shape)
			ulx[:,-1] = ulx[:,0]
			uly[-1,:] = uly[0,:]	


			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)
			Fx=F.dot(hxx_flux.T)
			Qy=hxx_flux.dot(Q)			
			#k1=phi_in+dt*(Fx*vis+Qy*vis)
			k1=(Fx*vis+Qy*vis)

			phi_temp_2=k1*dt*0.5+phi_in			
			ulx=phi_temp_2.dot(xdelt.T)
			uly=xdelt.dot(phi_temp_2)			

			ulx[:,0] = 3*np.ones(ulx[:,0].shape)
			uly[0,:] = 5*np.ones(uly[0,:].shape)
			ulx[:,-1] = ulx[:,0]
			uly[-1,:] = uly[0,:]	
			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)			
			Fx=F.dot(hxx_flux.T)
			Qy=hxx_flux.dot(Q)	
			k2=(Fx*vis+Qy*vis)


			phi_temp_3=k2*dt*0.5+phi_in			
			ulx=phi_temp_3.dot(xdelt.T)
			uly=xdelt.dot(phi_temp_3)			
			ulx[:,0] = 3*np.ones(ulx[:,0].shape)
			uly[0,:] = 5*np.ones(uly[0,:].shape)
			ulx[:,-1] = ulx[:,0]
			uly[-1,:] = uly[0,:]	
			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)	
			Fx=F.dot(hxx_flux.T)
			Qy=hxx_flux.dot(Q)						
			k3=(Fx*vis+Qy*vis)	


			phi_temp_4=k3*dt+phi_in
			ulx=phi_temp_4.dot(xdelt.T)
			uly=xdelt.dot(phi_temp_4)			
			ulx[:,0] = 3*np.ones(ulx[:,0].shape)
			uly[0,:] = 5*np.ones(uly[0,:].shape)
			ulx[:,-1] = ulx[:,0]
			uly[-1,:] = uly[0,:]	
			F=ulx.dot(lxx_flux.T)
			Q=lxx_flux.dot(uly)		
			Fx=F.dot(hxx_flux.T)
			Qy=hxx_flux.dot(Q)						
			k4=(Fx*vis+Qy*vis)	

			m = phi_in+dt*(k1+k2*2.0+k3*2.0+k4)/6.0
			phi_in_new = np.array(m)
			#set_trace()
			print 'iteration:',ii
			#print "phi_new i:s:",phi_in_new
			#print "phi -old is: ", phi_in
			print 'difference:',np.sum(np.abs(phi_in - phi_in_new)**2)**0.5
			print 'max_v:',np.max(phi_in_new)
			print 'min_v:',np.min(phi_in_new)
			aa= np.argmax(phi_in_new)
			time += dt


			if False:
			#if (ii==1) or (ii%3000==0):
				fig = plt.figure(0,figsize=(fig_width,fig_height))
				ax = fig.add_subplot(111, aspect='equal')		
				c=ax.pcolormesh(GX,GY,phi_in_new,shading='gouraud')			

				cbar= fig.colorbar(c)
				cbar.ax.tick_params(labelsize=gcafontSize)									

				c.set_clim(0, 5)
						 
				cbar.ax.set_ylabel(r'$/phi$', fontsize=1.2*gcafontSize)
				cl = plt.getp(cbar.ax, 'ymajorticklabels')
				plt.setp(cl, fontsize=gcafontSize) 						

				ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
				ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
				k='time= '+str(round(time-dt,3))+'s'
				ax.set_title(k, fontsize=1.2*gcafontSize)			
				ax.set_xticks([])
				ax.set_yticks([])
				fig.tight_layout()
				fig_name = 'sp_diff'+str(ii)+'.png'
				figure_path = '../report/figures/'
				fig_fullpath = figure_path + fig_name
				plt.savefig(fig_fullpath)
				plt.close()
				print fig_name+' saved!'	
