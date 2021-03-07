import sys
import os
import numpy as np 
import math as ma
import scipy
import copy
from pdb import set_trace
import matplotlib.pyplot as plt
from matplotlib import rc as matplotlibrc
from scipy.interpolate import griddata
import matplotlib.tri as tri
import umesh_reader








ONE_D=0
TWO_D=1
TEST_MAP=0
TEST_MAP2=0
TEST_2D =0

fig_name = 'one.png'
fig_width = 6
fig_height = 6
textFontSize   = 10
gcafontSize    = 20
lineWidth      = 2  
if TWO_D:
	fig_name = 'one.png'
	database={}

	database['x']=np.loadtxt('./temp/mesh_x.dat')
	database['y']=np.loadtxt('./temp/mesh_y.dat')
	database['phi']=np.loadtxt('./temp/solution_5000.dat')

	x = database['x']
	y = database['y']
	phi = database['phi']
	triang = tri.Triangulation(x.flatten().T,y.flatten().T)
#	plt.tricontourf(triang, phi.flatten().T)
#	plt.show()
	set_trace()	
	XX, YY = np.meshgrid(x, y)	
	n=50
	xg = np.linspace(x.min(),x.max(),n)
	yg = np.linspace(y.min(),y.max(),n)
	X,Y = np.meshgrid(xg,yg)
	# interpolate Z values on defined grid

	Z = griddata(np.vstack((x.flatten(),y.flatten())).T, phi.flatten().T, \
			        (X,Y), method='cubic').reshape(X.shape)
	# mask nan values, so they will not appear on plot
	Zm = np.ma.masked_where(np.isnan(Z),Z)	



	fig = plt.figure(0,figsize=(fig_width,fig_height))
	ax = fig.add_subplot(111, aspect='equal')		
	c=ax.pcolormesh(X,Y,Zm,shading='gouraud')	

	cbar= fig.colorbar(c)
	cbar.ax.tick_params(labelsize=gcafontSize)							

	c.set_clim(0, 5)

	cbar.ax.set_ylabel(r'$/phi$', fontsize=1.2*gcafontSize)
	cl = plt.getp(cbar.ax, 'ymajorticklabels')
	plt.setp(cl, fontsize=gcafontSize) 				

	ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
	ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
	#plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
	#plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)	


#	ax.set_xticks([])
#	ax.set_yticks([])
	fig.tight_layout()


	plt.savefig(fig_name)
	plt.show()
	plt.close()
	print fig_name+' saved!'	


if TEST_2D:
	fig_name = 'dy.png'
	database={}

	database['x']=np.loadtxt('mesh_x.dat')
	database['y']=np.loadtxt('mesh_y.dat')
	database['phi']=np.loadtxt('solution_1.dat')
	#phi_analytical = 
	x = database['x']
	y = database['y']
	phi = database['phi']	
	ana_x=-2*np.pi*np.cos(2*np.pi*x)*np.sin(2*np.pi*y)

	n=20
	xg = np.linspace(x.min(),x.max(),n)
	yg = np.linspace(y.min(),y.max(),n)
	X,Y = np.meshgrid(xg,yg)
	# interpolate Z values on defined grid
	Z = griddata(np.vstack((x.flatten(),y.flatten())).T, phi.flatten().T, \
			        (X,Y), method='cubic').reshape(X.shape)
	ZZ = griddata(np.vstack((x.flatten(),y.flatten())).T, ana_x.flatten().T, \
			        (X,Y), method='cubic').reshape(X.shape)	
	# mask nan values, so they will not appear on plot
	xxx = np.linspace(0,1,n)
	yyy = np.linspace(0,1,n)	
	XX,YY = np.meshgrid(x,y)	
	#ZZ=-np.cos(np.pi*XX)*np.sin(np.pi*YY)


	Zm = np.ma.masked_where(np.isnan(Z),Z)	
	ZZm = np.ma.masked_where(np.isnan(ZZ),ZZ)

	fig = plt.figure(0,figsize=(fig_width,fig_height))
	ax = fig.add_subplot(121, aspect='equal')
	ax2 = fig.add_subplot(122, aspect='equal')
	c=ax.pcolormesh(X,Y,Zm,shading='gouraud')	
	c2 = ax2.pcolormesh(X,Y,ZZm,shading='gouraud')	

	cbar= fig.colorbar(c)
	#c2bar = fig.colorbar(c2)

	cbar.ax.tick_params(labelsize=gcafontSize)							

	#c.set_clim(0, 3)

	cbar.ax.set_ylabel(r'$/phi$', fontsize=1.2*gcafontSize)
	cl = plt.getp(cbar.ax, 'ymajorticklabels')
	plt.setp(cl, fontsize=gcafontSize) 				

	ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
	ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
	#plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
	#plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)	


#	ax.set_xticks([])
#	ax.set_yticks([])
	fig.tight_layout()


	plt.savefig(fig_name)
	plt.show()
	plt.close()
	print fig_name+' saved!'	
	set_trace()

if ONE_D:

	database={}

	database['F_sol_1']=np.loadtxt('solution_500.dat')
	database['F_sol_2']=np.loadtxt('solution_1500.dat')
	database['F_sol_3']=np.loadtxt('solution_3000.dat')
	F_sol_1=database['F_sol_1']
	F_sol_2=database['F_sol_2']
	F_sol_3=database['F_sol_3']
	figname ='solution.png'
	fig = plt.figure(0,figsize=(fig_width,fig_height))#figsize=(fig_width,fig_height))
	ax = fig.add_subplot(111)	

	plt.plot(F_sol_1[:,0],F_sol_1[:,1],'--r',linewidth=lineWidth,label='FTCS')
	plt.plot(F_sol_2[:,0],F_sol_2[:,1],'--g',linewidth=lineWidth,label='FTCS')
	plt.plot(F_sol_3[:,0],F_sol_3[:,1],'--k',linewidth=lineWidth,label='FTCS')



	#handles, labels = ax.get_legend_handles_labels()
	#ax.legend(handles, labels)
	plt.ylabel('u',fontsize=gcafontSize*0.8)
	plt.xlabel('x',fontsize=gcafontSize*0.8)
	#	ax.set_title('Explicit FTCS')
	#plt.legend([aa,cc,bb], ('absolut error','ref_line : 2th order','ref_line : 1th order'),loc='best')
	plt.grid(True, which='both')		#	
	plt.tight_layout()
	plt.savefig(figname)




	plt.show()
	plt.close()
	print figname+' saved!'

	set_trace()



if TEST_MAP	:
	fig_name = 'one_map.png'
	database={}

	database['x']=np.loadtxt('mesh_x.dat')
	database['y']=np.loadtxt('mesh_y.dat')

	x = database['x']
	y = database['y']
	phi = np.ones(x.shape)*5				
	n=10
	xg = np.linspace(0,2,n)
	yg = np.linspace(0,2,n)
	X,Y = np.meshgrid(xg,yg)


	fig = plt.figure(0,figsize=(fig_width,fig_height))#figsize=(fig_width,fig_height))
	ax = fig.add_subplot(111)	

	ax.plot(x,y,'o')
		
	ax.axis('equal')
	ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
	ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
	#plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
	#plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)	


	#ax.set_xticks([])
	#ax.set_yticks([])
	fig.tight_layout()


	plt.savefig(fig_name)
	plt.show()
	plt.close()
	print fig_name+' saved!'		



if TEST_MAP2	:

	filename = "./mesh_folder/single11.msh"
	#filename = 'nine.msh'
	#filename = 'one.msh'
	#filename = 'two.msh'
	#filename  = 'one_skew.msh'
	#filename = 'quad.msh'
	mshfile_fullpath = filename#icemcfd_project_folder + filename



	part_names, xy_no, xy_fa, xy_cv, noofa, cvofa, faono, faocv, partofa = \
	        umesh_reader.read_unstructured_grid(mshfile_fullpath,node_reordering=True)
	     
	  #####################################################
	  ########## Plot Grid Labels / Connectivity ##########
	  ##################################################### 

	Plot_Node_Labels = False
	Plot_Face_Labels = False
	Plot_CV_Labels   = False
	mgplx = 0.05*np.abs(max(xy_no[:,0])-min(xy_no[:,0]))
	mgply = 0.05*np.abs(max(xy_no[:,1])-min(xy_no[:,1]))
	xlimits = [min(xy_no[:,0])-mgplx,max(xy_no[:,0])+mgplx]
	ylimits = [min(xy_no[:,1])-mgply,max(xy_no[:,1])+mgply] 

	fig = plt.figure(0,figsize=(fig_width,fig_height))
	ax = fig.add_subplot(111)
	ax.plot(xy_no[:,0],xy_no[:,1],'o',markersize=5,markerfacecolor='k') 
	node_color = 'k'
	centroid_color = 'r'  



	fig_name = 'nine_map.png'
	database={}

	database['x']=np.loadtxt('mesh_x.dat')
	database['y']=np.loadtxt('mesh_y.dat')

	x = database['x']
	y = database['y']
	phi = np.ones(x.shape)*5				
	n=10
	xg = np.linspace(0,2,n)
	yg = np.linspace(0,2,n)
	X,Y = np.meshgrid(xg,yg)


	fig = plt.figure(0,figsize=(fig_width,fig_height))#figsize=(fig_width,fig_height))
	ax = fig.add_subplot(111)	

	ax.plot(x[0:11],y[0:11],'ro')
	ax.plot(x[11:22],y[11:22],'go')
	ax.plot(x[22:33],y[22:33],'bo')
	ax.plot(x[33:44],y[33:44],'yo')			
	for inos_of_fa in noofa:
	   ax.plot(xy_no[inos_of_fa,0], xy_no[inos_of_fa,1], 'k-', linewidth = lineWidth) 		
	ax.axis('equal')
	ax.set_xlabel(r'$x$',fontsize=1.2*gcafontSize)
	ax.set_ylabel(r'$y$',fontsize=1.2*gcafontSize)
	#plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
	#plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)	


	#ax.set_xticks([])
	#ax.set_yticks([])
	fig.tight_layout()


	plt.savefig(fig_name)
	plt.show()
	plt.close()
	print fig_name+' saved!'		
