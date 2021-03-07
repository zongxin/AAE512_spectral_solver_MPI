import os
import sys
import numpy as np
import scipy as sp 
import scipy.sparse as scysparse
import scipy.sparse.linalg as splinalg
import pylab as plt
from pdb import set_trace
import umesh_reader
from matplotlib import rc as matplotlibrc


UD_FLIP= 0
#icemcfd_project_folder = 'mesh_folder/'

#filename = 'single11.msh'
#filename = 'nine.msh'
#filename = 'one.msh'
#filename = 'two.msh'
#filename  = 'one_skew.msh'
#filename = 'quad.msh'
#filename = '7_shape.msh'
#filename = 'fluent.msh'
#filename = 'upper.msh'
#filename = 'eye.msh'
filename  = 'final.msh'
mshfile_fullpath = filename#icemcfd_project_folder + filename



part_names, xy_no, xy_fa, xy_cv, noofa, cvofa, faono, faocv, partofa = \
        umesh_reader.read_unstructured_grid(mshfile_fullpath,node_reordering=True)
     
  #####################################################
  ########## Plot Grid Labels / Connectivity ##########
  ##################################################### 
if UD_FLIP==1:
  filename = 'down.msh'
  xy_no[:,1]=- xy_no[:,1]
  xy_fa[:,1]=- xy_fa[:,1]
  xy_cv[:,1]=- xy_cv[:,1]   
      
fig_width = 30*0.3
fig_height = 17*0.3
textFontSize   = 15
gcafontSize    = 32
lineWidth      = 2  
Plot_Node_Labels = False
Plot_Face_Labels = True
Plot_CV_Labels   = True 

# the following enables LaTeX typesetting, which will cause the plotting to take forever..  
# matplotlibrc('text.latex', preamble='\usepackage{color}')
# matplotlibrc('text',usetex=True)
# matplotlibrc('font', family='serif')  

mgplx = 0.05*np.abs(max(xy_no[:,0])-min(xy_no[:,0]))
mgply = 0.05*np.abs(max(xy_no[:,1])-min(xy_no[:,1]))
xlimits = [min(xy_no[:,0])-mgplx,max(xy_no[:,0])+mgplx]
ylimits = [min(xy_no[:,1])-mgply,max(xy_no[:,1])+mgply] 

fig = plt.figure(0,figsize=(fig_width,fig_height))
ax = fig.add_subplot(111)
ax.plot(xy_no[:,0],xy_no[:,1],'o',markersize=5,markerfacecolor='k') 
node_color = 'k'
centroid_color = 'r'  

for inos_of_fa in noofa:
   ax.plot(xy_no[inos_of_fa,0], xy_no[inos_of_fa,1], 'k-', linewidth = lineWidth) 

if Plot_Face_Labels:
  nfa = xy_fa.shape[0] # number of faces
  faces_indexes = range(0,nfa)
  for x_fa,y_fa,ifa in zip(xy_fa[:,0],xy_fa[:,1],faces_indexes):
    ax.text(x_fa,y_fa,repr(ifa),transform=ax.transData,color='k',
        verticalalignment='center',horizontalalignment='center',fontsize=textFontSize ) 

if Plot_Node_Labels:
  nno = xy_no.shape[0] # number of nodes
  node_indexes = range(0,nno)
  for xn,yn,ino in zip(xy_no[:,0],xy_no[:,1],node_indexes):
    ax.text(xn,yn,repr(ino),transform=ax.transData,color='r',
        verticalalignment='top',horizontalalignment='left',fontsize=textFontSize )  

if Plot_CV_Labels:
  ncv = xy_cv.shape[0]  # number of control volumes
  cv_indexes = range(0,ncv)
  for xcv,ycv,icv in zip(xy_cv[:,0],xy_cv[:,1],cv_indexes):
    ax.text(xcv,ycv,repr(icv),transform=ax.transData,color='b',
        verticalalignment='top',horizontalalignment='left',fontsize=textFontSize )  

ax.axis('equal')
ax.set_xlim(xlimits)
ax.set_ylim(ylimits)
ax.set_xlabel(r'$x$',fontsize=1.5*gcafontSize)
ax.set_ylabel(r'$y$',fontsize=1.5*gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)

fig_name = filename.split('.')[0]+'55555.png'
figure_path = '../report/figures/'
fig_fullpath = fig_name#figure_path + fig_name
plt.savefig(fig_fullpath)
plt.show()



