# ################################################################
#  
#    Active Particles on Curved Spaces (APCS)
#    
#    Author: Silke Henkes
#   
#    ICSMB, Department of Physics
#    University of Aberdeen
#    
#    (c) 2013, 2014
#    
#    This program cannot be used, copied, or modified without
#    explicit permission of the author.
#  
# ################################################################

#! /usr/bin/python

import sys, os, glob
import cPickle as pickle
import numpy as np
import scipy as sp
from scipy.io import savemat
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
#from matplotlib import rc
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

# setting global parameters
#matplotlib.rcParams['text.usetex'] = 'true'
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['font.size']=20
matplotlib.rcParams['legend.fontsize']=14


cdict = {'red':   [(0.0,  0.75, 0.75),
				   (0.3,  1.0, 1.0),
                   (0.5,  0.4, 0.0),
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.0),
				   (0.25,  0.0, 0.5),
                   (0.5, 1.0, 1.0),
                   (0.75,  0.5, 0.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (0.7, 1.0, 1.0),
                   (1.0,  0.25, 0.25)]}
                   

	

# This is the structured data file hierarchy. Replace as appropriate (do not go the Yaouen way and fully automatize ...)
basefolder = '/home/silke/Documents/CurrentProjects/Rastko/Runs/'
outfolder= '/home/silke/Documents/CurrentProjects/Rastko/analysis/'
JList=['10', '1', '0.1', '0.01']
vList=['0.005','0.01','0.02','0.05','0.1','0.2','0.5','1']
#JList=['10']
#vList=['1']
nbin=180

testmap=LinearSegmentedColormap('test',cdict,N=len(vList))
testmap2=LinearSegmentedColormap('test',cdict,N=len(JList))


# Profiles
# Set column to plot
usecolumn=1

profList=[r'$\theta$',r'$\rho$',r'$\sqrt{\langle v^2 \rangle}/v_0$','energy','pressure',r'$\alpha$',r'$\alpha_v$']
profName=['theta','rho','vrms','energy','pressure','alpha','alpha_v']

for j in range(len(JList)):
	plt.figure(figsize=(10,7),linewidth=2.0)
	for i in range(len(vList)):
		print vList[i],JList[j]
		ax=plt.gca()
		outfile=outfolder+'data/profiles_v0' + vList[i] + '_j' + JList[j] + '.dat'
		outfile2=outfolder + 'data/axis_v0' + vList[i] + '_j' + JList[j] + '.dat'
		# header='theta rho vel energy pressure alpha alpha_v'
		profiles=sp.loadtxt(outfile, unpack=True)[:,:] 
		isdata=[index for index,value in enumerate(profiles[1,:]) if (value >0)]
		if usecolumn==2:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata]/float(vList[i]),color=testmap(i), linestyle='solid',label=vList[i])
		else:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata],color=testmap(i), linestyle='solid',label=vList[i])
	if usecolumn>=5:
		plt.plot(profiles[0,isdata],profiles[0,isdata],'k-')
		plt.ylim(-0.5,0.5)
		plt.xlim(-np.pi/2,np.pi/2)
	plt.xlabel(profList[0]) 
	plt.ylabel(profList[usecolumn]) 
	plt.legend(loc=2,ncol=2)
	plt.title('Interaction strength ' + r'$J=$' + JList[j])
	
	filename=outfolder + 'pics/profile_' + profName[usecolumn] + '_J' + JList[j] +'.pdf'
	plt.savefig(filename)

for i in range(len(vList)):	
	plt.figure(figsize=(10,7),linewidth=2.0)
	for j in range(len(JList)):
		print vList[i],JList[j]
		ax=plt.gca()
		outfile=outfolder+'data/profiles_v0' + vList[i] + '_j' + JList[j] + '.dat'
		outfile2=outfolder + 'data/axis_v0' + vList[i] + '_j' + JList[j] + '.dat'
		# header='theta rho vel energy pressure alpha alpha_v'
		profiles=sp.loadtxt(outfile, unpack=True)[:,:] 
		isdata=[index for index,value in enumerate(profiles[1,:]) if (value >0)]
		if usecolumn==2:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata]/float(vList[i]),color=testmap2(j), linestyle='solid',label=JList[j])
		else:
			plt.plot(profiles[0,isdata],profiles[usecolumn,isdata],color=testmap2(j), linestyle='solid',label=JList[j])
	if usecolumn>=5:
		plt.plot(profiles[0,isdata],profiles[0,isdata],'k-')
		plt.ylim(-0.5,0.5)
		plt.xlim(-np.pi/2,np.pi/2)
	plt.xlabel(profList[0]) 
	plt.ylabel(profList[usecolumn]) 
	plt.legend(loc=2,ncol=2)
	plt.title('Velocity ' + r'$v_0=$' + vList[i])
	
	filename=outfolder + 'pics/profile_' + profName[usecolumn] + '_v' + vList[i] +'.pdf'
	plt.savefig(filename)
	
#for i in range(len(vList)):	
	#fig=plt.figure(figsize=(10,7),linewidth=2.0)
	#ax = fig.add_subplot(111, projection='3d')
	#for j in range(len(JList)):
		#print vList[i],JList[j]
		#ax=plt.gca()
		#outfile=outfolder+'data/profiles_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#outfile2=outfolder + 'data/axis_v0' + vList[i] + '_j' + JList[j] + '.dat'
		#axis=sp.loadtxt(outfile2, unpack=True)[:,:] 
		#ax.scatter(axis[0,:], axis[1,:], axis[2,:], zdir='z',color=testmap2(j), linestyle='solid',label=JList[j])
	#plt.xlabel('x') 
	#plt.ylabel('y') 
	##plt.xlim(-1,1)
	##plt.ylim(-1,1)
	##plt.zlim(-1,1)
	##ax.zlabel('z') 
	#ax.legend()
	#plt.legend()
	#plt.title('Velocity ' + r'$v_0=$' + vList[i])
	
	#filename=outfolder + 'pics/axis_' + profName[usecolumn] + '_v' + vList[i] +'.pdf'
	#plt.savefig(filename)
	


	
	
plt.show()	
		
		
		
		
		
		