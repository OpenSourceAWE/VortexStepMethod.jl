# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 11:49:33 2022

@author: oriol
"""


import numpy as np
import matplotlib.pyplot as plt

data_airf = np.loadtxt(r'./polars/clarky_maneia.csv',delimiter = ',')

plt_path = './plots/clarky/'


plt.figure(figsize=(6, 5))
plt.plot(data_airf[0:13,0],data_airf[0:13,1],marker= '.')
plt.plot(data_airf[12::,0],data_airf[12::,1],marker= '.')
plt.legend(['CFD','Polynomial approx.'])
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$C_l  $')
plt.grid()
plt.savefig(plt_path+'Clarky_Cl_alpha.png',dpi = 150)


plt.figure(figsize=(6, 5))
plt.plot(data_airf[0:13,0],data_airf[0:13,2],marker= '.')
plt.plot(data_airf[12::,0],data_airf[12::,2],marker= '.')
plt.legend(['CFD','Polynomial approx.'])
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$C_d  $')
plt.grid()
plt.savefig(plt_path+'Clarky_Cd_alpha.png',dpi = 150)

plt.figure(figsize=(6, 5))
plt.plot(data_airf[0:13,0],data_airf[0:13,1]/data_airf[0:13,2],marker= '.')
plt.plot(data_airf[12::,0],data_airf[12::,1]/data_airf[12::,2],marker= '.')
plt.legend(['CFD','Polynomial approx.'])
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$C_l/C_d  $')
plt.grid()
plt.savefig(plt_path+'Clarky_ClCd_alpha.png',dpi = 150)

plt.figure(figsize=(6, 5))
plt.plot(data_airf[0:13,2],data_airf[0:13,1],marker= '.')
plt.plot(data_airf[12::,2],data_airf[12::,1],marker= '.')
plt.legend(['CFD','Polynomial approx.'])
plt.xlabel(r'$C_d$')
plt.ylabel(r'$C_l  $')
plt.grid()
plt.savefig(plt_path+'Clarky_Cd_Cl.png',dpi = 150)


#%% CAMBER LINE
airf_path = r'./polars/clarky.dat'

a_file = open(airf_path)

lines = a_file.readlines()[2:]
lines.pop(61)

data = np.zeros((len(lines),2))
i = 0
for line in lines:
    data[i,:] = [float(i) for i in line.split(' ')[1:]]   
    i +=1
    
x_upp = data[0:61,:]
x_low = data[61::,:]

x_camber = np.array([(x_upp[i,:]+x_low[i,:])/2 for i in range(len(x_upp))])

aoa_0L = np.arctan((x_camber[45,1]-x_camber[44,1])/(x_camber[45,0]-x_camber[44,0]))*180/np.pi

plt.figure(figsize=(6, 2))
plt.plot(data[:,0],data[:,1],marker= '.')
plt.plot(x_camber[:,0],x_camber[:,1])
plt.xlabel(r'$x$')
plt.ylabel(r'$y  $')
plt.grid()
plt.savefig(plt_path+'Clarky_coord_camber.png',dpi = 150)