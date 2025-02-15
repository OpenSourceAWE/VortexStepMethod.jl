# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 15:01:05 2022

@author: oriol
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 00:26:08 2022

@author: oriol
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 18:09:16 2022

@author: oriol2
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
sys.path.insert(0, '../functions/')
import functions_VSM_LLT as VSM
plt.rcParams.update({'font.size': 14})
import pandas as pd

#%%  Input DATA
N = 61
#   Rectangular wing coordinates
chord = np.ones(N)*1
gamma = 20/180*np.pi
AR = 12
span = AR*max(chord)
Atot = span*max(chord)
beta = np.zeros(N)
offset = np.tan(gamma)*span/2
# beta = np.concatenate((np.linspace(-np.pi/6,0,round(N/2)),np.linspace(0,-np.pi/6,round(N/2))),axis =None)
twist = np.ones(N)*0
dist = 'lin'
coord = VSM.generate_coordinates_swept_wing(chord,offset, span, twist, beta, N,dist)
gamma = np.arctan(offset/(span/2))*180/np.pi
Atot = span*max(chord)

# Wind speed mag and direction  
Umag = 2
aoa = 18
aoa = aoa*np.pi/180
Uinf = np.array([np.cos(aoa),0,np.sin(aoa)])*Umag

#   Model and program specifics
ring_geo = '5fil'
model = 'VSM'

conv_crit = {
    'Niterations': 1500,
    'error' : 1e-5,
    'Relax_factor': 0.03
    }

Gamma0 = np.zeros(N-1)
polars_airf = pd.read_csv('./polars/NACA4415_XFOIL_Re3e6.csv', skiprows=9)
alpha_airf = np.arange(-10,30)
data_airf = np.zeros((len(polars_airf),4))
data_airf[:,0] = polars_airf['alpha']
data_airf[:,1] = polars_airf['CL']
data_airf[:,2] = polars_airf['CD']

polars_airf = np.loadtxt('./polars/NACA4415_CFD_Re3e6.csv',delimiter = ',')
data_airf = np.zeros((len(polars_airf),4))
data_airf[:,0] = polars_airf[:,0]
data_airf[:,1] = polars_airf[:,1]
#%% SOLVER  
start_time = time.time()
# Define system of vorticity
controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_general(
    coord, Uinf, N,ring_geo,model)
# Solve for Gamma

Fmag, Gamma,aero_coeffs= VSM.solve_lifting_line_system_matrix_approach_art_visc(
    ringvec, controlpoints, rings,Uinf,Gamma0,data_airf,conv_crit,model)
end_time = time.time()
print(end_time-start_time)

#%% OUTPUT Results

F_rel,F_gl,Ltot,Dtot,CL,CD,CS = VSM.output_results(Fmag,aero_coeffs,ringvec,Uinf,controlpoints,Atot)
print(CL,CD)
#%% PLOTS

fig = plt.figure( figsize = (6,5))
plt.plot(coord_L[:,1]/max(coord_L[:,1]),aero_coeffs[:,1], label = model) 
plt.title('Rectangular wing AR =' + str(AR)+r'$, \alpha = $' + str(aoa*180/np.pi))    
plt.xlabel(r'$y/s$')
plt.ylabel(r'$C_l$')
plt.legend()
plt.show()

fig = plt.figure(figsize = (6,5))
plt.plot(coord_L[:,1]/max(coord_L[:,1]),aero_coeffs[:,0]*180/np.pi, label = model) 
plt.title('Rectangular wing AR =' + str(AR)+r'$, \alpha = $' + str(aoa*180/np.pi))    
plt.xlabel(r'$y/s$')
plt.ylabel(r'$\alpha$')
plt.legend()
plt.show()

fig = plt.figure( figsize = (6,5))
plt.plot(coord_L[:,1]/max(coord_L[:,1]),aero_coeffs[:,2], label = model) 
plt.title('Rectangular wing AR =' + str(AR)+r'$, \alpha = $' + str(aoa*180/np.pi))    
plt.xlabel(r'$y/s$')
plt.ylabel(r'Viscous $C_d$')
plt.legend()
plt.show()

plt.figure(figsize = (6,4))
plt.plot(coord_L[:,1]/max(coord_L[:,1]),Gamma, label =  'No visc.') 

plt.xlabel(r'$y/s$')
plt.ylabel(r'$\Gamma$')
plt.legend()
plt.grid()

# VSM.plot_geometry(bladepanels,controlpoints,rings,F_gl,coord_L,'True')

#%% PLOT GEOMETRY

# plt_path = r'./plots/'
# plt.figure()
# ax = plt.axes(projection='3d')

# for panel in bladepanels:
#     coord = np.array([panel['p1'],panel['p2'],panel['p3'],panel['p4']])
#     ax.plot3D(coord[:,0], coord[:,1], coord[:,2], 'k',linestyle = '--')
    
# for ring in rings:
#     for filament in ring:
        
#         if filament['id'] == 'trailing_inf1' or filament['id'] == 'trailing_inf2':
#             coord = np.array([filament['x1'],filament['x1']+filament['dir']*1])
#             ax.plot3D(coord[:,0], coord[:,1], coord[:,2], 'gray')
#         else:
#             coord = np.array([filament['x1'],filament['x2']])
#             ax.plot3D(coord[:,0], coord[:,1], coord[:,2])