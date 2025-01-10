#%% Import libraries
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

#%% Read data
# Read documentation values

u_ghia = pd.read_csv('utils/Prof_u_Ghia_Re_100_a_10000.dat', sep='\s+')
y_ghia = u_ghia['y']

v_ghia = pd.read_csv('utils/Prof_v_Ghia_Re_100_a_10000.dat', sep='\s+')
x_ghia = v_ghia['x']
print(v_ghia)

# Read simulation values
u_simu = pd.read_csv('../src/u.dat', sep='\s+')
y_simu = u_simu['y']
v_simu = pd.read_csv('../src/v.dat', sep='\s+')
x_simu = v_simu['x']

#%% Analyse for Re=100

u_ghia_Re100 = u_ghia['Re=100']
u_simu_Re100 = u_simu['u']
v_ghia_Re100 = v_ghia['Re=100']
v_simu_Re100 = v_simu['v']


#%% Plot u velocity
plt.figure()
plt.plot(y_ghia, u_ghia_Re100, label='Ghia et al.')
plt.plot(y_simu, u_simu_Re100, label='Simulation')
plt.xlabel('y')
plt.ylabel('u')
plt.legend()
plt.title('u velocity for Re=100')
plt.show()


#%% Plot v velocity
plt.figure()
plt.plot(x_ghia, v_ghia_Re100, label='Ghia et al.')
plt.plot(x_simu, v_simu_Re100, label='Simulation')
plt.xlabel('x')
plt.ylabel('v')
plt.legend()
plt.title('v velocity for Re=100')
plt.show()
