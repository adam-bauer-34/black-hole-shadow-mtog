"""
Adam M. Bauer

Master file for finding the intensity profile of a cloud of gas 
around a Schwazschild black hole, with geodesic ray tracing. 

To switch between the radial infall or stationary model, see getGeodesicEvolutionPlus.py

Complementary files: GeodesicREZZ.py (general ray class), getGeodesicEvolutionREZZ.py (evolution equations)
getMiscFuncsREZZ.py (misc functions), getGeodesicEvolutionPlusREZZ.py (backwards evolution equations)

Figure files: IntensityPlot.py, eht_figure.py

To run: python MasterRadin.py

Copyright (C) 2021, Adam M. Bauer
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np 
import sympy as sy
import matplotlib.pyplot as plt 
from scipy.integrate import ode, solve_ivp
from matplotlib import font_manager
from matplotlib import ticker
import matplotlib as mpl
from src.getMiscFuncsREZZ import getIVsAndNames, getExpressions
from src.GeodesicREZZ import Geodesic
import random 

import matplotlib.colors as clrs 
import matplotlib 

cms = matplotlib.cm 

powers = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0] # choose range of powers in the accretion model, i.e., j_nu ~ r^(-powers)
zetas = [-0.20, -0.15, -0.10, -0.05, 0.05, 0.10, 0.15, 0.20]

accword = "radin" # options: "radin" "stat". chooses accretion model. 
    # "stat" is a stationary emitting gas model
    # "radin" is a radially free falling gas model (bondi accretion at zero temperature)

N = 2000 # total number of geodesics to ray trace
#b0 = 0 # initial impact parameter 
#Deltab = 0.003 # how much to increase impact paramater for each ray (i.e., b[i+1] = b[i] + Deltab)

for j in range(0, len(powers)):
    power = float(powers[j]) 

    for k in range(0, len(zetas)):
        zeta = float(zetas[k])
        
        if zeta == 0:
            theory = "GR"
        
        else: 
            theory = "REZZ"

        equation_list = getExpressions(theory, zeta) # fetch sympy expressions for various equations, imported from Mathematica

        Namearray, IVarray = getIVsAndNames(N, zeta) # make an array of names and initial values for geodesics

        print("we're on power " + str(power) + " and zeta equals " + str(zeta)) # tells you where you're at in the simulation

        for i in range(0, N):
            #print("We are on geodesic number " + str(i))
            Namearray[i] = Geodesic(IVarray[i], equation_list, zeta, accword, power, 0, 0) # make class instance for geodesic (i+1)
            #Namearray[i].getXArray() # get X array from forwards integration of geo(i+1)
            #Namearray[i].getYArray() # get Y array from forwards integration of geo(i+1)
            #Namearray[i].getEnergyArray() # get energy array from forwards integration of geo(i+1)
            #Namearray[i].getAngularMomentumArray() # get angular momentum array from forwards integration of geo(i+1)
            Namearray[i].getXArrayINT() # get X array from backwards integrated geo(i+1)
            Namearray[i].getYArrayINT() # get Y array from backwards integrated geo(i+1)
            Namearray[i].getEnergyArrayINT() # get energy array from forwards integration of geo(i+1)
            Namearray[i].getAngularMomentumArrayINT() # get angular momentum array from forwards integration of geo(i+1)

        ut_cam_fcn = equation_list[16]
        ut_cam = ut_cam_fcn(1000.) 
        freq_cam_inf = - IVarray[0][5] * ut_cam # = k_t u_cam^t, this cubed multiplied by the final intensity gives the intensity in the camera's frame! 
        
        #k = len(Namearray[0].intensity.t)
        #print(Namearray[0].intensity.y[8,k-1]*freq_cam_inf**3) # print central 

        # WRITE INTENSITY DATA TO CSV

        int_data = np.zeros((N,3))

        timesteps = []

        for i in range(0, N):
            tmp_array = []
            for k in range(0, len(Namearray[i].intensity.t)-1):
                tmp = Namearray[i].intensity.t[k+1] - Namearray[i].intensity.t[k]
                tmp_array.append(tmp)
            timesteps.append(tmp_array)

        for i in range(0,N):
            k = len(Namearray[i].intensity.t)
            int_data[i,0] = Namearray[i].b
            int_data[i,1] = Namearray[i].intensity.y[8,k-1]*freq_cam_inf**3 # invariant intensity! 
            int_data[i,2] = min(Namearray[i].intensity.y[5,:]) # minimum radius value ==> radius of closest approach

        PATH = "./data/" + accword + "/" + theory + "/" # ie "./data/accword/GR/"

        if theory == "GR":
            filename = "iData_" + accword + "_" + theory + "_p" + str(power) + "_N" + str(N) + ".csv"
        
        else: 
            filename = "iData_" + accword + "_" + theory + "_z" + str(zeta) + "_p" + str(power) + "_N" + str(N) + ".csv"

        np.savetxt(PATH+filename, int_data, delimiter=",")

"""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# INTENSITY PLOT - plots intensity curve as func of impact parameter

fig, ax = plt.subplots(1)

for i in range(0, N):
    k = len(Namearray[i].intensity.t)
    ax.plot(Namearray[i].b, Namearray[i].intensity.y[8, k-1]*freq_cam_inf**3, 'b.') # plot each rays integrated intensity

#ax.plot(Namearray[0].intensity.y[5,:], Namearray[0].intensity.y[8,:])

ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
#ax.set_xlim((0, 8))
#ax.set_ylim((-0.0005, 0.004))
ax.set_xlabel(r"$b \ (GM/c^{2})$", fontname='serif', fontstyle='normal', fontsize=20, labelpad=0)
ax.set_ylabel(r"$I \ (c^{2})$", fontname='serif', fontstyle='normal', fontsize=20, labelpad=-5)

#plt.savefig('./figures/intensity_INFALL.pdf')

# END INTENSITY PLOT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# GEODESIC PLOT - plots path of geodesics

fig, ax = plt.subplots(1)

circle1 = plt.Circle((0, 0), 2, color='k')
circle2 = plt.Circle((0, 0), 3, color='b', fill=False)

for i in range(0,N):
    label = r"$b^{2} = $" + str(Namearray[i].b)
    ax.plot(Namearray[i].XarrayINT, Namearray[i].YarrayINT, label=label, color=Namearray[i].color) # plot the position of each ray with a random colour
    # GET RID OF INT FOR FORWARDS INTEGRATED GEODESICS

#ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
ax.set_xlim((-30, 30))
ax.set_ylim((-30, 30))
ax.set_xlabel(r"$x \ (GM/c^{2})$", fontname='serif', fontstyle='normal', fontsize=20, labelpad=0)
ax.set_ylabel(r"$y \ (GM/c^{2})$", fontname='serif', fontstyle='normal', fontsize=20, labelpad=-5)
#ax.legend()

ax.add_artist(circle1) # add black hole of r =2
ax.add_artist(circle2) # add photon orbit at r = 3, for reference

ax.set_aspect('equal')

#plt.savefig('./figures/Geodesics_sgb_0dot2.pdf')

# END GEODESIC PLOT

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# START ENERGY/ANGULAR MOMENTUM PLOT - check plot for energy and angular momentum conservation

fig, ax = plt.subplots(2, sharex=True)

for i in range(0,N):
    ax[0].plot(Namearray[i].intensity.t, Namearray[i].EarrayINT, label = r"$b^{2} = $" + str(Namearray[i].b), color=Namearray[i].color, zorder=4) # plot each energy 
    ax[1].plot(Namearray[i].intensity.t, Namearray[i].LarrayINT, color=Namearray[i].color, zorder=4) # plot each angular momentum 
    # CHANGE intensity -> ray AND GET RID OF 'INT' ON Earray TO SEE FORWARDS INTEGRATED GEODESIC

ax[0].axhline(y=0, color='k')
ax[0].axvline(x=0, color='k')
#ax[0].set_xlim((0, 10))
ax[0].set_ylabel(r"$E$", fontname='serif', fontstyle='normal', fontsize=20, labelpad=-5)
#ax[0].legend()

ax[1].axhline(y=0, color='k')
ax[1].axvline(x=0, color='k')
#ax[1].set_xlim((0, 10))
ax[1].set_xlabel(r"$\lambda$", fontname='serif', fontstyle='normal', fontsize=20, labelpad=0)
ax[1].set_ylabel(r"$L$", fontname='serif', fontstyle='normal', fontsize=20, labelpad=-5)

#plt.savefig('./figures/ELvsT_sgb_0dot2.pdf')

# END ENERGY/ANGULAR MOMENTUM PLOT

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


plt.show()
"""
