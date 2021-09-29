# Adam M. Bauer
# This code makes a figure that plots our \vartheta
# as a function of emissivity power, \alpha, in GR and 
# \vartheta as a function of deformation parameter in sGB
# and RZ for a fixed value of \alpha using the Bondi 
# accretion model.

"""
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
import glob
import matplotlib.pyplot as plt 
import matplotlib as mpl
from scipy.stats import linregress
from scipy.integrate import simps 

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

cms = mpl.cm 

def getIntegration(integrand, xvals):
    tmp = 0

    for i in range(1, len(xvals)):
        tmp += 0.5 * (xvals[i] - xvals[i-1]) * (integrand[i] + integrand[i-1])
    
    return tmp

def getValue(func, integrand, xvals, tol, args): # root finder
    
    fullintegral = args[0]
    lowx = min(xvals)
    highx = max(xvals)
    mid = 0.5 * (lowx + highx)
    tracker = 0

    while 1:
        index = max(np.where(xvals <= mid)[0])

        test_xvals = xvals[0:index]
        trunc_integrand = integrand[0:index]

        test_val = func(trunc_integrand, test_xvals) * fullintegral**(-1) - 0.5
        tracker += 1

        #print(mid, test_val)

        if abs(test_val) < tol:
            return mid 
        elif test_val < 0:
            lowx = mid
            mid = 0.5 * (mid + highx)
        else: 
            highx = mid 
            mid = 0.5 * (mid + lowx)

        if tracker > 30:
            print('help i am stuck!')

def getCriticalB(zeta):
    return np.sqrt(27) - 4397*zeta*(2430 * np.sqrt(3) )**(-1) + zeta**2 * (113374080 * np.sqrt(3))**(-1) *336780431 + zeta**3 * 114832399336847 * (np.sqrt(3) * 37192366944000)**(-1) + zeta**4 * 125183193573305833463 * (np.sqrt(3) * 52057412164177920000)**(-1)

def getCrtiticalBREZZ(a1):
    return np.sqrt(27) - 4 * a1 * (3 * np.sqrt(3))**(-1) + 64 * a1**3 * (729 * np.sqrt(3))**(-1) + 512 * a1**4 * (59049 * np.sqrt(3))**(-1) - 11264 * a1**5 * (531441 * np.sqrt(3))**(-1) 

PATH = "../../data/bondi/"

GRfilenames_r10_g53 = "GR/iData_bondi_GR_p*_rcrit10.0_g53_N2000.csv"
GRfilenames_r10_g139 = "GR/iData_bondi_GR_p*_rcrit10.0_g139_N2000.csv"
GRfilenames_r30_g53 = "GR/iData_bondi_GR_p*_rcrit30.0_g53_N2000.csv"
GRfilenames_r30_g139 = "GR/iData_bondi_GR_p*_rcrit30.0_g139_N2000.csv"

SGBfilenames_p5_r10_g53 = "SGB/iData_bondi_SGB_z*_p" + str(5.0) + "_rcrit10.0_g53_N2000.csv"
SGBfilenames_p5_r10_g139 = "SGB/iData_bondi_SGB_z*_p" + str(5.0) + "_rcrit10.0_g139_N2000.csv"
SGBfilenames_p5_r30_g53 = "SGB/iData_bondi_SGB_z*_p" + str(5.0) + "_rcrit30.0_g53_N2000.csv"
SGBfilenames_p5_r30_g139 = "SGB/iData_bondi_SGB_z*_p" + str(5.0) + "_rcrit30.0_g139_N2000.csv"

SGBfilenames_p6_r10_g53 = "SGB/iData_bondi_SGB_z*_p" + str(6.0) + "_rcrit10.0_g53_N2000.csv"
SGBfilenames_p6_r10_g139 = "SGB/iData_bondi_SGB_z*_p" + str(6.0) + "_rcrit10.0_g139_N2000.csv"
SGBfilenames_p6_r30_g53 = "SGB/iData_bondi_SGB_z*_p" + str(6.0) + "_rcrit30.0_g53_N2000.csv"
SGBfilenames_p6_r30_g139 = "SGB/iData_bondi_SGB_z*_p" + str(6.0) + "_rcrit30.0_g139_N2000.csv"

SGBfilenames_p7_r10_g53 = "SGB/iData_bondi_SGB_z*_p" + str(7.0) + "_rcrit10.0_g53_N2000.csv"
SGBfilenames_p7_r10_g139 = "SGB/iData_bondi_SGB_z*_p" + str(7.0) + "_rcrit10.0_g139_N2000.csv"
SGBfilenames_p7_r30_g53 = "SGB/iData_bondi_SGB_z*_p" + str(7.0) + "_rcrit30.0_g53_N2000.csv"
SGBfilenames_p7_r30_g139 = "SGB/iData_bondi_SGB_z*_p" + str(7.0) + "_rcrit30.0_g139_N2000.csv"

SGBfilenames_p8_r10_g53 = "SGB/iData_bondi_SGB_z*_p" + str(8.0) + "_rcrit10.0_g53_N2000.csv"
SGBfilenames_p8_r10_g139 = "SGB/iData_bondi_SGB_z*_p" + str(8.0) + "_rcrit10.0_g139_N2000.csv"
SGBfilenames_p8_r30_g53 = "SGB/iData_bondi_SGB_z*_p" + str(8.0) + "_rcrit30.0_g53_N2000.csv"
SGBfilenames_p8_r30_g139 = "SGB/iData_bondi_SGB_z*_p" + str(8.0) + "_rcrit30.0_g139_N2000.csv"

pREZZfilenames6_r10_g53 = "REZZ/iData_bondi_REZZ_z0*_p" + str(6.0) + "_rcrit10.0_g53_N2000.csv"
pREZZfilenames6_r10_g139 = "REZZ/iData_bondi_REZZ_z0*_p" + str(6.0) + "_rcrit10.0_g139_N2000.csv"
pREZZfilenames6_r30_g53 = "REZZ/iData_bondi_REZZ_z0*_p" + str(6.0) + "_rcrit30.0_g53_N2000.csv"
pREZZfilenames6_r30_g139 = "REZZ/iData_bondi_REZZ_z0*_p" + str(6.0) + "_rcrit30.0_g139_N2000.csv"

nREZZfilenames6_r10_g53 = "REZZ/iData_bondi_REZZ_zneg*_p" + str(6.0) + "_rcrit10.0_g53_N2000.csv"
nREZZfilenames6_r10_g139 = "REZZ/iData_bondi_REZZ_zneg*_p" + str(6.0) + "_rcrit10.0_g139_N2000.csv"
nREZZfilenames6_r30_g53 = "REZZ/iData_bondi_REZZ_zneg*_p" + str(6.0) + "_rcrit30.0_g53_N2000.csv"
nREZZfilenames6_r30_g139 = "REZZ/iData_bondi_REZZ_zneg*_p" + str(6.0) + "_rcrit30.0_g139_N2000.csv"

pREZZfilenames8_r10_g53 = "REZZ/iData_bondi_REZZ_z0*_p" + str(8.0) + "_rcrit10.0_g53_N2000.csv"
pREZZfilenames8_r10_g139 = "REZZ/iData_bondi_REZZ_z0*_p" + str(8.0) + "_rcrit10.0_g139_N2000.csv"
pREZZfilenames8_r30_g53 = "REZZ/iData_bondi_REZZ_z0*_p" + str(8.0) + "_rcrit30.0_g53_N2000.csv"
pREZZfilenames8_r30_g139 = "REZZ/iData_bondi_REZZ_z0*_p" + str(8.0) + "_rcrit30.0_g139_N2000.csv"

nREZZfilenames8_r10_g53 = "REZZ/iData_bondi_REZZ_zneg*_p" + str(8.0) + "_rcrit10.0_g53_N2000.csv"
nREZZfilenames8_r10_g139 = "REZZ/iData_bondi_REZZ_zneg*_p" + str(8.0) + "_rcrit10.0_g139_N2000.csv"
nREZZfilenames8_r30_g53 = "REZZ/iData_bondi_REZZ_zneg*_p" + str(8.0) + "_rcrit30.0_g53_N2000.csv"
nREZZfilenames8_r30_g139 = "REZZ/iData_bondi_REZZ_zneg*_p" + str(8.0) + "_rcrit30.0_g139_N2000.csv"

zetas = [0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175]
powers = [5.0, 5.5, 6.0, 6.5, 7.0, 8.0]
a1s = [-0.2, -0.15, -0.1, -0.05, 0.05, 0.1, 0.15, 0.2]

GRfiles_r10_g53 = glob.glob(PATH+GRfilenames_r10_g53)
GRfiles_r10_g139 = glob.glob(PATH+GRfilenames_r10_g139)
GRfiles_r30_g53 = glob.glob(PATH+GRfilenames_r30_g53)
GRfiles_r30_g139 = glob.glob(PATH+GRfilenames_r30_g139)

SGBfiles_p5_r10_g53  = glob.glob(PATH+SGBfilenames_p5_r10_g53)
SGBfiles_p5_r10_g139 = glob.glob(PATH+SGBfilenames_p5_r10_g139)
SGBfiles_p5_r30_g53 = glob.glob(PATH+SGBfilenames_p5_r30_g53)
SGBfiles_p5_r30_g139  = glob.glob(PATH+SGBfilenames_p5_r30_g139)

SGBfiles_p6_r10_g53  = glob.glob(PATH+SGBfilenames_p6_r10_g53)
SGBfiles_p6_r10_g139 = glob.glob(PATH+SGBfilenames_p6_r10_g139)
SGBfiles_p6_r30_g53 = glob.glob(PATH+SGBfilenames_p6_r30_g53)
SGBfiles_p6_r30_g139  = glob.glob(PATH+SGBfilenames_p6_r30_g139)

SGBfiles_p7_r10_g53  = glob.glob(PATH+SGBfilenames_p7_r10_g53)
SGBfiles_p7_r10_g139 = glob.glob(PATH+SGBfilenames_p7_r10_g139)
SGBfiles_p7_r30_g53 = glob.glob(PATH+SGBfilenames_p7_r30_g53)
SGBfiles_p7_r30_g139  = glob.glob(PATH+SGBfilenames_p7_r30_g139)

SGBfiles_p8_r10_g53  = glob.glob(PATH+SGBfilenames_p8_r10_g53)
SGBfiles_p8_r10_g139 = glob.glob(PATH+SGBfilenames_p8_r10_g139)
SGBfiles_p8_r30_g53 = glob.glob(PATH+SGBfilenames_p8_r30_g53)
SGBfiles_p8_r30_g139  = glob.glob(PATH+SGBfilenames_p8_r30_g139)

pREZZfiles6_r10_g53 = glob.glob(PATH+pREZZfilenames6_r10_g53)
pREZZfiles6_r10_g139 = glob.glob(PATH+pREZZfilenames6_r10_g139)
pREZZfiles6_r30_g53 = glob.glob(PATH+pREZZfilenames6_r30_g53)
pREZZfiles6_r30_g139 = glob.glob(PATH+pREZZfilenames6_r30_g139)

nREZZfiles6_r10_g53 = glob.glob(PATH+nREZZfilenames6_r10_g53)
nREZZfiles6_r10_g139 = glob.glob(PATH+nREZZfilenames6_r10_g139)
nREZZfiles6_r30_g53 = glob.glob(PATH+nREZZfilenames6_r30_g53)
nREZZfiles6_r30_g139 = glob.glob(PATH+nREZZfilenames6_r30_g139)

pREZZfiles8_r10_g53 = glob.glob(PATH+pREZZfilenames8_r10_g53)
pREZZfiles8_r10_g139 = glob.glob(PATH+pREZZfilenames8_r10_g139)
pREZZfiles8_r30_g53 = glob.glob(PATH+pREZZfilenames8_r30_g53)
pREZZfiles8_r30_g139 = glob.glob(PATH+pREZZfilenames8_r30_g139)

nREZZfiles8_r10_g53 = glob.glob(PATH+nREZZfilenames8_r10_g53)
nREZZfiles8_r10_g139 = glob.glob(PATH+nREZZfilenames8_r10_g139)
nREZZfiles8_r30_g53 = glob.glob(PATH+nREZZfilenames8_r30_g53)
nREZZfiles8_r30_g139 = glob.glob(PATH+nREZZfilenames8_r30_g139)

GRfiles_r10_g53.sort() # put low p vals first
GRfiles_r10_g139.sort() # put low p vals first
GRfiles_r30_g53.sort() # put low p vals first
GRfiles_r30_g139.sort() # put low p vals first

SGBfiles_p5_r10_g53.sort()
SGBfiles_p5_r10_g139.sort()
SGBfiles_p5_r30_g53.sort()
SGBfiles_p5_r30_g139.sort()

SGBfiles_p6_r10_g53.sort()
SGBfiles_p6_r10_g139.sort()
SGBfiles_p6_r30_g53.sort()
SGBfiles_p6_r30_g139.sort()

SGBfiles_p7_r10_g53.sort()
SGBfiles_p7_r10_g139.sort()
SGBfiles_p7_r30_g53.sort()
SGBfiles_p7_r30_g139.sort()

SGBfiles_p8_r10_g53.sort()
SGBfiles_p8_r10_g139.sort()
SGBfiles_p8_r30_g53.sort()
SGBfiles_p8_r30_g139.sort()

pREZZfiles6_r10_g53.sort()
pREZZfiles6_r10_g139.sort()
pREZZfiles6_r30_g53.sort()
pREZZfiles6_r30_g139.sort()

nREZZfiles6_r10_g53.sort()
nREZZfiles6_r10_g139.sort()
nREZZfiles6_r30_g53.sort()
nREZZfiles6_r30_g139.sort()

nREZZfiles6_r10_g53.reverse()
nREZZfiles6_r10_g139.reverse()
nREZZfiles6_r30_g53.reverse()
nREZZfiles6_r30_g139.reverse()

pREZZfiles8_r10_g53.sort()
pREZZfiles8_r10_g139.sort()
pREZZfiles8_r30_g53.sort()
pREZZfiles8_r30_g139.sort()

nREZZfiles8_r10_g53.sort()
nREZZfiles8_r10_g139.sort()
nREZZfiles8_r30_g53.sort()
nREZZfiles8_r30_g139.sort()

nREZZfiles8_r10_g53.reverse()
nREZZfiles8_r10_g139.reverse()
nREZZfiles8_r30_g53.reverse()
nREZZfiles8_r30_g139.reverse()

Nfiles_GR = len(GRfiles_r10_g53)
GR_r10_g53_data_b = np.zeros((Nfiles_GR, 2000))
GR_r10_g53_data_I = np.zeros((Nfiles_GR, 2000))
GR_r10_g139_data_b = np.zeros((Nfiles_GR, 2000))
GR_r10_g139_data_I = np.zeros((Nfiles_GR, 2000))
GR_r30_g53_data_b = np.zeros((Nfiles_GR, 2000))
GR_r30_g53_data_I = np.zeros((Nfiles_GR, 2000))
GR_r30_g139_data_b = np.zeros((Nfiles_GR, 2000))
GR_r30_g139_data_I = np.zeros((Nfiles_GR, 2000))

for i in range(0, Nfiles_GR):
    GR_r10_g53_data_b[i,:] = np.genfromtxt(GRfiles_r10_g53[i], delimiter=",", usecols=0)
    GR_r10_g53_data_I[i,:] = np.genfromtxt(GRfiles_r10_g53[i], delimiter=",", usecols=1)

    GR_r10_g139_data_b[i,:] = np.genfromtxt(GRfiles_r10_g139[i], delimiter=",", usecols=0)
    GR_r10_g139_data_I[i,:] = np.genfromtxt(GRfiles_r10_g139[i], delimiter=",", usecols=1)

    GR_r30_g53_data_b[i,:] = np.genfromtxt(GRfiles_r30_g53[i], delimiter=",", usecols=0)
    GR_r30_g53_data_I[i,:] = np.genfromtxt(GRfiles_r30_g53[i], delimiter=",", usecols=1)

    GR_r30_g139_data_b[i,:] = np.genfromtxt(GRfiles_r30_g139[i], delimiter=",", usecols=0)
    GR_r30_g139_data_I[i,:] = np.genfromtxt(GRfiles_r30_g139[i], delimiter=",", usecols=1)

Ninterp = 5 * 10**4

all_b = np.linspace(0,15,Ninterp)
GR_r10_g53_interp_I = np.zeros((Nfiles_GR, Ninterp))
GR_r10_g139_interp_I = np.zeros((Nfiles_GR, Ninterp))
GR_r30_g53_interp_I = np.zeros((Nfiles_GR, Ninterp))
GR_r30_g139_interp_I = np.zeros((Nfiles_GR, Ninterp))

for i in range(0,Nfiles_GR):
    GR_r10_g53_interp_I[i,:] = np.interp(all_b, GR_r10_g53_data_b[i,:], GR_r10_g53_data_I[i,:])
    GR_r10_g139_interp_I[i,:] = np.interp(all_b, GR_r10_g139_data_b[i,:], GR_r10_g139_data_I[i,:])
    GR_r30_g53_interp_I[i,:] = np.interp(all_b, GR_r30_g53_data_b[i,:], GR_r30_g53_data_I[i,:])
    GR_r30_g139_interp_I[i,:] = np.interp(all_b, GR_r30_g139_data_b[i,:], GR_r30_g139_data_I[i,:])

Nfiles_SGB = len(SGBfiles_p5_r10_g53)
SGB_p5_r10_g53_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p5_r10_g53_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p5_r10_g139_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p5_r10_g139_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p5_r30_g53_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p5_r30_g53_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p5_r30_g139_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p5_r30_g139_data_I = np.zeros((Nfiles_SGB, 2000))

SGB_p6_r10_g53_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p6_r10_g53_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p6_r10_g139_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p6_r10_g139_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p6_r30_g53_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p6_r30_g53_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p6_r30_g139_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p6_r30_g139_data_I = np.zeros((Nfiles_SGB, 2000))

SGB_p7_r10_g53_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p7_r10_g53_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p7_r10_g139_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p7_r10_g139_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p7_r30_g53_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p7_r30_g53_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p7_r30_g139_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p7_r30_g139_data_I = np.zeros((Nfiles_SGB, 2000))

SGB_p8_r10_g53_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p8_r10_g53_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p8_r10_g139_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p8_r10_g139_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p8_r30_g53_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p8_r30_g53_data_I = np.zeros((Nfiles_SGB, 2000))
SGB_p8_r30_g139_data_b = np.zeros((Nfiles_SGB, 2000))
SGB_p8_r30_g139_data_I = np.zeros((Nfiles_SGB, 2000))

for i in range(0, Nfiles_SGB):
    SGB_p5_r10_g53_data_b[i,:] = np.genfromtxt(SGBfiles_p5_r10_g53[i], delimiter=",", usecols=0)
    SGB_p5_r10_g53_data_I[i,:] = np.genfromtxt(SGBfiles_p5_r10_g53[i], delimiter=",", usecols=1)
    SGB_p5_r10_g139_data_b[i,:] = np.genfromtxt(SGBfiles_p5_r10_g139[i], delimiter=",", usecols=0)
    SGB_p5_r10_g139_data_I[i,:] = np.genfromtxt(SGBfiles_p5_r10_g139[i], delimiter=",", usecols=1)
    SGB_p5_r30_g53_data_b[i,:] = np.genfromtxt(SGBfiles_p5_r30_g53[i], delimiter=",", usecols=0)
    SGB_p5_r30_g53_data_I[i,:] = np.genfromtxt(SGBfiles_p5_r30_g53[i], delimiter=",", usecols=1)
    SGB_p5_r30_g139_data_b[i,:] = np.genfromtxt(SGBfiles_p5_r30_g139[i], delimiter=",", usecols=0)
    SGB_p5_r30_g139_data_I[i,:] = np.genfromtxt(SGBfiles_p5_r30_g139[i], delimiter=",", usecols=1)

    SGB_p6_r10_g53_data_b[i,:] = np.genfromtxt(SGBfiles_p6_r10_g53[i], delimiter=",", usecols=0)
    SGB_p6_r10_g53_data_I[i,:] = np.genfromtxt(SGBfiles_p6_r10_g53[i], delimiter=",", usecols=1)
    SGB_p6_r10_g139_data_b[i,:] = np.genfromtxt(SGBfiles_p6_r10_g139[i], delimiter=",", usecols=0)
    SGB_p6_r10_g139_data_I[i,:] = np.genfromtxt(SGBfiles_p6_r10_g139[i], delimiter=",", usecols=1)
    SGB_p6_r30_g53_data_b[i,:] = np.genfromtxt(SGBfiles_p6_r30_g53[i], delimiter=",", usecols=0)
    SGB_p6_r30_g53_data_I[i,:] = np.genfromtxt(SGBfiles_p6_r30_g53[i], delimiter=",", usecols=1)
    SGB_p6_r30_g139_data_b[i,:] = np.genfromtxt(SGBfiles_p6_r30_g139[i], delimiter=",", usecols=0)
    SGB_p6_r30_g139_data_I[i,:] = np.genfromtxt(SGBfiles_p6_r30_g139[i], delimiter=",", usecols=1)

    SGB_p7_r10_g53_data_b[i,:] = np.genfromtxt(SGBfiles_p7_r10_g53[i], delimiter=",", usecols=0)
    SGB_p7_r10_g53_data_I[i,:] = np.genfromtxt(SGBfiles_p7_r10_g53[i], delimiter=",", usecols=1)
    SGB_p7_r10_g139_data_b[i,:] = np.genfromtxt(SGBfiles_p7_r10_g139[i], delimiter=",", usecols=0)
    SGB_p7_r10_g139_data_I[i,:] = np.genfromtxt(SGBfiles_p7_r10_g139[i], delimiter=",", usecols=1)
    SGB_p7_r30_g53_data_b[i,:] = np.genfromtxt(SGBfiles_p7_r30_g53[i], delimiter=",", usecols=0)
    SGB_p7_r30_g53_data_I[i,:] = np.genfromtxt(SGBfiles_p7_r30_g53[i], delimiter=",", usecols=1)
    SGB_p7_r30_g139_data_b[i,:] = np.genfromtxt(SGBfiles_p7_r30_g139[i], delimiter=",", usecols=0)
    SGB_p7_r30_g139_data_I[i,:] = np.genfromtxt(SGBfiles_p7_r30_g139[i], delimiter=",", usecols=1)

    SGB_p8_r10_g53_data_b[i,:] = np.genfromtxt(SGBfiles_p8_r10_g53[i], delimiter=",", usecols=0)
    SGB_p8_r10_g53_data_I[i,:] = np.genfromtxt(SGBfiles_p8_r10_g53[i], delimiter=",", usecols=1)
    SGB_p8_r10_g139_data_b[i,:] = np.genfromtxt(SGBfiles_p8_r10_g139[i], delimiter=",", usecols=0)
    SGB_p8_r10_g139_data_I[i,:] = np.genfromtxt(SGBfiles_p8_r10_g139[i], delimiter=",", usecols=1)
    SGB_p8_r30_g53_data_b[i,:] = np.genfromtxt(SGBfiles_p8_r30_g53[i], delimiter=",", usecols=0)
    SGB_p8_r30_g53_data_I[i,:] = np.genfromtxt(SGBfiles_p8_r30_g53[i], delimiter=",", usecols=1)
    SGB_p8_r30_g139_data_b[i,:] = np.genfromtxt(SGBfiles_p8_r30_g139[i], delimiter=",", usecols=0)
    SGB_p8_r30_g139_data_I[i,:] = np.genfromtxt(SGBfiles_p8_r30_g139[i], delimiter=",", usecols=1)

SGB_p5_r10_g53_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p5_r10_g139_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p5_r30_g53_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p5_r30_g139_interp_I = np.zeros((Nfiles_SGB, Ninterp))

SGB_p6_r10_g53_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p6_r10_g139_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p6_r30_g53_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p6_r30_g139_interp_I = np.zeros((Nfiles_SGB, Ninterp))

SGB_p7_r10_g53_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p7_r10_g139_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p7_r30_g53_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p7_r30_g139_interp_I = np.zeros((Nfiles_SGB, Ninterp))

SGB_p8_r10_g53_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p8_r10_g139_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p8_r30_g53_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB_p8_r30_g139_interp_I = np.zeros((Nfiles_SGB, Ninterp))

for i in range(0, Nfiles_SGB):
    SGB_p5_r10_g53_interp_I[i,:] = np.interp(all_b, SGB_p5_r10_g53_data_b[i,:], SGB_p5_r10_g53_data_I[i,:])
    SGB_p5_r10_g139_interp_I[i,:] = np.interp(all_b, SGB_p5_r10_g139_data_b[i,:], SGB_p5_r10_g139_data_I[i,:])
    SGB_p5_r30_g53_interp_I[i,:] = np.interp(all_b, SGB_p5_r30_g53_data_b[i,:], SGB_p5_r30_g53_data_I[i,:])
    SGB_p5_r30_g139_interp_I[i,:] = np.interp(all_b, SGB_p5_r30_g139_data_b[i,:], SGB_p5_r30_g139_data_I[i,:])

    SGB_p6_r10_g53_interp_I[i,:] = np.interp(all_b, SGB_p6_r10_g53_data_b[i,:], SGB_p6_r10_g53_data_I[i,:])
    SGB_p6_r10_g139_interp_I[i,:] = np.interp(all_b, SGB_p6_r10_g139_data_b[i,:], SGB_p6_r10_g139_data_I[i,:])
    SGB_p6_r30_g53_interp_I[i,:] = np.interp(all_b, SGB_p6_r30_g53_data_b[i,:], SGB_p6_r30_g53_data_I[i,:])
    SGB_p6_r30_g139_interp_I[i,:] = np.interp(all_b, SGB_p6_r30_g139_data_b[i,:], SGB_p6_r30_g139_data_I[i,:])

    SGB_p7_r10_g53_interp_I[i,:] = np.interp(all_b, SGB_p7_r10_g53_data_b[i,:], SGB_p7_r10_g53_data_I[i,:])
    SGB_p7_r10_g139_interp_I[i,:] = np.interp(all_b, SGB_p7_r10_g139_data_b[i,:], SGB_p7_r10_g139_data_I[i,:])
    SGB_p7_r30_g53_interp_I[i,:] = np.interp(all_b, SGB_p7_r30_g53_data_b[i,:], SGB_p7_r30_g53_data_I[i,:])
    SGB_p7_r30_g139_interp_I[i,:] = np.interp(all_b, SGB_p7_r30_g139_data_b[i,:], SGB_p7_r30_g139_data_I[i,:])

    SGB_p8_r10_g53_interp_I[i,:] = np.interp(all_b, SGB_p8_r10_g53_data_b[i,:], SGB_p8_r10_g53_data_I[i,:])
    SGB_p8_r10_g139_interp_I[i,:]= np.interp(all_b, SGB_p8_r10_g139_data_b[i,:], SGB_p8_r10_g139_data_I[i,:])
    SGB_p8_r30_g53_interp_I[i,:] = np.interp(all_b, SGB_p8_r30_g53_data_b[i,:], SGB_p8_r30_g53_data_I[i,:])
    SGB_p8_r30_g139_interp_I[i,:] = np.interp(all_b, SGB_p8_r30_g139_data_b[i,:], SGB_p8_r30_g139_data_I[i,:])

Nfiles_REZZ = len(pREZZfiles6_r10_g53)
pREZZ_p6_r10_g53_data_b = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p6_r10_g53_data_I = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p6_r10_g139_data_b = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p6_r10_g139_data_I = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p6_r30_g53_data_b = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p6_r30_g53_data_I = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p6_r30_g139_data_b = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p6_r30_g139_data_I = np.zeros((Nfiles_REZZ, 2000))

nREZZ_p6_r10_g53_data_b = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p6_r10_g53_data_I = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p6_r10_g139_data_b = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p6_r10_g139_data_I = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p6_r30_g53_data_b = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p6_r30_g53_data_I = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p6_r30_g139_data_b = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p6_r30_g139_data_I = np.zeros((Nfiles_REZZ, 2000))

pREZZ_p8_r10_g53_data_b = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p8_r10_g53_data_I = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p8_r10_g139_data_b = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p8_r10_g139_data_I = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p8_r30_g53_data_b = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p8_r30_g53_data_I = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p8_r30_g139_data_b = np.zeros((Nfiles_REZZ, 2000))
pREZZ_p8_r30_g139_data_I = np.zeros((Nfiles_REZZ, 2000))

nREZZ_p8_r10_g53_data_b = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p8_r10_g53_data_I = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p8_r10_g139_data_b = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p8_r10_g139_data_I = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p8_r30_g53_data_b = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p8_r30_g53_data_I = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p8_r30_g139_data_b = np.zeros((Nfiles_REZZ, 2000))
nREZZ_p8_r30_g139_data_I = np.zeros((Nfiles_REZZ, 2000))

for i in range(0, Nfiles_REZZ):
    pREZZ_p6_r10_g53_data_b[i,:] = np.genfromtxt(pREZZfiles6_r10_g53[i], delimiter=",", usecols=0)
    pREZZ_p6_r10_g53_data_I[i,:] = np.genfromtxt(pREZZfiles6_r10_g53[i], delimiter=",", usecols=1)
    pREZZ_p6_r10_g139_data_b[i,:] = np.genfromtxt(pREZZfiles6_r10_g139[i], delimiter=",", usecols=0)
    pREZZ_p6_r10_g139_data_I[i,:] = np.genfromtxt(pREZZfiles6_r10_g139[i], delimiter=",", usecols=1)
    pREZZ_p6_r30_g53_data_b[i,:] = np.genfromtxt(pREZZfiles6_r30_g53[i], delimiter=",", usecols=0)
    pREZZ_p6_r30_g53_data_I[i,:] = np.genfromtxt(pREZZfiles6_r30_g53[i], delimiter=",", usecols=1)
    pREZZ_p6_r30_g139_data_b[i,:] = np.genfromtxt(pREZZfiles6_r30_g139[i], delimiter=",", usecols=0)
    pREZZ_p6_r30_g139_data_I[i,:] = np.genfromtxt(pREZZfiles6_r30_g139[i], delimiter=",", usecols=1)

    nREZZ_p6_r10_g53_data_b[i,:] = np.genfromtxt(nREZZfiles6_r10_g53[i], delimiter=",", usecols=0)
    nREZZ_p6_r10_g53_data_I[i,:] = np.genfromtxt(nREZZfiles6_r10_g53[i], delimiter=",", usecols=1)
    nREZZ_p6_r10_g139_data_b[i,:] = np.genfromtxt(nREZZfiles6_r10_g139[i], delimiter=",", usecols=0)
    nREZZ_p6_r10_g139_data_I[i,:] = np.genfromtxt(nREZZfiles6_r10_g139[i], delimiter=",", usecols=1)
    nREZZ_p6_r30_g53_data_b[i,:] = np.genfromtxt(nREZZfiles6_r30_g53[i], delimiter=",", usecols=0)
    nREZZ_p6_r30_g53_data_I[i,:] = np.genfromtxt(nREZZfiles6_r30_g53[i], delimiter=",", usecols=1)
    nREZZ_p6_r30_g139_data_b[i,:] = np.genfromtxt(nREZZfiles6_r30_g139[i], delimiter=",", usecols=0)
    nREZZ_p6_r30_g139_data_I[i,:] = np.genfromtxt(nREZZfiles6_r30_g139[i], delimiter=",", usecols=1)

    pREZZ_p8_r10_g53_data_b[i,:] = np.genfromtxt(pREZZfiles8_r10_g53[i], delimiter=",", usecols=0)
    pREZZ_p8_r10_g53_data_I[i,:] = np.genfromtxt(pREZZfiles8_r10_g53[i], delimiter=",", usecols=1)
    pREZZ_p8_r10_g139_data_b[i,:] = np.genfromtxt(pREZZfiles8_r10_g139[i], delimiter=",", usecols=0)
    pREZZ_p8_r10_g139_data_I[i,:] = np.genfromtxt(pREZZfiles8_r10_g139[i], delimiter=",", usecols=1)
    pREZZ_p8_r30_g53_data_b[i,:] = np.genfromtxt(pREZZfiles8_r30_g53[i], delimiter=",", usecols=0)
    pREZZ_p8_r30_g53_data_I[i,:] = np.genfromtxt(pREZZfiles8_r30_g53[i], delimiter=",", usecols=1)
    pREZZ_p8_r30_g139_data_b[i,:] = np.genfromtxt(pREZZfiles8_r30_g139[i], delimiter=",", usecols=0)
    pREZZ_p8_r30_g139_data_I[i,:] = np.genfromtxt(pREZZfiles8_r30_g139[i], delimiter=",", usecols=1)

    nREZZ_p8_r10_g53_data_b[i,:] = np.genfromtxt(nREZZfiles8_r10_g53[i], delimiter=",", usecols=0)
    nREZZ_p8_r10_g53_data_I[i,:] = np.genfromtxt(nREZZfiles8_r10_g53[i], delimiter=",", usecols=1)
    nREZZ_p8_r10_g139_data_b[i,:] = np.genfromtxt(nREZZfiles8_r10_g139[i], delimiter=",", usecols=0)
    nREZZ_p8_r10_g139_data_I[i,:] = np.genfromtxt(nREZZfiles8_r10_g139[i], delimiter=",", usecols=1)
    nREZZ_p8_r30_g53_data_b[i,:] = np.genfromtxt(nREZZfiles8_r30_g53[i], delimiter=",", usecols=0)
    nREZZ_p8_r30_g53_data_I[i,:] = np.genfromtxt(nREZZfiles8_r30_g53[i], delimiter=",", usecols=1)
    nREZZ_p8_r30_g139_data_b[i,:] = np.genfromtxt(nREZZfiles8_r30_g139[i], delimiter=",", usecols=0)
    nREZZ_p8_r30_g139_data_I[i,:] = np.genfromtxt(nREZZfiles8_r30_g139[i], delimiter=",", usecols=1)

pREZZ_p6_r10_g53_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
pREZZ_p6_r10_g139_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
pREZZ_p6_r30_g53_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
pREZZ_p6_r30_g139_interp_I = np.zeros((Nfiles_REZZ, Ninterp))

nREZZ_p6_r10_g53_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
nREZZ_p6_r10_g139_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
nREZZ_p6_r30_g53_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
nREZZ_p6_r30_g139_interp_I = np.zeros((Nfiles_REZZ, Ninterp))

pREZZ_p8_r10_g53_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
pREZZ_p8_r10_g139_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
pREZZ_p8_r30_g53_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
pREZZ_p8_r30_g139_interp_I = np.zeros((Nfiles_REZZ, Ninterp))

nREZZ_p8_r10_g53_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
nREZZ_p8_r10_g139_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
nREZZ_p8_r30_g53_interp_I = np.zeros((Nfiles_REZZ, Ninterp))
nREZZ_p8_r30_g139_interp_I = np.zeros((Nfiles_REZZ, Ninterp))

for i in range(0, Nfiles_REZZ):
    pREZZ_p6_r10_g53_interp_I[i,:] = np.interp(all_b, pREZZ_p6_r10_g53_data_b[i,:], pREZZ_p6_r10_g53_data_I[i,:])
    pREZZ_p6_r10_g139_interp_I[i,:] = np.interp(all_b, pREZZ_p6_r10_g139_data_b[i,:], pREZZ_p6_r10_g139_data_I[i,:])
    pREZZ_p6_r30_g53_interp_I[i,:] = np.interp(all_b, pREZZ_p6_r30_g53_data_b[i,:], pREZZ_p6_r30_g53_data_I[i,:])
    pREZZ_p6_r30_g139_interp_I[i,:] = np.interp(all_b, pREZZ_p6_r30_g139_data_b[i,:], pREZZ_p6_r30_g139_data_I[i,:])

    nREZZ_p6_r10_g53_interp_I[i,:] = np.interp(all_b, nREZZ_p6_r10_g53_data_b[i,:], nREZZ_p6_r10_g53_data_I[i,:])
    nREZZ_p6_r10_g139_interp_I[i,:] = np.interp(all_b, nREZZ_p6_r10_g139_data_b[i,:], nREZZ_p6_r10_g139_data_I[i,:])
    nREZZ_p6_r30_g53_interp_I[i,:] = np.interp(all_b, nREZZ_p6_r30_g53_data_b[i,:], nREZZ_p6_r30_g53_data_I[i,:])
    nREZZ_p6_r30_g139_interp_I[i,:] = np.interp(all_b, nREZZ_p6_r30_g139_data_b[i,:], nREZZ_p6_r30_g139_data_I[i,:])

    pREZZ_p8_r10_g53_interp_I[i,:] = np.interp(all_b, pREZZ_p8_r10_g53_data_b[i,:], pREZZ_p8_r10_g53_data_I[i,:])
    pREZZ_p8_r10_g139_interp_I[i,:] = np.interp(all_b, pREZZ_p8_r10_g139_data_b[i,:], pREZZ_p8_r10_g139_data_I[i,:])
    pREZZ_p8_r30_g53_interp_I[i,:] = np.interp(all_b, pREZZ_p8_r30_g53_data_b[i,:], pREZZ_p8_r30_g53_data_I[i,:])
    pREZZ_p8_r30_g139_interp_I[i,:] = np.interp(all_b, pREZZ_p8_r30_g139_data_b[i,:], pREZZ_p8_r30_g139_data_I[i,:])

    nREZZ_p8_r10_g53_interp_I[i,:] = np.interp(all_b, nREZZ_p8_r10_g53_data_b[i,:], nREZZ_p8_r10_g53_data_I[i,:])
    nREZZ_p8_r10_g139_interp_I[i,:]= np.interp(all_b, nREZZ_p8_r10_g139_data_b[i,:], nREZZ_p8_r10_g139_data_I[i,:])
    nREZZ_p8_r30_g53_interp_I[i,:] = np.interp(all_b, nREZZ_p8_r30_g53_data_b[i,:], nREZZ_p8_r30_g53_data_I[i,:])
    nREZZ_p8_r30_g139_interp_I[i,:] = np.interp(all_b, nREZZ_p8_r30_g139_data_b[i,:], nREZZ_p8_r30_g139_data_I[i,:])

GRint_r10_g53 = np.zeros((Nfiles_GR, Ninterp))
GRint_r10_g139 = np.zeros((Nfiles_GR, Ninterp))
GRint_r30_g53 = np.zeros((Nfiles_GR, Ninterp))
GRint_r30_g139 = np.zeros((Nfiles_GR, Ninterp))

GRintalt_r10_g53 = np.zeros((Nfiles_GR, Ninterp))
GRintalt_r10_g139 = np.zeros((Nfiles_GR, Ninterp))
GRintalt_r30_g53 = np.zeros((Nfiles_GR, Ninterp))
GRintalt_r30_g139 = np.zeros((Nfiles_GR, Ninterp))

for i in range(0, Nfiles_GR):
    GRint_r10_g53[i,:] = np.pi * 2 * all_b**3 * GR_r10_g53_interp_I[i,:]
    GRint_r10_g139[i,:] = np.pi * 2 * all_b**3 * GR_r10_g139_interp_I[i,:]
    GRint_r30_g53[i,:] = np.pi * 2 * all_b**3 * GR_r30_g53_interp_I[i,:]
    GRint_r30_g139[i,:] = np.pi * 2 * all_b**3 * GR_r30_g139_interp_I[i,:]
    
    GRintalt_r10_g53[i,:] = 2 * np.pi * all_b * GR_r10_g53_interp_I[i,:]
    GRintalt_r10_g139[i,:] = 2 * np.pi * all_b * GR_r10_g139_interp_I[i,:]
    GRintalt_r30_g53[i,:] = 2 * np.pi * all_b * GR_r30_g53_interp_I[i,:]
    GRintalt_r30_g139[i,:] = 2 * np.pi * all_b * GR_r30_g139_interp_I[i,:]

SGBint_p5_r10_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p5_r10_g139 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p5_r30_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p5_r30_g139 = np.zeros((Nfiles_SGB, Ninterp))

SGBint_p6_r10_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p6_r10_g139 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p6_r30_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p6_r30_g139 = np.zeros((Nfiles_SGB, Ninterp))

SGBint_p7_r10_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p7_r10_g139 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p7_r30_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p7_r30_g139 = np.zeros((Nfiles_SGB, Ninterp))

SGBint_p8_r10_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p8_r10_g139 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p8_r30_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBint_p8_r30_g139 = np.zeros((Nfiles_SGB, Ninterp))

SGBintalt_p5_r10_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p5_r10_g139 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p5_r30_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p5_r30_g139 = np.zeros((Nfiles_SGB, Ninterp))

SGBintalt_p6_r10_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p6_r10_g139 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p6_r30_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p6_r30_g139 = np.zeros((Nfiles_SGB, Ninterp))

SGBintalt_p7_r10_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p7_r10_g139 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p7_r30_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p7_r30_g139 = np.zeros((Nfiles_SGB, Ninterp))

SGBintalt_p8_r10_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p8_r10_g139 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p8_r30_g53 = np.zeros((Nfiles_SGB, Ninterp))
SGBintalt_p8_r30_g139 = np.zeros((Nfiles_SGB, Ninterp))

for i in range(0, Nfiles_SGB):
    SGBint_p5_r10_g53[i,:] = np.pi * 2 * all_b**3 * SGB_p5_r10_g53_interp_I[i,:]
    SGBint_p5_r10_g139[i,:] = np.pi * 2 * all_b**3 * SGB_p5_r10_g139_interp_I[i,:]
    SGBint_p5_r30_g53[i,:] = np.pi * 2 * all_b**3 * SGB_p5_r30_g53_interp_I[i,:]
    SGBint_p5_r30_g139[i,:] = np.pi * 2 * all_b**3 * SGB_p5_r30_g139_interp_I[i,:]

    SGBint_p6_r10_g53[i,:] = np.pi * 2 * all_b**3 * SGB_p6_r10_g53_interp_I[i,:]
    SGBint_p6_r10_g139[i,:] = np.pi * 2 * all_b**3 * SGB_p6_r10_g139_interp_I[i,:]
    SGBint_p6_r30_g53[i,:] = np.pi * 2 * all_b**3 * SGB_p6_r30_g53_interp_I[i,:]
    SGBint_p6_r30_g139[i,:] = np.pi * 2 * all_b**3 * SGB_p6_r30_g139_interp_I[i,:]

    SGBint_p7_r10_g53[i,:] = np.pi * 2 * all_b**3 * SGB_p7_r10_g53_interp_I[i,:]
    SGBint_p7_r10_g139[i,:] = np.pi * 2 * all_b**3 * SGB_p7_r10_g139_interp_I[i,:]
    SGBint_p7_r30_g53[i,:] = np.pi * 2 * all_b**3 * SGB_p7_r30_g53_interp_I[i,:]
    SGBint_p7_r30_g139[i,:] = np.pi * 2 * all_b**3 * SGB_p7_r30_g139_interp_I[i,:]

    SGBint_p8_r10_g53[i,:] = np.pi * 2 * all_b**3 * SGB_p8_r10_g53_interp_I[i,:]
    SGBint_p8_r10_g139[i,:] = np.pi * 2 * all_b**3 * SGB_p8_r10_g139_interp_I[i,:]
    SGBint_p8_r30_g53[i,:] = np.pi * 2 * all_b**3 * SGB_p8_r30_g53_interp_I[i,:]
    SGBint_p8_r30_g139[i,:] = np.pi * 2 * all_b**3 * SGB_p8_r30_g139_interp_I[i,:]
    
    SGBintalt_p5_r10_g53[i,:] = np.pi * 2 * all_b * SGB_p5_r10_g53_interp_I[i,:]
    SGBintalt_p5_r10_g139[i,:] = np.pi * 2 * all_b * SGB_p5_r10_g139_interp_I[i,:]
    SGBintalt_p5_r30_g53[i,:] = np.pi * 2 * all_b * SGB_p5_r30_g53_interp_I[i,:]
    SGBintalt_p5_r30_g139[i,:] = np.pi * 2 * all_b * SGB_p5_r30_g139_interp_I[i,:]

    SGBintalt_p6_r10_g53[i,:] = np.pi * 2 * all_b * SGB_p6_r10_g53_interp_I[i,:]
    SGBintalt_p6_r10_g139[i,:] = np.pi * 2 * all_b * SGB_p6_r10_g139_interp_I[i,:]
    SGBintalt_p6_r30_g53[i,:] = np.pi * 2 * all_b * SGB_p6_r30_g53_interp_I[i,:]
    SGBintalt_p6_r30_g139[i,:] = np.pi * 2 * all_b * SGB_p6_r30_g139_interp_I[i,:]

    SGBintalt_p7_r10_g53[i,:] = np.pi * 2 * all_b * SGB_p7_r10_g53_interp_I[i,:]
    SGBintalt_p7_r10_g139[i,:] = np.pi * 2 * all_b * SGB_p7_r10_g139_interp_I[i,:]
    SGBintalt_p7_r30_g53[i,:] = np.pi * 2 * all_b * SGB_p7_r30_g53_interp_I[i,:]
    SGBintalt_p7_r30_g139[i,:] = np.pi * 2 * all_b * SGB_p7_r30_g139_interp_I[i,:]

    SGBintalt_p8_r10_g53[i,:] = np.pi * 2 * all_b * SGB_p8_r10_g53_interp_I[i,:]
    SGBintalt_p8_r10_g139[i,:] = np.pi * 2 * all_b * SGB_p8_r10_g139_interp_I[i,:]
    SGBintalt_p8_r30_g53[i,:] = np.pi * 2 * all_b * SGB_p8_r30_g53_interp_I[i,:]
    SGBintalt_p8_r30_g139[i,:] = np.pi * 2 * all_b * SGB_p8_r30_g139_interp_I[i,:]

pREZZint_p6_r10_g53 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZint_p6_r10_g139 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZint_p6_r30_g53 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZint_p6_r30_g139 = np.zeros((Nfiles_REZZ, Ninterp))

nREZZint_p6_r10_g53 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZint_p6_r10_g139 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZint_p6_r30_g53 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZint_p6_r30_g139 = np.zeros((Nfiles_REZZ, Ninterp))

pREZZint_p8_r10_g53 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZint_p8_r10_g139 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZint_p8_r30_g53 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZint_p8_r30_g139 = np.zeros((Nfiles_REZZ, Ninterp))

nREZZint_p8_r10_g53 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZint_p8_r10_g139 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZint_p8_r30_g53 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZint_p8_r30_g139 = np.zeros((Nfiles_REZZ, Ninterp))

pREZZintalt_p6_r10_g53 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZintalt_p6_r10_g139 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZintalt_p6_r30_g53 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZintalt_p6_r30_g139 = np.zeros((Nfiles_REZZ, Ninterp))

nREZZintalt_p6_r10_g53 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZintalt_p6_r10_g139 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZintalt_p6_r30_g53 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZintalt_p6_r30_g139 = np.zeros((Nfiles_REZZ, Ninterp))

pREZZintalt_p8_r10_g53 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZintalt_p8_r10_g139 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZintalt_p8_r30_g53 = np.zeros((Nfiles_REZZ, Ninterp))
pREZZintalt_p8_r30_g139 = np.zeros((Nfiles_REZZ, Ninterp))

nREZZintalt_p8_r10_g53 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZintalt_p8_r10_g139 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZintalt_p8_r30_g53 = np.zeros((Nfiles_REZZ, Ninterp))
nREZZintalt_p8_r30_g139 = np.zeros((Nfiles_REZZ, Ninterp))

for i in range(0, Nfiles_REZZ):
    pREZZint_p6_r10_g53[i,:] = np.pi * 2 * all_b**3 * pREZZ_p6_r10_g53_interp_I[i,:]
    pREZZint_p6_r10_g139[i,:] = np.pi * 2 * all_b**3 * pREZZ_p6_r10_g139_interp_I[i,:]
    pREZZint_p6_r30_g53[i,:] = np.pi * 2 * all_b**3 * pREZZ_p6_r30_g53_interp_I[i,:]
    pREZZint_p6_r30_g139[i,:] = np.pi * 2 * all_b**3 * pREZZ_p6_r30_g139_interp_I[i,:]

    nREZZint_p6_r10_g53[i,:] = np.pi * 2 * all_b**3 * nREZZ_p6_r10_g53_interp_I[i,:]
    nREZZint_p6_r10_g139[i,:] = np.pi * 2 * all_b**3 * nREZZ_p6_r10_g139_interp_I[i,:]
    nREZZint_p6_r30_g53[i,:] = np.pi * 2 * all_b**3 * nREZZ_p6_r30_g53_interp_I[i,:]
    nREZZint_p6_r30_g139[i,:] = np.pi * 2 * all_b**3 * nREZZ_p6_r30_g139_interp_I[i,:]

    pREZZint_p8_r10_g53[i,:] = np.pi * 2 * all_b**3 * pREZZ_p8_r10_g53_interp_I[i,:]
    pREZZint_p8_r10_g139[i,:] = np.pi * 2 * all_b**3 * pREZZ_p8_r10_g139_interp_I[i,:]
    pREZZint_p8_r30_g53[i,:] = np.pi * 2 * all_b**3 * pREZZ_p8_r30_g53_interp_I[i,:]
    pREZZint_p8_r30_g139[i,:] = np.pi * 2 * all_b**3 * pREZZ_p8_r30_g139_interp_I[i,:]

    nREZZint_p8_r10_g53[i,:] = np.pi * 2 * all_b**3 * nREZZ_p8_r10_g53_interp_I[i,:]
    nREZZint_p8_r10_g139[i,:] = np.pi * 2 * all_b**3 * nREZZ_p8_r10_g139_interp_I[i,:]
    nREZZint_p8_r30_g53[i,:] = np.pi * 2 * all_b**3 * nREZZ_p8_r30_g53_interp_I[i,:]
    nREZZint_p8_r30_g139[i,:] = np.pi * 2 * all_b**3 * nREZZ_p8_r30_g139_interp_I[i,:]

    pREZZintalt_p6_r10_g53[i,:] = np.pi * 2 * all_b * pREZZ_p6_r10_g53_interp_I[i,:]
    pREZZintalt_p6_r10_g139[i,:] = np.pi * 2 * all_b * pREZZ_p6_r10_g139_interp_I[i,:]
    pREZZintalt_p6_r30_g53[i,:] = np.pi * 2 * all_b * pREZZ_p6_r30_g53_interp_I[i,:]
    pREZZintalt_p6_r30_g139[i,:] = np.pi * 2 * all_b * pREZZ_p6_r30_g139_interp_I[i,:]

    nREZZintalt_p6_r10_g53[i,:] = np.pi * 2 * all_b * nREZZ_p6_r10_g53_interp_I[i,:]
    nREZZintalt_p6_r10_g139[i,:] = np.pi * 2 * all_b * nREZZ_p6_r10_g139_interp_I[i,:]
    nREZZintalt_p6_r30_g53[i,:] = np.pi * 2 * all_b * nREZZ_p6_r30_g53_interp_I[i,:]
    nREZZintalt_p6_r30_g139[i,:] = np.pi * 2 * all_b * nREZZ_p6_r30_g139_interp_I[i,:]

    pREZZintalt_p8_r10_g53[i,:] = np.pi * 2 * all_b * pREZZ_p8_r10_g53_interp_I[i,:]
    pREZZintalt_p8_r10_g139[i,:] = np.pi * 2 * all_b * pREZZ_p8_r10_g139_interp_I[i,:]
    pREZZintalt_p8_r30_g53[i,:] = np.pi * 2 * all_b * pREZZ_p8_r30_g53_interp_I[i,:]
    pREZZintalt_p8_r30_g139[i,:] = np.pi * 2 * all_b * pREZZ_p8_r30_g139_interp_I[i,:]

    nREZZintalt_p8_r10_g53[i,:] = np.pi * 2 * all_b * nREZZ_p8_r10_g53_interp_I[i,:]
    nREZZintalt_p8_r10_g139[i,:] = np.pi * 2 * all_b * nREZZ_p8_r10_g139_interp_I[i,:]
    nREZZintalt_p8_r30_g53[i,:] = np.pi * 2 * all_b * nREZZ_p8_r30_g53_interp_I[i,:]
    nREZZintalt_p8_r30_g139[i,:] = np.pi * 2 * all_b * nREZZ_p8_r30_g139_interp_I[i,:]

GRrchar_r10_g53 = np.zeros(Nfiles_GR)
GRrchar_r10_g139 = np.zeros(Nfiles_GR)
GRrchar_r30_g53 = np.zeros(Nfiles_GR)
GRrchar_r30_g139 = np.zeros(Nfiles_GR)

GRrcharalt_r10_g53 = np.zeros(Nfiles_GR)
GRrcharalt_r10_g139 = np.zeros(Nfiles_GR)
GRrcharalt_r30_g53 = np.zeros(Nfiles_GR)
GRrcharalt_r30_g139 = np.zeros(Nfiles_GR)

tol = 5 * 10**(-4)

for i in range(0, Nfiles_GR):
    GRrchar_r10_g53[i] = (getIntegration(GRint_r10_g53[i,:], all_b) * (getIntegration(GRintalt_r10_g53[i,:], all_b))**(-1))**(0.5)
    GRrchar_r10_g139[i] = (getIntegration(GRint_r10_g139[i,:], all_b) * (getIntegration(GRintalt_r10_g139[i,:], all_b))**(-1))**(0.5)
    GRrchar_r30_g53[i] = (getIntegration(GRint_r30_g53[i,:], all_b) * (getIntegration(GRintalt_r30_g53[i,:], all_b))**(-1))**(0.5)
    GRrchar_r30_g139[i] = (getIntegration(GRint_r30_g139[i,:], all_b) * (getIntegration(GRintalt_r30_g139[i,:], all_b))**(-1))**(0.5)

    GRrcharalt_r10_g53[i] = getValue(getIntegration, GRintalt_r10_g53[i,:], all_b, tol, [getIntegration(GRintalt_r10_g53[i,:], all_b)] )
    GRrcharalt_r10_g139[i] = getValue(getIntegration, GRintalt_r10_g139[i,:], all_b, tol, [getIntegration(GRintalt_r10_g139[i,:], all_b)] )
    GRrcharalt_r30_g53[i] = getValue(getIntegration, GRintalt_r30_g53[i,:], all_b, tol, [getIntegration(GRintalt_r30_g53[i,:], all_b)] )
    GRrcharalt_r30_g139[i] = getValue(getIntegration, GRintalt_r30_g139[i,:], all_b, tol, [getIntegration(GRintalt_r30_g139[i,:], all_b)] )

SGBrchar_p5_r10_g53 = np.zeros(Nfiles_SGB)
SGBrchar_p5_r10_g139 = np.zeros(Nfiles_SGB)
SGBrchar_p5_r30_g53 = np.zeros(Nfiles_SGB)
SGBrchar_p5_r30_g139 = np.zeros(Nfiles_SGB)

SGBrchar_p6_r10_g53 = np.zeros(Nfiles_SGB)
SGBrchar_p6_r10_g139 = np.zeros(Nfiles_SGB)
SGBrchar_p6_r30_g53 = np.zeros(Nfiles_SGB)
SGBrchar_p6_r30_g139 = np.zeros(Nfiles_SGB)

SGBrchar_p7_r10_g53 = np.zeros(Nfiles_SGB)
SGBrchar_p7_r10_g139 = np.zeros(Nfiles_SGB)
SGBrchar_p7_r30_g53 = np.zeros(Nfiles_SGB)
SGBrchar_p7_r30_g139 = np.zeros(Nfiles_SGB)

SGBrchar_p8_r10_g53 = np.zeros(Nfiles_SGB)
SGBrchar_p8_r10_g139 = np.zeros(Nfiles_SGB)
SGBrchar_p8_r30_g53 = np.zeros(Nfiles_SGB)
SGBrchar_p8_r30_g139 = np.zeros(Nfiles_SGB)

SGBrcharalt_p5_r10_g53 = np.zeros(Nfiles_SGB)
SGBrcharalt_p5_r10_g139 = np.zeros(Nfiles_SGB)
SGBrcharalt_p5_r30_g53 = np.zeros(Nfiles_SGB)
SGBrcharalt_p5_r30_g139 = np.zeros(Nfiles_SGB)

SGBrcharalt_p6_r10_g53 = np.zeros(Nfiles_SGB)
SGBrcharalt_p6_r10_g139 = np.zeros(Nfiles_SGB)
SGBrcharalt_p6_r30_g53 = np.zeros(Nfiles_SGB)
SGBrcharalt_p6_r30_g139 = np.zeros(Nfiles_SGB)

SGBrcharalt_p7_r10_g53 = np.zeros(Nfiles_SGB)
SGBrcharalt_p7_r10_g139 = np.zeros(Nfiles_SGB)
SGBrcharalt_p7_r30_g53 = np.zeros(Nfiles_SGB)
SGBrcharalt_p7_r30_g139 = np.zeros(Nfiles_SGB)

SGBrcharalt_p8_r10_g53 = np.zeros(Nfiles_SGB)
SGBrcharalt_p8_r10_g139 = np.zeros(Nfiles_SGB)
SGBrcharalt_p8_r30_g53 = np.zeros(Nfiles_SGB)
SGBrcharalt_p8_r30_g139 = np.zeros(Nfiles_SGB)

for i in range(0, Nfiles_SGB):
    SGBrchar_p5_r10_g53[i] = (getIntegration(SGBint_p5_r10_g53[i,:], all_b) * (getIntegration(SGBintalt_p5_r10_g53[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p5_r10_g139[i] = (getIntegration(SGBint_p5_r10_g139[i,:], all_b) * (getIntegration(SGBintalt_p5_r10_g139[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p5_r30_g53[i] = (getIntegration(SGBint_p5_r30_g53[i,:], all_b) * (getIntegration(SGBintalt_p5_r30_g53[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p5_r30_g139[i] = (getIntegration(SGBint_p5_r30_g139[i,:], all_b) * (getIntegration(SGBintalt_p5_r30_g139[i,:], all_b))**(-1))**(0.5)

    SGBrchar_p6_r10_g53[i] = (getIntegration(SGBint_p6_r10_g53[i,:], all_b) * (getIntegration(SGBintalt_p6_r10_g53[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p6_r10_g139[i] = (getIntegration(SGBint_p6_r10_g139[i,:], all_b) * (getIntegration(SGBintalt_p6_r10_g139[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p6_r30_g53[i] = (getIntegration(SGBint_p6_r30_g53[i,:], all_b) * (getIntegration(SGBintalt_p6_r30_g53[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p6_r30_g139[i] = (getIntegration(SGBint_p6_r30_g139[i,:], all_b) * (getIntegration(SGBintalt_p6_r30_g139[i,:], all_b))**(-1))**(0.5)

    SGBrchar_p7_r10_g53[i] = (getIntegration(SGBint_p7_r10_g53[i,:], all_b) * (getIntegration(SGBintalt_p7_r10_g53[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p7_r10_g139[i] = (getIntegration(SGBint_p7_r10_g139[i,:], all_b) * (getIntegration(SGBintalt_p7_r10_g139[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p7_r30_g53[i] = (getIntegration(SGBint_p7_r30_g53[i,:], all_b) * (getIntegration(SGBintalt_p7_r30_g53[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p7_r30_g139[i] = (getIntegration(SGBint_p7_r30_g139[i,:], all_b) * (getIntegration(SGBintalt_p7_r30_g139[i,:], all_b))**(-1))**(0.5)

    SGBrchar_p8_r10_g53[i] = (getIntegration(SGBint_p8_r10_g53[i,:], all_b) * (getIntegration(SGBintalt_p8_r10_g53[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p8_r10_g139[i] = (getIntegration(SGBint_p8_r10_g139[i,:], all_b) * (getIntegration(SGBintalt_p8_r10_g139[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p8_r30_g53[i] = (getIntegration(SGBint_p8_r30_g53[i,:], all_b) * (getIntegration(SGBintalt_p8_r30_g53[i,:], all_b))**(-1))**(0.5)
    SGBrchar_p8_r30_g139[i] = (getIntegration(SGBint_p8_r30_g139[i,:], all_b) * (getIntegration(SGBintalt_p8_r30_g139[i,:], all_b))**(-1))**(0.5)
    
    SGBrcharalt_p5_r10_g53[i] = getValue(getIntegration, SGBintalt_p5_r10_g53[i,:], all_b, tol, [getIntegration(SGBintalt_p5_r10_g53[i,:], all_b)])
    SGBrcharalt_p5_r10_g139[i] = getValue(getIntegration, SGBintalt_p5_r10_g139[i,:], all_b, tol, [getIntegration(SGBintalt_p5_r10_g139[i,:], all_b)])
    SGBrcharalt_p5_r30_g53[i] = getValue(getIntegration, SGBintalt_p5_r30_g53[i,:], all_b, tol, [getIntegration(SGBintalt_p5_r30_g53[i,:], all_b)])
    SGBrcharalt_p5_r30_g139[i] = getValue(getIntegration, SGBintalt_p5_r30_g139[i,:], all_b, tol, [getIntegration(SGBintalt_p5_r30_g139[i,:], all_b)])

    SGBrcharalt_p6_r10_g53[i] = getValue(getIntegration, SGBintalt_p6_r10_g53[i,:], all_b, tol, [getIntegration(SGBintalt_p6_r10_g53[i,:], all_b)])
    SGBrcharalt_p6_r10_g139[i] = getValue(getIntegration, SGBintalt_p6_r10_g139[i,:], all_b, tol, [getIntegration(SGBintalt_p6_r10_g139[i,:], all_b)])
    SGBrcharalt_p6_r30_g53[i] = getValue(getIntegration, SGBintalt_p6_r30_g53[i,:], all_b, tol, [getIntegration(SGBintalt_p6_r30_g53[i,:], all_b)])
    SGBrcharalt_p6_r30_g139[i] = getValue(getIntegration, SGBintalt_p6_r30_g139[i,:], all_b, tol, [getIntegration(SGBintalt_p6_r30_g139[i,:], all_b)])

    SGBrcharalt_p7_r10_g53[i] = getValue(getIntegration, SGBintalt_p7_r10_g53[i,:], all_b, tol, [getIntegration(SGBintalt_p7_r10_g53[i,:], all_b)])
    SGBrcharalt_p7_r10_g139[i] = getValue(getIntegration, SGBintalt_p7_r10_g139[i,:], all_b, tol, [getIntegration(SGBintalt_p7_r10_g139[i,:], all_b)])
    SGBrcharalt_p7_r30_g53[i] = getValue(getIntegration, SGBintalt_p7_r30_g53[i,:], all_b, tol, [getIntegration(SGBintalt_p7_r30_g53[i,:], all_b)])
    SGBrcharalt_p7_r30_g139[i] = getValue(getIntegration, SGBintalt_p7_r30_g139[i,:], all_b, tol, [getIntegration(SGBintalt_p7_r30_g139[i,:], all_b)])

    SGBrcharalt_p8_r10_g53[i] = getValue(getIntegration, SGBintalt_p8_r10_g53[i,:], all_b, tol, [getIntegration(SGBintalt_p8_r10_g53[i,:], all_b)])
    SGBrcharalt_p8_r10_g139[i] = getValue(getIntegration, SGBintalt_p8_r10_g139[i,:], all_b, tol, [getIntegration(SGBintalt_p8_r10_g139[i,:], all_b)])
    SGBrcharalt_p8_r30_g53[i] = getValue(getIntegration, SGBintalt_p8_r30_g53[i,:], all_b, tol, [getIntegration(SGBintalt_p8_r30_g53[i,:], all_b)])
    SGBrcharalt_p8_r30_g139[i] = getValue(getIntegration, SGBintalt_p8_r30_g139[i,:], all_b, tol, [getIntegration(SGBintalt_p8_r30_g139[i,:], all_b)])

pREZZrchar_p6_r10_g53 = np.zeros(Nfiles_REZZ)
pREZZrchar_p6_r10_g139 = np.zeros(Nfiles_REZZ)
pREZZrchar_p6_r30_g53 = np.zeros(Nfiles_REZZ)
pREZZrchar_p6_r30_g139 = np.zeros(Nfiles_REZZ)

nREZZrchar_p6_r10_g53 = np.zeros(Nfiles_REZZ)
nREZZrchar_p6_r10_g139 = np.zeros(Nfiles_REZZ)
nREZZrchar_p6_r30_g53 = np.zeros(Nfiles_REZZ)
nREZZrchar_p6_r30_g139 = np.zeros(Nfiles_REZZ)

pREZZrchar_p8_r10_g53 = np.zeros(Nfiles_REZZ)
pREZZrchar_p8_r10_g139 = np.zeros(Nfiles_REZZ)
pREZZrchar_p8_r30_g53 = np.zeros(Nfiles_REZZ)
pREZZrchar_p8_r30_g139 = np.zeros(Nfiles_REZZ)

nREZZrchar_p8_r10_g53 = np.zeros(Nfiles_REZZ)
nREZZrchar_p8_r10_g139 = np.zeros(Nfiles_REZZ)
nREZZrchar_p8_r30_g53 = np.zeros(Nfiles_REZZ)
nREZZrchar_p8_r30_g139 = np.zeros(Nfiles_REZZ)

for i in range(0, Nfiles_REZZ):
    pREZZrchar_p6_r10_g53[i] = (getIntegration(pREZZint_p6_r10_g53[i,:], all_b) * (getIntegration(pREZZintalt_p6_r10_g53[i,:], all_b))**(-1))**(0.5)
    pREZZrchar_p6_r10_g139[i] = (getIntegration(pREZZint_p6_r10_g139[i,:], all_b) * (getIntegration(pREZZintalt_p6_r10_g139[i,:], all_b))**(-1))**(0.5)
    pREZZrchar_p6_r30_g53[i] = (getIntegration(pREZZint_p6_r30_g53[i,:], all_b) * (getIntegration(pREZZintalt_p6_r30_g53[i,:], all_b))**(-1))**(0.5)
    pREZZrchar_p6_r30_g139[i] = (getIntegration(pREZZint_p6_r30_g139[i,:], all_b) * (getIntegration(pREZZintalt_p6_r30_g139[i,:], all_b))**(-1))**(0.5)

    nREZZrchar_p6_r10_g53[i] = (getIntegration(nREZZint_p6_r10_g53[i,:], all_b) * (getIntegration(nREZZintalt_p6_r10_g53[i,:], all_b))**(-1))**(0.5)
    nREZZrchar_p6_r10_g139[i] = (getIntegration(nREZZint_p6_r10_g139[i,:], all_b) * (getIntegration(nREZZintalt_p6_r10_g139[i,:], all_b))**(-1))**(0.5)
    nREZZrchar_p6_r30_g53[i] = (getIntegration(nREZZint_p6_r30_g53[i,:], all_b) * (getIntegration(nREZZintalt_p6_r30_g53[i,:], all_b))**(-1))**(0.5)
    nREZZrchar_p6_r30_g139[i] = (getIntegration(nREZZint_p6_r30_g139[i,:], all_b) * (getIntegration(nREZZintalt_p6_r30_g139[i,:], all_b))**(-1))**(0.5)

    pREZZrchar_p8_r10_g53[i] = (getIntegration(pREZZint_p8_r10_g53[i,:], all_b) * (getIntegration(pREZZintalt_p8_r10_g53[i,:], all_b))**(-1))**(0.5)
    pREZZrchar_p8_r10_g139[i] = (getIntegration(pREZZint_p8_r10_g139[i,:], all_b) * (getIntegration(pREZZintalt_p8_r10_g139[i,:], all_b))**(-1))**(0.5)
    pREZZrchar_p8_r30_g53[i] = (getIntegration(pREZZint_p8_r30_g53[i,:], all_b) * (getIntegration(pREZZintalt_p8_r30_g53[i,:], all_b))**(-1))**(0.5)
    pREZZrchar_p8_r30_g139[i] = (getIntegration(pREZZint_p8_r30_g139[i,:], all_b) * (getIntegration(pREZZintalt_p8_r30_g139[i,:], all_b))**(-1))**(0.5)

    nREZZrchar_p8_r10_g53[i] = (getIntegration(nREZZint_p8_r10_g53[i,:], all_b) * (getIntegration(nREZZintalt_p8_r10_g53[i,:], all_b))**(-1))**(0.5)
    nREZZrchar_p8_r10_g139[i] = (getIntegration(nREZZint_p8_r10_g139[i,:], all_b) * (getIntegration(nREZZintalt_p8_r10_g139[i,:], all_b))**(-1))**(0.5)
    nREZZrchar_p8_r30_g53[i] = (getIntegration(nREZZint_p8_r30_g53[i,:], all_b) * (getIntegration(nREZZintalt_p8_r30_g53[i,:], all_b))**(-1))**(0.5)
    nREZZrchar_p8_r30_g139[i] = (getIntegration(nREZZint_p8_r30_g139[i,:], all_b) * (getIntegration(nREZZintalt_p8_r30_g139[i,:], all_b))**(-1))**(0.5)

REZZrchar_p6_r10_g53 = np.hstack((nREZZrchar_p6_r10_g53, pREZZrchar_p6_r10_g53))
REZZrchar_p6_r10_g139 = np.hstack((nREZZrchar_p6_r10_g139, pREZZrchar_p6_r10_g139))
REZZrchar_p6_r30_g53 = np.hstack((nREZZrchar_p6_r30_g53, pREZZrchar_p6_r30_g53))
REZZrchar_p6_r30_g139 = np.hstack((nREZZrchar_p6_r30_g139, pREZZrchar_p6_r30_g139))

REZZrchar_p8_r10_g53 = np.hstack((nREZZrchar_p8_r10_g53, pREZZrchar_p8_r10_g53))
REZZrchar_p8_r10_g139 = np.hstack((nREZZrchar_p8_r10_g139, pREZZrchar_p8_r10_g139))
REZZrchar_p8_r30_g53 = np.hstack((nREZZrchar_p8_r30_g53, pREZZrchar_p8_r30_g53))
REZZrchar_p8_r30_g139 = np.hstack((nREZZrchar_p8_r30_g139, pREZZrchar_p8_r30_g139))

GRtheta_r10_g53 = np.zeros(Nfiles_GR)
GRtheta_r10_g139 = np.zeros(Nfiles_GR)
GRtheta_r30_g53 = np.zeros(Nfiles_GR)
GRtheta_r30_g139 = np.zeros(Nfiles_GR)

GRthetaalt_r10_g53 = np.zeros(Nfiles_GR)
GRthetaalt_r10_g139 = np.zeros(Nfiles_GR)
GRthetaalt_r30_g53 = np.zeros(Nfiles_GR)
GRthetaalt_r30_g139 = np.zeros(Nfiles_GR)

for i in range(0, Nfiles_GR):
    GRtheta_r10_g53[i] = GRrchar_r10_g53[i] * getCriticalB(0)**(-1) # GR
    GRtheta_r10_g139[i] = GRrchar_r10_g139[i] * getCriticalB(0)**(-1) # GR
    GRtheta_r30_g53[i] = GRrchar_r30_g53[i] * getCriticalB(0)**(-1) # GR
    GRtheta_r30_g139[i] = GRrchar_r30_g139[i] * getCriticalB(0)**(-1) # GR

    GRthetaalt_r10_g53[i] = GRrcharalt_r10_g53[i] * getCriticalB(0)**(-1) 
    GRthetaalt_r10_g139[i] = GRrcharalt_r10_g139[i] * getCriticalB(0)**(-1) 
    GRthetaalt_r30_g53[i] = GRrcharalt_r30_g53[i] * getCriticalB(0)**(-1) 
    GRthetaalt_r30_g139[i] = GRrcharalt_r30_g139[i] * getCriticalB(0)**(-1) 

SGBtheta_p5_r10_g53 = np.zeros(Nfiles_SGB)
SGBtheta_p5_r10_g139 = np.zeros(Nfiles_SGB)
SGBtheta_p5_r30_g53 = np.zeros(Nfiles_SGB)
SGBtheta_p5_r30_g139 = np.zeros(Nfiles_SGB)

SGBtheta_p6_r10_g53 = np.zeros(Nfiles_SGB)
SGBtheta_p6_r10_g139 = np.zeros(Nfiles_SGB)
SGBtheta_p6_r30_g53 = np.zeros(Nfiles_SGB)
SGBtheta_p6_r30_g139 = np.zeros(Nfiles_SGB)

SGBtheta_p7_r10_g53 = np.zeros(Nfiles_SGB)
SGBtheta_p7_r10_g139 = np.zeros(Nfiles_SGB)
SGBtheta_p7_r30_g53 = np.zeros(Nfiles_SGB)
SGBtheta_p7_r30_g139 = np.zeros(Nfiles_SGB)

SGBtheta_p8_r10_g53 = np.zeros(Nfiles_SGB)
SGBtheta_p8_r10_g139 = np.zeros(Nfiles_SGB)
SGBtheta_p8_r30_g53 = np.zeros(Nfiles_SGB)
SGBtheta_p8_r30_g139 = np.zeros(Nfiles_SGB)

SGBthetaalt_p5_r10_g53 = np.zeros(Nfiles_SGB)
SGBthetaalt_p5_r10_g139 = np.zeros(Nfiles_SGB)
SGBthetaalt_p5_r30_g53 = np.zeros(Nfiles_SGB)
SGBthetaalt_p5_r30_g139 = np.zeros(Nfiles_SGB)

SGBthetaalt_p6_r10_g53 = np.zeros(Nfiles_SGB)
SGBthetaalt_p6_r10_g139 = np.zeros(Nfiles_SGB)
SGBthetaalt_p6_r30_g53 = np.zeros(Nfiles_SGB)
SGBthetaalt_p6_r30_g139 = np.zeros(Nfiles_SGB)

SGBthetaalt_p7_r10_g53 = np.zeros(Nfiles_SGB)
SGBthetaalt_p7_r10_g139 = np.zeros(Nfiles_SGB)
SGBthetaalt_p7_r30_g53 = np.zeros(Nfiles_SGB)
SGBthetaalt_p7_r30_g139 = np.zeros(Nfiles_SGB)

SGBthetaalt_p8_r10_g53 = np.zeros(Nfiles_SGB)
SGBthetaalt_p8_r10_g139 = np.zeros(Nfiles_SGB)
SGBthetaalt_p8_r30_g53 = np.zeros(Nfiles_SGB)
SGBthetaalt_p8_r30_g139 = np.zeros(Nfiles_SGB)

for i in range(0, Nfiles_SGB):
    SGBtheta_p5_r10_g53[i] = SGBrchar_p5_r10_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p5_r10_g139[i] = SGBrchar_p5_r10_g139[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p5_r30_g53[i] = SGBrchar_p5_r30_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p5_r30_g139[i] = SGBrchar_p5_r30_g139[i] * getCriticalB(zetas[i])**(-1)

    SGBtheta_p6_r10_g53[i] = SGBrchar_p6_r10_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p6_r10_g139[i] = SGBrchar_p6_r10_g139[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p6_r30_g53[i] = SGBrchar_p6_r30_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p6_r30_g139[i] = SGBrchar_p6_r30_g139[i] * getCriticalB(zetas[i])**(-1)

    SGBtheta_p7_r10_g53[i] = SGBrchar_p7_r10_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p7_r10_g139[i] = SGBrchar_p7_r10_g139[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p7_r30_g53[i] = SGBrchar_p7_r30_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p7_r30_g139[i] = SGBrchar_p7_r30_g139[i] * getCriticalB(zetas[i])**(-1)

    SGBtheta_p8_r10_g53[i] = SGBrchar_p8_r10_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p8_r10_g139[i] = SGBrchar_p8_r10_g139[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p8_r30_g53[i] = SGBrchar_p8_r30_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBtheta_p8_r30_g139[i] = SGBrchar_p8_r30_g139[i] * getCriticalB(zetas[i])**(-1)

    SGBthetaalt_p5_r10_g53[i] = SGBrcharalt_p5_r10_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p5_r10_g139[i] = SGBrcharalt_p5_r10_g139[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p5_r30_g53[i] = SGBrcharalt_p5_r30_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p5_r30_g139[i] = SGBrcharalt_p5_r30_g139[i] * getCriticalB(zetas[i])**(-1)

    SGBthetaalt_p6_r10_g53[i] = SGBrcharalt_p6_r10_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p6_r10_g139[i] = SGBrcharalt_p6_r10_g139[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p6_r30_g53[i] = SGBrcharalt_p6_r30_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p6_r30_g139[i] = SGBrcharalt_p6_r30_g139[i] * getCriticalB(zetas[i])**(-1)

    SGBthetaalt_p7_r10_g53[i] = SGBrcharalt_p7_r10_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p7_r10_g139[i] = SGBrcharalt_p7_r10_g139[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p7_r30_g53[i] = SGBrcharalt_p7_r30_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p7_r30_g139[i] = SGBrcharalt_p7_r30_g139[i] * getCriticalB(zetas[i])**(-1)

    SGBthetaalt_p8_r10_g53[i] = SGBrcharalt_p8_r10_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p8_r10_g139[i] = SGBrcharalt_p8_r10_g139[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p8_r30_g53[i] = SGBrcharalt_p8_r30_g53[i] * getCriticalB(zetas[i])**(-1)
    SGBthetaalt_p8_r30_g139[i] = SGBrcharalt_p8_r30_g139[i] * getCriticalB(zetas[i])**(-1)

REZZtheta_p6_r10_g53 = np.zeros(2 * Nfiles_REZZ)
REZZtheta_p6_r10_g139 = np.zeros(2 * Nfiles_REZZ)
REZZtheta_p6_r30_g53 = np.zeros(2 * Nfiles_REZZ)
REZZtheta_p6_r30_g139 = np.zeros(2 * Nfiles_REZZ)

REZZtheta_p8_r10_g53 = np.zeros(2 * Nfiles_REZZ)
REZZtheta_p8_r10_g139 = np.zeros(2 * Nfiles_REZZ)
REZZtheta_p8_r30_g53 = np.zeros(2 * Nfiles_REZZ)
REZZtheta_p8_r30_g139 = np.zeros(2 * Nfiles_REZZ)

for i in range(0, 2 * Nfiles_REZZ):
    REZZtheta_p6_r10_g53[i] = REZZrchar_p6_r10_g53[i] * getCrtiticalBREZZ(a1s[i])**(-1)
    REZZtheta_p6_r10_g139[i] = REZZrchar_p6_r10_g139[i] * getCrtiticalBREZZ(a1s[i])**(-1)
    REZZtheta_p6_r30_g53[i] = REZZrchar_p6_r30_g53[i] * getCrtiticalBREZZ(a1s[i])**(-1)
    REZZtheta_p6_r30_g139[i] = REZZrchar_p6_r30_g139[i] * getCrtiticalBREZZ(a1s[i])**(-1)

    REZZtheta_p8_r10_g53[i] = REZZrchar_p8_r10_g53[i] * getCrtiticalBREZZ(a1s[i])**(-1)
    REZZtheta_p8_r10_g139[i] = REZZrchar_p8_r10_g139[i] * getCrtiticalBREZZ(a1s[i])**(-1)
    REZZtheta_p8_r30_g53[i] = REZZrchar_p8_r30_g53[i] * getCrtiticalBREZZ(a1s[i])**(-1)
    REZZtheta_p8_r30_g139[i] = REZZrchar_p8_r30_g139[i] * getCrtiticalBREZZ(a1s[i])**(-1)

lnGRtheta_r10_g53 = np.log(GRtheta_r10_g53)
lnGRtheta_r10_g139 = np.log(GRtheta_r10_g139)
lnGRtheta_r30_g53 = np.log(GRtheta_r30_g53)
lnGRtheta_r30_g139 = np.log(GRtheta_r30_g139)

lnGRthetaalt_r10_g53 = np.log(GRthetaalt_r10_g53)
lnGRthetaalt_r10_g139 = np.log(GRthetaalt_r10_g139)
lnGRthetaalt_r30_g53 = np.log(GRthetaalt_r30_g53)
lnGRthetaalt_r30_g139 = np.log(GRthetaalt_r30_g139)

lnSGBtheta_p5_r10_g53 = np.log(SGBtheta_p5_r10_g53)
lnSGBtheta_p5_r10_g139 = np.log(SGBtheta_p5_r10_g139)
lnSGBtheta_p5_r30_g53 = np.log(SGBtheta_p5_r30_g53)
lnSGBtheta_p5_r30_g139 = np.log(SGBtheta_p5_r30_g139)

lnSGBthetaalt_p5_r10_g53 = np.log(SGBthetaalt_p5_r10_g53)
lnSGBthetaalt_p5_r10_g139 = np.log(SGBthetaalt_p5_r10_g139)
lnSGBthetaalt_p5_r30_g53 = np.log(SGBthetaalt_p5_r30_g53)
lnSGBthetaalt_p5_r30_g139 = np.log(SGBthetaalt_p5_r30_g139)

lnSGBtheta_p6_r10_g53 = np.log(SGBtheta_p6_r10_g53)
lnSGBtheta_p6_r10_g139 = np.log(SGBtheta_p6_r10_g139)
lnSGBtheta_p6_r30_g53 = np.log(SGBtheta_p6_r30_g53)
lnSGBtheta_p6_r30_g139 = np.log(SGBtheta_p6_r30_g139)

lnSGBthetaalt_p6_r10_g53 = np.log(SGBthetaalt_p6_r10_g53)
lnSGBthetaalt_p6_r10_g139 = np.log(SGBthetaalt_p6_r10_g139)
lnSGBthetaalt_p6_r30_g53 = np.log(SGBthetaalt_p6_r30_g53)
lnSGBthetaalt_p6_r30_g139 = np.log(SGBthetaalt_p6_r30_g139)

lnSGBtheta_p7_r10_g53 = np.log(SGBtheta_p7_r10_g53)
lnSGBtheta_p7_r10_g139 = np.log(SGBtheta_p7_r10_g139)
lnSGBtheta_p7_r30_g53 = np.log(SGBtheta_p7_r30_g53)
lnSGBtheta_p7_r30_g139 = np.log(SGBtheta_p7_r30_g139)

lnSGBthetaalt_p7_r10_g53 = np.log(SGBthetaalt_p7_r10_g53)
lnSGBthetaalt_p7_r10_g139 = np.log(SGBthetaalt_p7_r10_g139)
lnSGBthetaalt_p7_r30_g53 = np.log(SGBthetaalt_p7_r30_g53)
lnSGBthetaalt_p7_r30_g139 = np.log(SGBthetaalt_p7_r30_g139)

lnSGBtheta_p8_r10_g53 = np.log(SGBtheta_p8_r10_g53)
lnSGBtheta_p8_r10_g139 = np.log(SGBtheta_p8_r10_g139)
lnSGBtheta_p8_r30_g53 = np.log(SGBtheta_p8_r30_g53)
lnSGBtheta_p8_r30_g139 = np.log(SGBtheta_p8_r30_g139)

lnSGBthetaalt_p8_r10_g53 = np.log(SGBthetaalt_p8_r10_g53)
lnSGBthetaalt_p8_r10_g139 = np.log(SGBthetaalt_p8_r10_g139)
lnSGBthetaalt_p8_r30_g53 = np.log(SGBthetaalt_p8_r30_g53)
lnSGBthetaalt_p8_r30_g139 = np.log(SGBthetaalt_p8_r30_g139)

lnREZZtheta_p6_r10_g53 = np.log(REZZtheta_p6_r10_g53)
lnREZZtheta_p6_r10_g139 = np.log(REZZtheta_p6_r10_g139)
lnREZZtheta_p6_r30_g53 = np.log(REZZtheta_p6_r30_g53)
lnREZZtheta_p6_r30_g139 = np.log(REZZtheta_p6_r30_g139)

lnREZZtheta_p8_r10_g53 = np.log(REZZtheta_p8_r10_g53)
lnREZZtheta_p8_r10_g139 = np.log(REZZtheta_p8_r10_g139)
lnREZZtheta_p8_r30_g53 = np.log(REZZtheta_p8_r30_g53)
lnREZZtheta_p8_r30_g139 = np.log(REZZtheta_p8_r30_g139)

intpowers = np.arange(5,8,0.1)
inta1s = np.arange(-0.2,0.2,0.05)

lnGRtheta_r10_g53_interp = np.interp(intpowers, powers, lnGRtheta_r10_g53)
lnGRtheta_r10_g139_interp = np.interp(intpowers, powers, lnGRtheta_r10_g139)
lnGRtheta_r30_g53_interp = np.interp(intpowers, powers, lnGRtheta_r30_g53)
lnGRtheta_r30_g139_interp = np.interp(intpowers, powers, lnGRtheta_r30_g139)

lnGRthetaalt_r10_g53_interp = np.interp(intpowers, powers, lnGRthetaalt_r10_g53)
lnGRthetaalt_r10_g139_interp = np.interp(intpowers, powers, lnGRthetaalt_r10_g139)
lnGRthetaalt_r30_g53_interp = np.interp(intpowers, powers, lnGRthetaalt_r30_g53)
lnGRthetaalt_r30_g139_interp = np.interp(intpowers, powers, lnGRthetaalt_r30_g139)

lnREZZtheta_p6_r10_g53_interp = np.interp(inta1s, a1s, lnREZZtheta_p6_r10_g53)
lnREZZtheta_p6_r10_g139_interp = np.interp(inta1s, a1s, lnREZZtheta_p6_r10_g139)
lnREZZtheta_p6_r30_g53_interp = np.interp(inta1s, a1s, lnREZZtheta_p6_r30_g53)
lnREZZtheta_p6_r30_g139_interp = np.interp(inta1s, a1s, lnREZZtheta_p6_r30_g139)

lnREZZtheta_p8_r10_g53_interp = np.interp(inta1s, a1s, lnREZZtheta_p8_r10_g53)
lnREZZtheta_p8_r10_g139_interp = np.interp(inta1s, a1s, lnREZZtheta_p8_r10_g139)
lnREZZtheta_p8_r30_g53_interp = np.interp(inta1s, a1s, lnREZZtheta_p8_r30_g53)
lnREZZtheta_p8_r30_g139_interp = np.interp(inta1s, a1s, lnREZZtheta_p8_r30_g139)

fig, ax = plt.subplots(1,3, sharey=True, figsize=(12,6))

fntsze = 22

ax[0].plot(intpowers, lnGRtheta_r10_g53_interp, color='k', label=r"$\tilde{R} = 10$, $\gamma = 5/3$")
ax[0].plot(intpowers, lnGRtheta_r10_g139_interp, color='b', label=r"$\tilde{R} = 10$, $\gamma = 13/9$")
ax[0].plot(intpowers, lnGRtheta_r30_g53_interp, color='r', label=r"$\tilde{R} \gg 1$, $\gamma = 5/3$")
ax[0].plot(intpowers, lnGRtheta_r30_g139_interp, color='g', label=r"$\tilde{R} \gg 1$, $\gamma = 13/9$")

ax[0].set_xlabel(r"$\alpha$", fontsize=fntsze)
ax[0].set_ylabel(r"$\ln \ \vartheta$", fontsize=fntsze)
ax[0].set_title("GR", fontsize=fntsze)

ax[1].plot(zetas, lnSGBtheta_p6_r10_g53, color='k', label=r"$\tilde{R} = 10$, $\gamma = 5/3$")
ax[1].plot(zetas, lnSGBtheta_p6_r10_g139, color='b', label=r"$\tilde{R} = 10$, $\gamma = 13/9$")
ax[1].plot(zetas, lnSGBtheta_p6_r30_g53, color='r', label=r"$\tilde{R} \gg 1$, $\gamma = 5/3$")
ax[1].plot(zetas, lnSGBtheta_p6_r30_g139, color='g', label=r"$\tilde{R} = \gg 1$, $\gamma = 13/9$")
ax[1].plot(zetas, lnSGBtheta_p8_r10_g53, color='k', linestyle="--", label=r"$\tilde{R} = 10$, $\gamma = 5/3$")
ax[1].plot(zetas, lnSGBtheta_p8_r10_g139, color='b', linestyle="--", label=r"$\tilde{R} = 10$, $\gamma = 13/9$")
ax[1].plot(zetas, lnSGBtheta_p8_r30_g53, color='r', linestyle="--", label=r"$\tilde{R} \gg 1$, $\gamma = 5/3$")
ax[1].plot(zetas, lnSGBtheta_p8_r30_g139, color='g', linestyle="--", label=r"$\tilde{R} \gg 1$, $\gamma = 13/9$")
ax[1].set_title(r"sGB", fontsize=fntsze)
ax[1].set_xlabel(r"$\zeta$", fontsize=fntsze)

ax[2].plot(inta1s, lnREZZtheta_p6_r10_g53_interp, color='k', label=r"$\tilde{R} = 10$, $\gamma = 5/3$")
ax[2].plot(inta1s, lnREZZtheta_p6_r10_g139_interp, color='b', label=r"$\tilde{R} = 10$, $\gamma = 13/9$")
ax[2].plot(inta1s, lnREZZtheta_p6_r30_g53_interp, color='r', label=r"$\tilde{R} \gg 1$, $\gamma = 5/3$")
ax[2].plot(inta1s, lnREZZtheta_p6_r30_g139_interp, color='g', label=r"$\tilde{R} = \gg 1$, $\gamma = 13/9$")
ax[2].plot(inta1s, lnREZZtheta_p8_r10_g53_interp, color='k', linestyle="--")
ax[2].plot(inta1s, lnREZZtheta_p8_r10_g139_interp, color='b', linestyle="--")
ax[2].plot(inta1s, lnREZZtheta_p8_r30_g53_interp, color='r', linestyle="--")
ax[2].plot(inta1s, lnREZZtheta_p8_r30_g139_interp, color='g', linestyle="--")
ax[2].set_title(r"RZ", fontsize=fntsze)
ax[2].set_xlabel(r"$a_{1}$", fontsize=fntsze)
ax[2].legend(fontsize=16, frameon=False)

for i in range(0,3):
    ax[i].tick_params(axis='x', labelsize=16)
    ax[i].tick_params(axis='y', labelsize=16)

fig.savefig("./pfig_theta_bondi.png",dpi=400)

plt.show()