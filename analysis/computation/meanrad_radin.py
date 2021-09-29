# Adam M. Bauer
# This code makes a figure that plots our \vartheta
# as a function of emissivity power, \alpha, in GR and 
# \vartheta as a function of deformation parameter in sGB
# and RZ for a fixed value of \alpha using the radial infall 
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
import random

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

cms = mpl.cm 

def getIntegration(integrand, xvals):
    tmp = 0

    for i in range(1, len(xvals)):
        tmp += 0.5 * (xvals[i] - xvals[i-1]) * (integrand[i] + integrand[i-1])
    
    return tmp

def getValue(func, integrand, xvals, tol, args):
    
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
    return np.sqrt(27) - 4397*zeta*(2430 * np.sqrt(3) )**(-1) + zeta**2 * (113374080 * np.sqrt(3))**(-1) *336780431 + zeta**3 * 114832399336847 * (np.sqrt(3) * 37192366944000)**(-1) + zeta**4 * 125183193573305833463 * (np.sqrt(3) * 52057412164177920000)**(-1) + zeta**5 * 47568239982564118029349 * (np.sqrt(3) * 126499511558952345600000)**(-1) - zeta**6 * 926458082627090815686668481877 * (np.sqrt(3) * 265588254508251628634112000000)**(-1)

def getCrtiticalBREZZ(a1):
    return np.sqrt(27) - 4 * a1 * (3 * np.sqrt(3))**(-1) + 64 * a1**3 * (729 * np.sqrt(3))**(-1) + 512 * a1**4 * (59049 * np.sqrt(3))**(-1) - 11264 * a1**5 * (531441 * np.sqrt(3))**(-1) 

PATH = "../../data/radin/"
GRfilenames = "GR/iData_radin_GR*_N2000.csv"
SGBfilenames5 = "SGB/iData_radin_SGB_z*_p" + str(5.0) + "*N2000.csv"
SGBfilenames55 = "SGB/iData_radin_SGB_z*_p" + str(5.5) + "*N2000.csv"
SGBfilenames6 = "SGB/iData_radin_SGB_z*_p" + str(6.0) + "*N2000.csv"
SGBfilenames65 = "SGB/iData_radin_SGB_z*_p" + str(6.5) + "*N2000.csv"
SGBfilenames7 = "SGB/iData_radin_SGB_z*_p" + str(7.0) + "*N2000.csv"
SGBfilenames75 = "SGB/iData_radin_SGB_z*_p" + str(7.5) + "*N2000.csv"
SGBfilenames8 = "SGB/iData_radin_SGB_z*_p" + str(8.0) + "*N2000.csv"

pREZZfilenames6 = "REZZ/iData_radin_REZZ_z0*_p" + str(6.0) + "*N2000.csv"
nREZZfilenames6 = "REZZ/iData_radin_REZZ_zneg*_p" + str(6.0) + "*N2000.csv"

pREZZfilenames8 = "REZZ/iData_radin_REZZ_z0*_p" + str(8.0) + "*N2000.csv"
nREZZfilenames8 = "REZZ/iData_radin_REZZ_zneg*_p" + str(8.0) + "*N2000.csv"

zetas = [0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175]
powers = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
a1s = [-0.2, -0.15, -0.1, -0.05, 0.05, 0.1, 0.15, 0.2]

GRfiles = glob.glob(PATH+GRfilenames)
SGBfiles5 = glob.glob(PATH+SGBfilenames5)
SGBfiles55 = glob.glob(PATH+SGBfilenames55)
SGBfiles6 = glob.glob(PATH+SGBfilenames6)
SGBfiles65 = glob.glob(PATH+SGBfilenames65)
SGBfiles7 = glob.glob(PATH+SGBfilenames7)
SGBfiles75 = glob.glob(PATH+SGBfilenames75)
SGBfiles8 = glob.glob(PATH+SGBfilenames8)
pREZZfiles6 = glob.glob(PATH+pREZZfilenames6)
nREZZfiles6 = glob.glob(PATH+nREZZfilenames6)
pREZZfiles8 = glob.glob(PATH+pREZZfilenames8)
nREZZfiles8 = glob.glob(PATH+nREZZfilenames8)

GRfiles.sort() # put low p vals first
SGBfiles5.sort() 
SGBfiles55.sort() 
SGBfiles6.sort() 
SGBfiles65.sort()
SGBfiles7.sort()
SGBfiles75.sort()
SGBfiles8.sort()
pREZZfiles6.sort()
nREZZfiles6.sort()
nREZZfiles6.reverse()
pREZZfiles8.sort()
nREZZfiles8.sort()
nREZZfiles8.reverse()

Nfiles_GR = len(GRfiles)

GR_data_b = np.zeros((Nfiles_GR, 2000))
GR_data_I = np.zeros((Nfiles_GR, 2000))

for i in range(0, Nfiles_GR):
    tmp_GR_data_b = np.genfromtxt(GRfiles[i], delimiter=",", usecols=0)
    tmp_GR_data_I = np.genfromtxt(GRfiles[i], delimiter=",", usecols=1)
    GR_data_b[i,:] = tmp_GR_data_b[:]
    GR_data_I[i,:] = tmp_GR_data_I[:]

Ninterp = 10000

GR_b = np.linspace(0,15,Ninterp)
GR_interp_I = np.zeros((Nfiles_GR, Ninterp))
for i in range(0, Nfiles_GR):
    GR_interp_I[i,:] = np.interp(GR_b, GR_data_b[i,:], GR_data_I[i,:])

Nfiles_SGB = len(SGBfiles5)
SGB5_data_b = np.zeros((Nfiles_SGB, 2000))
SGB5_data_I = np.zeros((Nfiles_SGB, 2000))
SGB55_data_b = np.zeros((Nfiles_SGB, 2000))
SGB55_data_I = np.zeros((Nfiles_SGB, 2000))
SGB6_data_b = np.zeros((Nfiles_SGB, 2000))
SGB6_data_I = np.zeros((Nfiles_SGB, 2000))
SGB65_data_b = np.zeros((Nfiles_SGB, 2000))
SGB65_data_I = np.zeros((Nfiles_SGB, 2000))
SGB7_data_b = np.zeros((Nfiles_SGB, 2000))
SGB7_data_I = np.zeros((Nfiles_SGB, 2000))
SGB75_data_b = np.zeros((Nfiles_SGB, 2000))
SGB75_data_I = np.zeros((Nfiles_SGB, 2000))
SGB8_data_b = np.zeros((Nfiles_SGB, 2000))
SGB8_data_I = np.zeros((Nfiles_SGB, 2000))

for i in range(0, Nfiles_SGB):
    tmp_SGB5_data_b = np.genfromtxt(SGBfiles5[i], delimiter=",", usecols=0)
    tmp_SGB5_data_I = np.genfromtxt(SGBfiles5[i], delimiter=",", usecols=1)
    SGB5_data_b[i,:] = tmp_SGB5_data_b[:]
    SGB5_data_I[i,:] = tmp_SGB5_data_I[:]

    tmp_SGB55_data_b = np.genfromtxt(SGBfiles55[i], delimiter=",", usecols=0)
    tmp_SGB55_data_I = np.genfromtxt(SGBfiles55[i], delimiter=",", usecols=1)
    SGB55_data_b[i,:] = tmp_SGB55_data_b[:]
    SGB55_data_I[i,:] = tmp_SGB55_data_I[:]

    tmp_SGB6_data_b = np.genfromtxt(SGBfiles6[i], delimiter=",", usecols=0)
    tmp_SGB6_data_I = np.genfromtxt(SGBfiles6[i], delimiter=",", usecols=1)
    SGB6_data_b[i,:] = tmp_SGB6_data_b[:]
    SGB6_data_I[i,:] = tmp_SGB6_data_I[:]

    tmp_SGB65_data_b = np.genfromtxt(SGBfiles65[i], delimiter=",", usecols=0)
    tmp_SGB65_data_I = np.genfromtxt(SGBfiles65[i], delimiter=",", usecols=1)
    SGB65_data_b[i,:] = tmp_SGB65_data_b[:]
    SGB65_data_I[i,:] = tmp_SGB65_data_I[:]

    tmp_SGB7_data_b = np.genfromtxt(SGBfiles7[i], delimiter=",", usecols=0)
    tmp_SGB7_data_I = np.genfromtxt(SGBfiles7[i], delimiter=",", usecols=1)
    SGB7_data_b[i,:] = tmp_SGB7_data_b[:]
    SGB7_data_I[i,:] = tmp_SGB7_data_I[:]

    tmp_SGB75_data_b = np.genfromtxt(SGBfiles75[i], delimiter=",", usecols=0)
    tmp_SGB75_data_I = np.genfromtxt(SGBfiles75[i], delimiter=",", usecols=1)
    SGB75_data_b[i,:] = tmp_SGB75_data_b[:]
    SGB75_data_I[i,:] = tmp_SGB75_data_I[:]

    tmp_SGB8_data_b = np.genfromtxt(SGBfiles8[i], delimiter=",", usecols=0)
    tmp_SGB8_data_I = np.genfromtxt(SGBfiles8[i], delimiter=",", usecols=1)
    SGB8_data_b[i,:] = tmp_SGB8_data_b[:]
    SGB8_data_I[i,:] = tmp_SGB8_data_I[:]

SGB5_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB55_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB6_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB65_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB7_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB75_interp_I = np.zeros((Nfiles_SGB, Ninterp))
SGB8_interp_I = np.zeros((Nfiles_SGB, Ninterp))


for i in range(0, Nfiles_SGB):
    SGB5_interp_I[i,:] = np.interp(GR_b, SGB5_data_b[i,:], SGB5_data_I[i,:])
    SGB55_interp_I[i,:] = np.interp(GR_b, SGB55_data_b[i,:], SGB55_data_I[i,:])
    SGB6_interp_I[i,:] = np.interp(GR_b, SGB6_data_b[i,:], SGB6_data_I[i,:])
    SGB65_interp_I[i,:] = np.interp(GR_b, SGB65_data_b[i,:], SGB65_data_I[i,:])
    SGB7_interp_I[i,:] = np.interp(GR_b, SGB7_data_b[i,:], SGB7_data_I[i,:])
    SGB75_interp_I[i,:] = np.interp(GR_b, SGB75_data_b[i,:], SGB75_data_I[i,:])
    SGB8_interp_I[i,:] = np.interp(GR_b, SGB8_data_b[i,:], SGB8_data_I[i,:])

Nfiles_pREZZ = len(pREZZfiles6)
pREZZ6_data_b = np.zeros((Nfiles_pREZZ, 2000))
pREZZ6_data_I = np.zeros((Nfiles_pREZZ, 2000))
nREZZ6_data_b = np.zeros((Nfiles_pREZZ, 2000))
nREZZ6_data_I = np.zeros((Nfiles_pREZZ, 2000))

pREZZ8_data_b = np.zeros((Nfiles_pREZZ, 2000))
pREZZ8_data_I = np.zeros((Nfiles_pREZZ, 2000))
nREZZ8_data_b = np.zeros((Nfiles_pREZZ, 2000))
nREZZ8_data_I = np.zeros((Nfiles_pREZZ, 2000))

for i in range(0, Nfiles_pREZZ):
    tmp_pREZZ6_data_b = np.genfromtxt(pREZZfiles6[i], delimiter=",", usecols=0)
    tmp_pREZZ6_data_I = np.genfromtxt(pREZZfiles6[i], delimiter=",", usecols=1)
    pREZZ6_data_b[i,:] = tmp_pREZZ6_data_b[:]
    pREZZ6_data_I[i,:] = tmp_pREZZ6_data_I[:]

    tmp_nREZZ6_data_b = np.genfromtxt(nREZZfiles6[i], delimiter=",", usecols=0)
    tmp_nREZZ6_data_I = np.genfromtxt(nREZZfiles6[i], delimiter=",", usecols=1)
    nREZZ6_data_b[i,:] = tmp_nREZZ6_data_b[:]
    nREZZ6_data_I[i,:] = tmp_nREZZ6_data_I[:]

    tmp_pREZZ8_data_b = np.genfromtxt(pREZZfiles8[i], delimiter=",", usecols=0)
    tmp_pREZZ8_data_I = np.genfromtxt(pREZZfiles8[i], delimiter=",", usecols=1)
    pREZZ8_data_b[i,:] = tmp_pREZZ8_data_b[:]
    pREZZ8_data_I[i,:] = tmp_pREZZ8_data_I[:]

    tmp_nREZZ8_data_b = np.genfromtxt(nREZZfiles8[i], delimiter=",", usecols=0)
    tmp_nREZZ8_data_I = np.genfromtxt(nREZZfiles8[i], delimiter=",", usecols=1)
    nREZZ8_data_b[i,:] = tmp_nREZZ8_data_b[:]
    nREZZ8_data_I[i,:] = tmp_nREZZ8_data_I[:]

pREZZ6_interp_I = np.zeros((Nfiles_pREZZ, Ninterp))
nREZZ6_interp_I = np.zeros((Nfiles_pREZZ, Ninterp))
pREZZ8_interp_I = np.zeros((Nfiles_pREZZ, Ninterp))
nREZZ8_interp_I = np.zeros((Nfiles_pREZZ, Ninterp))

for i in range(0, Nfiles_pREZZ):
    pREZZ6_interp_I[i,:] = np.interp(GR_b, pREZZ6_data_b[i,:], pREZZ6_data_I[i,:])
    nREZZ6_interp_I[i,:] = np.interp(GR_b, nREZZ6_data_b[i,:], nREZZ6_data_I[i,:])
    pREZZ8_interp_I[i,:] = np.interp(GR_b, pREZZ8_data_b[i,:], pREZZ8_data_I[i,:])
    nREZZ8_interp_I[i,:] = np.interp(GR_b, nREZZ8_data_b[i,:], nREZZ8_data_I[i,:])

GR_data_integrands = np.zeros((Nfiles_GR, Ninterp))
GR_data_integrandsalt = np.zeros((Nfiles_GR, Ninterp))

for i in range(0, Nfiles_GR):
    GR_data_integrands[i,:] = np.pi * 2 * GR_b**3 * GR_interp_I[i,:]
    GR_data_integrandsalt[i,:] = 2 * np.pi * GR_b * GR_interp_I[i,:]

SGB5_integrands = np.zeros((Nfiles_SGB, Ninterp))
SGB5_integrandsalt = np.zeros((Nfiles_SGB, Ninterp))

SGB55_integrands = np.zeros((Nfiles_SGB, Ninterp))
SGB55_integrandsalt = np.zeros((Nfiles_SGB, Ninterp))

SGB6_integrands = np.zeros((Nfiles_SGB, Ninterp))
SGB6_integrandsalt = np.zeros((Nfiles_SGB, Ninterp))

SGB65_integrands = np.zeros((Nfiles_SGB, Ninterp))
SGB65_integrandsalt = np.zeros((Nfiles_SGB, Ninterp))

SGB7_integrands = np.zeros((Nfiles_SGB, Ninterp))
SGB7_integrandsalt = np.zeros((Nfiles_SGB, Ninterp))

SGB75_integrands = np.zeros((Nfiles_SGB, Ninterp))
SGB75_integrandsalt = np.zeros((Nfiles_SGB, Ninterp))

SGB8_integrands = np.zeros((Nfiles_SGB, Ninterp))
SGB8_integrandsalt = np.zeros((Nfiles_SGB, Ninterp))


for i in range(0, Nfiles_SGB):
    SGB5_integrands[i,:] = np.pi * 2 * GR_b**3 * SGB5_interp_I[i,:]
    SGB5_integrandsalt[i,:] = 2 * np.pi *GR_b * SGB5_interp_I[i,:]

    SGB55_integrands[i,:] = np.pi * 2 * GR_b**3 * SGB55_interp_I[i,:]
    SGB55_integrandsalt[i,:] = 2 * np.pi * GR_b * SGB55_interp_I[i,:]

    SGB6_integrands[i,:] = np.pi * 2 * GR_b**3 * SGB6_interp_I[i,:]
    SGB6_integrandsalt[i,:] = 2 * np.pi *GR_b * SGB6_interp_I[i,:]

    SGB65_integrands[i,:] = np.pi * 2 * GR_b**3 * SGB65_interp_I[i,:]
    SGB65_integrandsalt[i,:] = 2 * np.pi *GR_b * SGB65_interp_I[i,:]

    SGB7_integrands[i,:] = np.pi * 2 * GR_b**3 * SGB7_interp_I[i,:]
    SGB7_integrandsalt[i,:] = 2 * np.pi *GR_b * SGB7_interp_I[i,:]

    SGB75_integrands[i,:] = np.pi * 2 * GR_b**3 * SGB75_interp_I[i,:]
    SGB75_integrandsalt[i,:] = 2 * np.pi *GR_b * SGB75_interp_I[i,:]

    SGB8_integrands[i,:] = np.pi * 2 * GR_b**3 * SGB8_interp_I[i,:]
    SGB8_integrandsalt[i,:] = 2 * np.pi *GR_b * SGB8_interp_I[i,:]

pREZZ6_integrands = np.zeros((Nfiles_pREZZ, Ninterp))
pREZZ6_integrandsalt = np.zeros((Nfiles_pREZZ, Ninterp))
nREZZ6_integrands = np.zeros((Nfiles_pREZZ, Ninterp))
nREZZ6_integrandsalt = np.zeros((Nfiles_pREZZ, Ninterp))
pREZZ8_integrands = np.zeros((Nfiles_pREZZ, Ninterp))
pREZZ8_integrandsalt = np.zeros((Nfiles_pREZZ, Ninterp))
nREZZ8_integrands = np.zeros((Nfiles_pREZZ, Ninterp))
nREZZ8_integrandsalt = np.zeros((Nfiles_pREZZ, Ninterp))

for i in range(0, Nfiles_pREZZ):
    pREZZ6_integrands[i,:] = np.pi * 2 * GR_b**3 * pREZZ6_interp_I[i,:]
    pREZZ6_integrandsalt[i,:] = 2 * np.pi * GR_b * pREZZ6_interp_I[i,:]

    nREZZ6_integrands[i,:] = np.pi * 2 * GR_b**3 * nREZZ6_interp_I[i,:]
    nREZZ6_integrandsalt[i,:] = 2 * np.pi * GR_b * nREZZ6_interp_I[i,:]
    
    pREZZ8_integrands[i,:] = np.pi * 2 * GR_b**3 * pREZZ8_interp_I[i,:]
    pREZZ8_integrandsalt[i,:] = 2 * np.pi * GR_b * pREZZ8_interp_I[i,:]
    
    nREZZ8_integrands[i,:] = np.pi * 2 * GR_b**3 * nREZZ8_interp_I[i,:]
    nREZZ8_integrandsalt[i,:] = 2 * np.pi * GR_b * nREZZ8_interp_I[i,:]

GR_data_rchar = np.zeros(Nfiles_GR)
GR_data_rcharalt = np.zeros(Nfiles_GR)

GRtol = 1 * 10**(-1)
tol = 1 * 10**(-1)

for i in range(0, Nfiles_GR):
    GR_data_rchar[i] = (getIntegration(GR_data_integrands[i,:], GR_b) * (getIntegration(GR_data_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    GR_data_rcharalt[i] = getValue(getIntegration, GR_data_integrandsalt[i,:], GR_b, GRtol, [getIntegration(GR_data_integrandsalt[i,:], GR_b)] )

SGB5_rchar = np.zeros(Nfiles_SGB)
SGB5_rcharalt = np.zeros(Nfiles_SGB)

SGB55_rchar = np.zeros(Nfiles_SGB)
SGB55_rcharalt = np.zeros(Nfiles_SGB)

SGB6_rchar = np.zeros(Nfiles_SGB)
SGB6_rcharalt = np.zeros(Nfiles_SGB)

SGB65_rchar = np.zeros(Nfiles_SGB)
SGB65_rcharalt = np.zeros(Nfiles_SGB)

SGB7_rchar = np.zeros(Nfiles_SGB)
SGB7_rcharalt = np.zeros(Nfiles_SGB)

SGB75_rchar = np.zeros(Nfiles_SGB)
SGB75_rcharalt = np.zeros(Nfiles_SGB)

SGB8_rchar = np.zeros(Nfiles_SGB)
SGB8_rcharalt = np.zeros(Nfiles_SGB)

for i in range(0, Nfiles_SGB):
    SGB5_rchar[i] = (getIntegration(SGB5_integrands[i,:], GR_b) * (getIntegration(SGB5_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    SGB5_rcharalt[i] = getValue(getIntegration, SGB5_integrandsalt[i,:], GR_b, tol, [getIntegration(SGB5_integrandsalt[i,:], GR_b)] )

    SGB55_rchar[i] = (getIntegration(SGB55_integrands[i,:], GR_b) * (getIntegration(SGB55_integrandsalt[i,:], GR_b))**(-1) )**(0.5)
    SGB55_rcharalt[i] = getValue(getIntegration, SGB55_integrandsalt[i,:], GR_b, tol, [getIntegration(SGB55_integrandsalt[i,:], GR_b)] )

    SGB6_rchar[i] = (getIntegration(SGB6_integrands[i,:], GR_b) * (getIntegration(SGB6_integrandsalt[i,:], GR_b))**(-1) )**(0.5)
    SGB6_rcharalt[i] = getValue(getIntegration, SGB6_integrandsalt[i,:], GR_b, tol, [getIntegration(SGB6_integrandsalt[i,:], GR_b)] )

    SGB65_rchar[i] = (getIntegration(SGB65_integrands[i,:], GR_b) * (getIntegration(SGB65_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    SGB65_rcharalt[i] = getValue(getIntegration, SGB65_integrandsalt[i,:], GR_b, tol, [getIntegration(SGB65_integrandsalt[i,:], GR_b)] )

    SGB7_rchar[i] = (getIntegration(SGB7_integrands[i,:], GR_b) * (getIntegration(SGB7_integrandsalt[i,:], GR_b))**(-1) )**(0.5)
    SGB7_rcharalt[i] = getValue(getIntegration, SGB7_integrandsalt[i,:], GR_b, tol, [getIntegration(SGB7_integrandsalt[i,:], GR_b)] )

    SGB75_rchar[i] = (getIntegration(SGB75_integrands[i,:], GR_b) * (getIntegration(SGB75_integrandsalt[i,:], GR_b))**(-1) )**(0.5)
    SGB75_rcharalt[i] = getValue(getIntegration, SGB75_integrandsalt[i,:], GR_b, tol, [getIntegration(SGB75_integrandsalt[i,:], GR_b)] )

    SGB8_rchar[i] = (getIntegration(SGB8_integrands[i,:], GR_b) * (getIntegration(SGB8_integrandsalt[i,:], GR_b))**(-1) )**(0.5)
    SGB8_rcharalt[i] = getValue(getIntegration, SGB8_integrandsalt[i,:], GR_b, tol, [getIntegration(SGB8_integrandsalt[i,:], GR_b)] )

pREZZ6_rchar = np.zeros(Nfiles_pREZZ)
pREZZ6_rcharalt = np.zeros(Nfiles_pREZZ)

nREZZ6_rchar = np.zeros(Nfiles_pREZZ)
nREZZ6_rcharalt = np.zeros(Nfiles_pREZZ)

pREZZ8_rchar = np.zeros(Nfiles_pREZZ)
pREZZ8_rcharalt = np.zeros(Nfiles_pREZZ)

nREZZ8_rchar = np.zeros(Nfiles_pREZZ)
nREZZ8_rcharalt = np.zeros(Nfiles_pREZZ)

for i in range(0, Nfiles_pREZZ):
    pREZZ6_rchar[i] = (getIntegration(pREZZ6_integrands[i,:], GR_b) * (getIntegration(pREZZ6_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    pREZZ6_rcharalt[i] = getValue(getIntegration, pREZZ6_integrandsalt[i,:], GR_b, tol, [getIntegration(pREZZ6_integrandsalt[i,:], GR_b)] )

    nREZZ6_rchar[i] = (getIntegration(nREZZ6_integrands[i,:], GR_b) * (getIntegration(nREZZ6_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    nREZZ6_rcharalt[i] = getValue(getIntegration, nREZZ6_integrandsalt[i,:], GR_b, tol, [getIntegration(nREZZ6_integrandsalt[i,:], GR_b)] )

    pREZZ8_rchar[i] = (getIntegration(pREZZ8_integrands[i,:], GR_b) * (getIntegration(pREZZ8_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    pREZZ8_rcharalt[i] = getValue(getIntegration, pREZZ8_integrandsalt[i,:], GR_b, tol, [getIntegration(pREZZ8_integrandsalt[i,:], GR_b)] )

    nREZZ8_rchar[i] = (getIntegration(nREZZ8_integrands[i,:], GR_b) * (getIntegration(nREZZ8_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    nREZZ8_rcharalt[i] = getValue(getIntegration, nREZZ8_integrandsalt[i,:], GR_b, tol, [getIntegration(nREZZ8_integrandsalt[i,:], GR_b)] )

GR_data_ratios = np.zeros(Nfiles_GR)
GR_data_ratiosalt = np.zeros(Nfiles_GR)

for i in range(0, Nfiles_GR):
    GR_data_ratios[i] = GR_data_rchar[i] * getCriticalB(0)**(-1) # GR
    GR_data_ratiosalt[i] = GR_data_rcharalt[i] * getCriticalB(0)**(-1) 

print(GR_data_ratios)

SGB5_ratios = np.zeros(Nfiles_SGB)
SGB5_ratiosalt = np.zeros(Nfiles_SGB)

SGB55_ratios = np.zeros(Nfiles_SGB)
SGB55_ratiosalt = np.zeros(Nfiles_SGB)

SGB6_ratios = np.zeros(Nfiles_SGB)
SGB6_ratiosalt = np.zeros(Nfiles_SGB)

SGB65_ratios = np.zeros(Nfiles_SGB)
SGB65_ratiosalt = np.zeros(Nfiles_SGB)

SGB7_ratios = np.zeros(Nfiles_SGB)
SGB7_ratiosalt = np.zeros(Nfiles_SGB)

SGB75_ratios = np.zeros(Nfiles_SGB)
SGB75_ratiosalt = np.zeros(Nfiles_SGB)

SGB8_ratios = np.zeros(Nfiles_SGB)
SGB8_ratiosalt = np.zeros(Nfiles_SGB)

for i in range(0, Nfiles_SGB):
    SGB5_ratios[i] = SGB5_rchar[i] * getCriticalB(zetas[i])**(-1) 
    SGB5_ratiosalt[i] = SGB5_rcharalt[i] * getCriticalB(zetas[i])**(-1) 

    SGB55_ratios[i] = SGB55_rchar[i] * getCriticalB(zetas[i])**(-1) 
    SGB55_ratiosalt[i] = SGB55_rcharalt[i] * getCriticalB(zetas[i])**(-1) 

    SGB6_ratios[i] = SGB6_rchar[i] * getCriticalB(zetas[i])**(-1) 
    SGB6_ratiosalt[i] = SGB6_rcharalt[i] * getCriticalB(zetas[i])**(-1) 

    SGB65_ratios[i] = SGB65_rchar[i] * getCriticalB(zetas[i])**(-1) 
    SGB65_ratiosalt[i] = SGB65_rcharalt[i] * getCriticalB(zetas[i])**(-1) 

    SGB7_ratios[i] = SGB7_rchar[i] * getCriticalB(zetas[i])**(-1) 
    SGB7_ratiosalt[i] = SGB7_rcharalt[i] * getCriticalB(zetas[i])**(-1) 

    SGB75_ratios[i] = SGB75_rchar[i] * getCriticalB(zetas[i])**(-1) 
    SGB75_ratiosalt[i] = SGB75_rcharalt[i] * getCriticalB(zetas[i])**(-1) 

    SGB8_ratios[i] = SGB8_rchar[i] * getCriticalB(zetas[i])**(-1) 
    SGB8_ratiosalt[i] = SGB8_rcharalt[i] * getCriticalB(zetas[i])**(-1) 

REZZ6_rchars = np.hstack((nREZZ6_rchar, pREZZ6_rchar))
REZZ8_rchars = np.hstack((nREZZ8_rchar, pREZZ8_rchar))

REZZ6_ratios = np.zeros(Nfiles_pREZZ * 2)
REZZ8_ratios = np.zeros(Nfiles_pREZZ * 2)

for i in range(0, 2 * Nfiles_pREZZ):
    REZZ6_ratios[i] = REZZ6_rchars[i] * getCrtiticalBREZZ(a1s[i])**(-1)
    REZZ8_ratios[i] = REZZ8_rchars[i] * getCrtiticalBREZZ(a1s[i])**(-1)

"""
reg5 = linregress(zetas, SGB5_ratios)
reg55 = linregress(zetas, SGB55_ratios)
reg6 = linregress(zetas, SGB6_rchar* GR_data_rchar[2]**(-1))
reg65 = linregress(zetas, SGB6_ratios* GR_data_ratios[2]**(-1))
reg7 = linregress(zetas, SGB7_ratios)
reg75 = linregress(zetas, SGB75_ratios)
reg8 = linregress(zetas, SGB8_ratios)

regcoeffs = [reg5.slope, reg55.slope, reg6.slope, reg65.slope, reg7.slope, reg75.slope, reg8.slope]

finalreg = linregress(powers, regcoeffs)
print(GR_data_rchar[2], reg6.slope)

print(reg6.slope)
print(- 4397*(2430 * np.sqrt(3) * np.sqrt(27))**(-1))
print(reg65.slope)


#fig, ax = plt.subplots(1)

#ax.plot(powers, regcoeffs)
#ax.set_xlabel(r"Power of $j_{\nu}$", fontsize=16)
#ax.set_ylabel(r"Values of $A$ such that $\vartheta = A \zeta$", fontsize=16)

#fig.savefig("regresscoeff_vs_power.pdf")

"""

GR_data_lnratios = np.zeros(Nfiles_GR)
GR_data_lnratiosalt = np.zeros(Nfiles_GR)
lnpowers = np.zeros(Nfiles_GR)

for i in range(0, Nfiles_GR):
    GR_data_lnratios[i] = np.log(GR_data_ratios[i])
    GR_data_lnratiosalt[i] = np.log(GR_data_ratiosalt[i])
    lnpowers[i] = np.log(powers[i])

SGB5_lnratios = np.zeros(Nfiles_SGB)
SGB5_lnratiosalt = np.zeros(Nfiles_SGB)

SGB6_lnratios = np.zeros(Nfiles_SGB)
SGB6_lnratiosalt = np.zeros(Nfiles_SGB)

SGB8_lnratios = np.zeros(Nfiles_SGB)
SGB8_lnratiosalt = np.zeros(Nfiles_SGB)

lnzetas = np.zeros(Nfiles_SGB)

for i in range(0, Nfiles_SGB):
    SGB5_lnratios[i] = np.log(SGB5_ratios[i])
    SGB5_lnratiosalt[i] = np.log(SGB5_ratiosalt[i])

    SGB6_lnratios[i] = np.log(SGB6_ratios[i])
    SGB6_lnratiosalt[i] = np.log(SGB6_ratiosalt[i])

    SGB8_lnratios[i] = np.log(SGB8_ratios[i])
    SGB8_lnratiosalt[i] = np.log(SGB8_ratiosalt[i])

    lnzetas[i] = np.log(zetas[i])

REZZ6_lnratios = np.log(REZZ6_ratios)
REZZ8_lnratios = np.log(REZZ8_ratios)

intpowers = np.arange(5,8,0.1)
intzetas = np.arange(0.025,0.175,0.05)
inta1s = np.arange(-0.2, 0.2, 0.05)

GR_data_lnratios_interp = np.interp(intpowers, powers, GR_data_lnratios)
GR_data_lnratiosalt_interp = np.interp(intpowers, powers, GR_data_lnratiosalt)

SGB5_lnratios_interp = np.interp(intzetas, zetas, SGB5_lnratios)
SGB5_lnratiosalt_interp = np.interp(intzetas, zetas, SGB5_lnratiosalt)

SGB6_lnratios_interp = np.interp(intzetas, zetas, SGB6_lnratios)
SGB6_lnratiosalt_interp = np.interp(intzetas, zetas, SGB6_lnratiosalt)

SGB8_lnratios_interp = np.interp(intzetas, zetas, SGB8_lnratios)
SGB8_lnratiosalt_interp = np.interp(intzetas, zetas, SGB8_lnratiosalt)

REZZ6_lnratios_interp = np.interp(inta1s, a1s, REZZ6_lnratios)
REZZ8_lnratios_interp = np.interp(inta1s, a1s, REZZ8_lnratios)


linregressREZZ6 = linregress(inta1s, REZZ6_lnratios_interp)
linregressREZZ8 = linregress(inta1s, REZZ8_lnratios_interp)
print(linregressREZZ6.slope, linregressREZZ8.slope)
"""
regr = linregress(intpowers, GR_data_lnratios_interp)
reg3 = linregress(intzetas, SGB5_lnratios_interp)
reg45 = linregress(intzetas, SGB6_lnratios_interp)
reg6 = linregress(intzetas, SGB8_lnratios_interp)

print(regr.slope, regr.stderr, reg3.slope, reg3.stderr, reg45.slope, reg45.stderr, reg6.slope, reg6.stderr)
"""

fig, ax = plt.subplots(1,3, sharey=True, figsize=(12,6))

fntsze = 22

ax[0].plot(intpowers, GR_data_lnratios_interp, color='k', label="2nd Moment")
#ax[0].plot(intpowers, GR_data_lnratiosalt_interp, color='r', label=r"Medain $b$")
ax[0].set_xlabel(r"$\alpha$", fontsize=fntsze)
ax[0].set_ylabel(r"$\ln \ \vartheta$", fontsize=fntsze)
ax[0].set_title("GR", fontsize=fntsze)

ax[1].plot(intzetas, SGB6_lnratios_interp, color='r', linestyle='--', label="2nd Moment")
ax[1].plot(intzetas, SGB8_lnratios_interp, color='b')
#ax[2-1].plot(intzetas, SGB6_lnratiosalt_interp, color='r', label=r"Medain $b$")
ax[1].set_xlabel(r"$ \zeta $", fontsize=fntsze)
#ax[2].set_ylabel(r"$\ln \ \vartheta$", fontsize=16)
ax[1].set_title(r"sGB", fontsize=fntsze)

"""
ax[2].plot(intzetas, SGB8_lnratios_interp, color='k', label="2nd Moment")
#ax[3-1].plot(intzetas, SGB8_lnratiosalt_interp, color='r', label=r"Medain $b$")
ax[2].set_xlabel(r"$ \zeta $", fontsize=fntsze)
#ax[3].set_ylabel(r"$\ln \ \vartheta$", fontsize=16)
ax[2].set_title(r"$\alpha = 8$", fontsize=fntsze)
#ax[3-1].legend(fontsize=18)
"""

ax[2].plot(inta1s, REZZ6_lnratios_interp, color='r', linestyle='--', label=r"$\alpha = 6$")
ax[2].plot(inta1s, REZZ8_lnratios_interp, color='b', label=r"$\alpha = 8$")
ax[2].set_title("RZ", fontsize=fntsze)
ax[2].set_xlabel(r"$a_{1}$", fontsize=fntsze)
ax[2].legend(fontsize=20, frameon=False)

for i in range(0,3):
    ax[i].tick_params(axis='x', labelsize=16)
    ax[i].tick_params(axis='y', labelsize=16)

fig.savefig("./pfig_theta_radin.png", dpi=400)

"""
# second figure for Nico 

rel3ratios = abs(SGB3_ratios - GR_data_ratios[0]) * (GR_data_ratios[0])**(-1)
rel45ratios = abs(SGB45_ratios - GR_data_ratios[3]) * (GR_data_ratios[3])**(-1)
rel6ratios = abs(SGB6_ratios - GR_data_ratios[6]) * (GR_data_ratios[6])**(-1)

rel3ln = np.zeros(len(rel3ratios))
rel45ln = np.zeros(len(rel45ratios))
rel6ln = np.zeros(len(rel6ratios))

for i in range(0, len(rel3ratios)):
    rel3ln[i] = np.log(rel3ratios[i])
    rel45ln[i] = np.log(rel45ratios[i])
    rel6ln[i] = np.log(rel6ratios[i])

#print(len(rel3ratios), len(rel3ln), len(zetas))

rel3regress = linregress(lnzetas, rel3ln)
rel45regress = linregress(lnzetas, rel45ln)
rel6regress = linregress(lnzetas, rel6ln)

print(rel3regress.slope, rel3regress.stderr, rel45regress.slope, rel45regress.stderr, rel6regress.slope, rel6regress.stderr)


#print(rel3ratios)

fig, ax = plt.subplots(1,3, sharey=True, figsize=(12,6))

ax[0].plot(zetas, rel3ratios, 'b.')
ax[0].set_xlabel(r"$\zeta$", fontsize=16)
ax[0].set_ylabel(r"$|\vartheta - \vartheta_{GR}|/\vartheta_{GR}$", fontsize=16)
ax[0].text(0.05,-0.004, r"$\alpha = 5$", fontsize=16)
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_title(r"$\alpha = 5$", fontsize=16)
ax[1].plot(zetas, rel45ratios, 'b.')
ax[1].set_xlabel(r"$\zeta$", fontsize=16)
ax[1].text(0.05,-0.004, r"$\alpha = 6.5$", fontsize=16)
ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_title(r"$\alpha = 6.5$", fontsize=16)
ax[2].plot(zetas, rel6ratios, 'b.')
ax[2].set_xlabel(r"$\zeta$", fontsize=16)
ax[2].text(0.05,-0.004, r"$\alpha = 8$", fontsize=16)
ax[2].set_yscale('log')
ax[2].set_xscale('log')
ax[2].set_title(r"$\alpha = 8$", fontsize=16)
fig.savefig("nicocheck_newzetarange.pdf")


plt.show()
"""

"""
GRtol = 1 * 10**(-3)
tol = 1 * 10**(-3)
partials = np.zeros((Nfiles_GR, Ninterp))
ratios = np.zeros((Nfiles_GR,Ninterp))

for i in range(0, Nfiles_GR):
    print("file" + str(i))
    full = (getIntegration(GR_data_integrands[i,:], GR_b) * (getIntegration(GR_data_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    for j in range(100, Ninterp):
        print("num" + str(j))
        partials[i,j] = (getIntegration(GR_data_integrands[i,0:j], GR_b[0:j]) * (getIntegration(GR_data_integrandsalt[i,0:j], GR_b[0:j]))**(-1))**(0.5)
        ratios[i,j] = partials[i,j] * full**(-1)

sgbpartials = np.zeros((Nfiles_SGB, Ninterp))
sgbratios = np.zeros((Nfiles_SGB,Ninterp))

for i in range(0, Nfiles_SGB):
    print("file" + str(i))
    sgbfull = (getIntegration(SGB6_integrands[i,:], GR_b) * (getIntegration(SGB6_integrandsalt[i,:], GR_b))**(-1))**(0.5)
    for j in range(100, Ninterp):
        print("num" + str(j))
        sgbpartials[i,j] = (getIntegration(SGB6_integrands[i,0:j], GR_b[0:j]) * (getIntegration(SGB6_integrandsalt[i,0:j], GR_b[0:j]))**(-1))**(0.5)
        sgbratios[i,j] = sgbpartials[i,j] * sgbfull**(-1)

fig, ax = plt.subplots(1,2, sharey=True)

for i in range(0,Nfiles_GR):
    label1=str(powers[i])
    color1 = (random.random(), random.random(), random.random())
    ax[0].plot(GR_b[100:], ratios[i,100:], color=color1, label=label1)

for i in range(0,Nfiles_SGB):
    label2=str(zetas[i])
    color2 = (random.random(), random.random(), random.random())
    ax[1].plot(GR_b[100:], sgbratios[i,100:], color=color2, label=label2)


ax[0].legend()
ax[1].legend()
ax[0].set_xlabel(r'upper critical impact parameter, $\sigma$')
ax[1].set_xlabel(r'upper critical impact parameter, $\sigma$')
ax[1].set_title(r"$\alpha = 6$")
ax[0].set_ylabel(r'$r_{char}$')

fig.savefig('percentfull.pdf')
plt.show()

"""
fig, ax = plt.subplots(2,2, figsize=(12,6))
ax[0,0].plot(zetas, SGB6_rchar * GR_data_rchar[2]**(-1))
ax[0,0].set_ylabel(r'$r_{char}/r_{char,GR}$', fontsize=16)
ax[0,0].set_xlabel(r'$\zeta$', fontsize=16)

crits = []
for i in range(0, len(zetas)):
    crits.append(getCriticalB(zetas[i]) * np.sqrt(27)**(-1) )

ax[0,1].plot(zetas, crits)
ax[0,1].set_ylabel(r"$b_{crit}/b_{crit,GR}$", fontsize=16)
ax[0,1].set_xlabel(r"$\zeta$", fontsize=16)

ax[1,0].plot(zetas, SGB6_ratios * GR_data_ratios[2]**(-1))
ax[1,0].set_ylabel(r"$\vartheta/\vartheta_{GR}$", fontsize=16)
ax[1,0].set_xlabel(r"$\zeta$", fontsize=16)

ax[1,1].plot(zetas, SGB6_rchar * GR_data_rchar[2]**(-1), color='r', label=r"$r_{char}/r_{char,GR}$")
ax[1,1].plot(zetas, crits, color='k', label=r"$b_{crit}/b_{crit,GR}$")
ax[1,1].set_xlabel(r"$\zeta$", fontsize=16)
ax[1,1].legend()

fig.tight_layout()

plt.show()