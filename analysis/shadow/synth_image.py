"""
Adam M. Bauer
adammb4@illinois.edu

This code is meant to make a "EHT style" image of the 
black hole shadow. The trick is the following: in MasterIntensityCalculator.py,
we noted using a check when the ray crossed x = 23 MG/c^2. This was used to find
the y position of the ray at that point, and we also noted its intensity at the 
camera. This creates a sort of "screen" at x = 23 MG/c^2. The curvature effects at
those radii are not strong, and so we can approximate the spacetime as flat,
thus giving us a CCD like image.

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
import matplotlib.pyplot as plt 
from matplotlib import font_manager
from matplotlib import ticker
import matplotlib as mpl

import matplotlib.colors as clrs 

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

cms = mpl.cm 

# GR radin p = 6
PATHrad = '../../data/radin/'
b_radin_gr = np.genfromtxt(PATHrad + "GR/iData_radin_GR_p6.0_N2000.csv", delimiter=',', usecols=0)
intensity_radin_gr = np.genfromtxt(PATHrad + "GR/iData_radin_GR_p6.0_N2000.csv", delimiter=',', usecols=1)

# sgb radin p = 6; zeta = 0.175
b_radin_sgb = np.genfromtxt(PATHrad + "SGB/iData_radin_SGB_z0.175_p6.0_N2000.csv", delimiter=',', usecols=0)
intensity_radin_sgb = np.genfromtxt(PATHrad + "SGB/iData_radin_SGB_z0.175_p6.0_N2000.csv", delimiter=',', usecols=1)

# rezz radin p = 6; a_1 = 0.2 
b_radin_rezz = np.genfromtxt(PATHrad + "REZZ/iData_radin_REZZ_z0.20_p6.0_N2000.csv", delimiter=',', usecols=0)
intensity_radin_rezz = np.genfromtxt(PATHrad + "REZZ/iData_radin_REZZ_z0.20_p6.0_N2000.csv", delimiter=',', usecols=1)

# GR bondi p = 6; rtilde = 10; gamma = 5/3
PATHbondi = '../../data/bondi/'
b_bondi_gr = np.genfromtxt(PATHbondi + "GR/iData_bondi_GR_p6.0_rcrit10.0_g53_N2000.csv", delimiter=',', usecols=0)
intensity_bondi_gr = np.genfromtxt(PATHbondi + "GR/iData_bondi_GR_p6.0_rcrit10.0_g53_N2000.csv", delimiter=',', usecols=1)

# sgb bondi p = 6; rtilde = 10; gamma = 5/3; zeta = 0.175
b_bondi_sgb = np.genfromtxt(PATHbondi + "SGB/iData_bondi_SGB_z0.175_p6.0_rcrit10.0_g53_N2000.csv", delimiter=',', usecols=0)
intensity_bondi_sgb = np.genfromtxt(PATHbondi + "SGB/iData_bondi_SGB_z0.175_p6.0_rcrit10.0_g53_N2000.csv", delimiter=',', usecols=1)

# rezz bondi p = 6; rtilde = 10; gamma = 5/3; zeta = 0.175
b_bondi_rezz = np.genfromtxt(PATHbondi + "REZZ/iData_bondi_REZZ_z0.20_p6.0_rcrit10.0_g53_N2000.csv", delimiter=',', usecols=0)
intensity_bondi_rezz = np.genfromtxt(PATHbondi + "REZZ/iData_bondi_REZZ_z0.20_p6.0_rcrit10.0_g53_N2000.csv", delimiter=',', usecols=1)

fig, ax = plt.subplots(2, 3, figsize=(13,9))

colorscale = 8 * 10**4

#fig.suptitle("Black Hole Shadow", fontsize=22)

#### GR CASE

ax[0,0].set_title(r"Radial Infall, GR", fontsize=17)

radin_gr_normfactor = max(np.log(intensity_radin_gr/min(intensity_radin_gr)))
for i in range(0, len(b_radin_gr)):

    order = len(b_radin_gr) - i
    
    tmpratio = np.log(intensity_radin_gr[i]/min(intensity_radin_gr))
    tmpcirc = plt.Circle((0,0), b_radin_gr[i], color=mpl.cm.afmhot(tmpratio/radin_gr_normfactor), zorder=order, fill=True)
    ax[0,0].add_patch(tmpcirc)

ax[0,0].set_xlim((-10,10))
ax[0,0].set_ylim((-10,10))

# for major ticks
ax[0,0].set_xticks([-10,-5,0,5,10])
# for minor ticks
ax[0,0].set_xticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[0,0].set_yticks([-10,-5,0,5,10])
# for minor ticks
ax[0,0].set_yticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[0,0].set_aspect('equal')

#### 

#### sgb radin PLOT

ax[0,1].set_title(r"Radial Infall, sGB", fontsize=17)

radin_sgb_normfactor = max(np.log(intensity_radin_sgb/min(intensity_radin_sgb)))
for i in range(0, len(b_radin_sgb)):

    tmpratio = np.log(intensity_radin_sgb[i]/min(intensity_radin_sgb))
    color = cms.afmhot(tmpratio/radin_sgb_normfactor)
    order = len(b_radin_sgb) - i

    tmpcirc = plt.Circle((0,0), b_radin_sgb[i], color=color, zorder=order, fill=True)
    ax[0,1].add_patch(tmpcirc)

ax[0,1].set_xlim((-10,10))
ax[0,1].set_ylim((-10,10))

# for major ticks
ax[0,1].set_xticks([-10,-5,0,5,10])
# for minor ticks
ax[0,1].set_xticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[0,1].set_yticks([-10,-5,0,5,10])
# for minor ticks
ax[0,1].set_yticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[0,1].set_aspect('equal')

####

#### rezz radin PLOT

ax[0,2].set_title(r"Radial Infall, RZ", fontsize=17)

radin_rz_normfactor = max(np.log(intensity_radin_rezz/min(intensity_radin_rezz)))
for i in range(0, len(b_radin_rezz)):

    tmpratio = np.log(intensity_radin_rezz[i]/min(intensity_radin_rezz))
    color = cms.afmhot(tmpratio/radin_rz_normfactor)
    order = len(b_radin_rezz) - i

    tmpcirc = plt.Circle((0,0), b_radin_rezz[i], color=color, zorder=order, fill=True)
    ax[0,2].add_patch(tmpcirc)

ax[0,2].set_xlim((-10,10))
ax[0,2].set_ylim((-10,10))

# for major ticks
ax[0,2].set_xticks([-10,-5,0,5,10])
# for minor ticks
ax[0,2].set_xticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[0,2].set_yticks([-10,-5,0,5,10])
# for minor ticks
ax[0,2].set_yticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[0,2].set_aspect('equal')

####

#### GR BONDI PLOT

ax[1,0].set_title(r"Bondi Accretion, GR", fontsize=17)

bondi_gr_normfactor = max(np.log(intensity_bondi_gr/min(intensity_bondi_gr)))
# BONDI PLOT
for i in range(0, len(b_bondi_gr)):

    order = len(b_bondi_gr) - i

    tmpratio = np.log(intensity_bondi_gr[i]/min(intensity_bondi_gr))
    tmpcirc = plt.Circle((0,0), b_bondi_gr[i], color=mpl.cm.afmhot(tmpratio/bondi_gr_normfactor), zorder=order, fill=True)
    ax[1,0].add_patch(tmpcirc)

ax[1,0].set_xlim((-10,10))
ax[1,0].set_ylim((-10,10))

# for major ticks
ax[1,0].set_xticks([-10,-5,0,5,10])
# for minor ticks
ax[1,0].set_xticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[1,0].set_yticks([-10,-5,0,5,10])
# for minor ticks
ax[1,0].set_yticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[1,0].set_aspect('equal')

#### 

#### ZETA = 0.2 PLOT

ax[1,1].set_title(r"Bondi Accretion, sGB", fontsize=17)

bondi_sgb_normfactor = max(np.log(intensity_bondi_sgb/min(intensity_bondi_sgb)))
for i in range(0, len(b_bondi_sgb)):

    order = len(b_bondi_sgb) - i

    tmpratio = np.log(intensity_bondi_sgb[i]/min(intensity_bondi_sgb))
    tmpcirc = plt.Circle((0,0), b_bondi_sgb[i], color=mpl.cm.afmhot(tmpratio/bondi_sgb_normfactor), zorder=order, fill=True)
    ax[1,1].add_patch(tmpcirc)

ax[1,1].set_xlim((-10,10))
ax[1,1].set_ylim((-10,10))

# for major ticks
ax[1,1].set_xticks([-10,-5,0,5,10])
# for minor ticks
ax[1,1].set_xticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[1,1].set_yticks([-10,-5,0,5,10])
# for minor ticks
ax[1,1].set_yticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[1,1].set_aspect('equal')

for i in range(0,2):
    for j in range(0,3):
        ax[i,j].tick_params(axis='both', labelsize=15)

#### ZETA = 0.2 PLOT

ax[1,2].set_title(r"Bondi Accretion, RZ", fontsize=17)
bondi_rezz_normfactor = max(np.log(intensity_bondi_rezz/min(intensity_bondi_rezz)))

for i in range(0, len(b_bondi_rezz)):

    order = len(b_bondi_rezz) - i
    tmpratio = np.log(intensity_bondi_rezz[i]/min(intensity_bondi_rezz))

    #tmpcirc = plt.Circle((0,0), b_bondi_rezz[i], color=mpl.cm.afmhot(max(intensity_bondi_rezz)**(-1)*intensity_bondi_rezz[i]), zorder=order, fill=True)
    tmpcirc = plt.Circle((0,0), b_bondi_rezz[i], color=mpl.cm.afmhot(tmpratio/bondi_rezz_normfactor), zorder=order, fill=True)
    ax[1,2].add_patch(tmpcirc)

ax[1,2].set_xlim((-10,10))
ax[1,2].set_ylim((-10,10))

# for major ticks
ax[1,2].set_xticks([-10,-5,0,5,10])
# for minor ticks
ax[1,2].set_xticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[1,2].set_yticks([-10,-5,0,5,10])
# for minor ticks
ax[1,2].set_yticks([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10], minor=True)

ax[1,2].set_aspect('equal')

fig.tight_layout(pad=5)

fig.savefig('./pfig_EHTimages.pdf')

plt.show()