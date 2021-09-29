"""
Adam M. Bauer
adammb4@illinos.edu

This code has just the equation of motion for geodesics
that are sent from the camera towards the supermassive
black hole. 

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

def getGeodesicEvolution(a, var, expressions):

    #if a == 0:
        #print("In we go...")

    # compute christoffels FROM MATHEMATICA BB
    # these are all in upper,lower,lower notation, as i'm pretty sure this is all i will need in this problem.
    # note var[4] = t, var[5] = r, var[6] = theta, var[7] = phi.

    kti = var[0]
    kri = var[1]
    kthi = var[2]
    kphi = var[3]
    #ti = var[4]
    ri = var[5]
    #thi = var[6]
    #phi = var[7]

    # expressions[i] referes to which entry in equation_list (initialized in MasterIntensityCalculator.py) corresponds to that equation

    chris_ttr_fcn = expressions[8] 
    chris_rtt_fcn = expressions[9]
    chris_rrr_fcn = expressions[10]
    chris_rthth_fcn = expressions[11]
    chris_rpp_fcn = expressions[12]
    #chris_thrth_fcn = expressions[13]
    chris_prp_fcn = expressions[14]
    chris_pthp = 0

    # evaluate the above functions at the point we're currently at, ri

    chris_ttr = chris_ttr_fcn(ri)
    chris_rtt = chris_rtt_fcn(ri)
    chris_rrr = chris_rrr_fcn(ri)
    chris_rthth = chris_rthth_fcn(ri)
    chris_rpp = chris_rpp_fcn(ri)
    #chris_thrth = chris_thrth_fcn(ri)
    chris_prp = chris_prp_fcn(ri)

    # geodesic equation ripped straight out of Sean Carrolls' masterclass of a book

    dkt = - 2 * chris_ttr * kri * kti 
    dkr = - chris_rpp * kphi**2 - chris_rthth * kthi**2 - chris_rrr * kri**2 - chris_rtt * kti**2
    dkth = 0
    dkp = - 2*chris_prp * kri * kphi - 2*chris_pthp * kphi * kthi
    dt = kti
    dr = kri
    dth = kthi
    dp = kphi

    return [dkt, dkr, dkth, dkp, dt, dr, dth, dp]