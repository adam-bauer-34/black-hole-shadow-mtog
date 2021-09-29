"""
Adam M. Bauer
adammb4@illinois.edu

This is the equation of motion, solved in Geodesic.py, for the ray being traced back to
the camera. It also has an intensity profile. 

Args: 
- a (dummy arg for integrator)
- var (variables being integrated)
- expressions: the list of equations that are used in the equation of motion for the goedesic and for intensity
- zeta: modified gravity parameter
- accword: "accretion keyword", tells the code which model to use. options are
    - "stat": stationary emitting gas model
    - "radin": radial infall model, zero temperature limit of bondi flow
    - "bondi": bondi accretion model
- power: the power of the emissivity, i.e., j_nu ~ r^(-power)

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

def getGeodesicEvolutionPlus(a, var, expressions, zeta, accword, power, rcrit, gamma):

    #if a == 0:
        #print("And out we move...")

    # compute christoffels FROM MATHEMATICA YUH
    # these are all in upper,lower,lower notation, as i'm pretty sure this is all i will need in this problem.

    kti = var[0]
    kri = var[1]
    kthi = var[2]
    kphi = var[3]
    #ti = var[4]
    ri = var[5]
    #thi = var[6]
    #phi = var[7]

    # again, import functions from equation_list

    chris_ttr_fcn = expressions[8]
    chris_rtt_fcn = expressions[9]
    chris_rrr_fcn = expressions[10]
    chris_rthth_fcn = expressions[11]
    chris_rpp_fcn = expressions[12]
    #chris_thrth_fcn = expressions[13]
    chris_prp_fcn = expressions[14]
    chris_pthp = 0

    # eval at ri 

    chris_ttr = chris_ttr_fcn(ri)
    chris_rtt = chris_rtt_fcn(ri)
    chris_rrr = chris_rrr_fcn(ri)
    chris_rthth = chris_rthth_fcn(ri)
    chris_rpp = chris_rpp_fcn(ri)
    #chris_thrth = chris_thrth_fcn(ri)
    chris_prp = chris_prp_fcn(ri)

    # geodesic equation with LAMBDA -> - LAMBDA, so all equations got multiplied by a minus sign. this is because we're integating backwards in reference to 
    # the forward integrated geodesic, which follows the evolution of getGeodesicEvolution. Therefore, if your radial velocity as found here is negative, it's actually positve,
    # because all of these differentials are d(thing)/d(-lambda). so if dr = u and u is negative, then this means dr/d(-lambda) = u is negative, so -u is positive, therefore 
    # dr/dlambda is positive. confusing but it works LOL 

    dkt = 2 * chris_ttr * kri * kti 
    dkr = chris_rpp * kphi**2 + chris_rthth * kthi**2 + chris_rrr * kri**2 + chris_rtt * kti**2
    dkth = 0
    dkp = 2*chris_prp * kri * kphi + 2*chris_pthp * kphi * kthi
    dt = -kti
    dr = -kri
    dth = -kthi
    dp = -kphi

    # the intensity hardly grows outside of r > 25, so no point in integrating past this point

    if ri < 25:

        # these if statements assign which accretion model we're considering, and fetch the appropriate equations (or data files) from equation_list

        if accword == "stat":

            u_r = 0
            u_t_stat_fcn = expressions[17]
            u_t = -1 * u_t_stat_fcn(ri)

        if accword == "radin":

            u_r_fcn = expressions[15]
            u_r = u_r_fcn(ri)
            u_t = -1.
        
        if accword == "bondi":

            gdrr = expressions[1]
            gutt = expressions[4]
            PATH = "./data/vels/"
            # this needs to be overhauled somewhat, hopefully i'll one day be able to run getBondiVels right in this code to fully streamline things. but radial infall is currently the goal. 
            # so this is a #FutureMeProblem

            if gamma < 1.5:
                g = "139"
            else:
                g = "53"

            if rcrit < 20:
                rcrit = "10.0"
            else:
                rcrit = "30.0"

            if zeta == 0: 
                filename = "vData_GR_g" + g + "_rcrit" + rcrit + ".csv"
                bondi_data = np.genfromtxt(PATH+filename, delimiter=",") # Import radius and velocity values from a pre-run code that found the velocities

            else: 
                filename = "vData_SGB_z" + str(zeta) + "_g" + g + "_rcrit" + rcrit + ".csv"
                bondi_data = np.genfromtxt(PATH+filename, delimiter=",") # Import radius and velocity values from a pre-run code that found the velocities

            table_r = bondi_data[:,0] # R values from table
            table_u = bondi_data[:,1] # ur values from table

            arb_r = var[5] # current r value

            #u_index = max(np.where(table_r < arb_r)[0]) # gets array of indices where the condition is true, then takes the maximum as the index of the velocity in the table
            condition_index = sum(1 for i in table_r if i < arb_r) # gives number of array vals, starting from beginning, that satisfy
            u_index = condition_index - 1 # use above to find array index of desired element, just minus one since the above starts counting at 1 not 0
            
            ur = (table_u[u_index+1] - table_u[u_index]) * (table_r[u_index+1] - table_r[u_index])**(-1) * (arb_r - table_r[u_index]) + table_u[u_index]    # linear interp for velocity value 
            u_r = gdrr(ri) * ur 
            u_t = - np.sqrt(-1 * (1 + gdrr(ri) * ur**2) * gutt(ri)**(-1)) # from u^{\mu} u_{\mu} = -1 


        freq = - ( kti * u_t + kri * u_r )

        dI = freq**(-2) * 4.**(-2) * np.pi**(-2) * (var[5])**(-power) # power law emissivity, where the power is an argument of the geodesic evolution 

    else:
        dI = 0

    return [dkt, dkr, dkth, dkp, dt, dr, dth, dp, dI] 