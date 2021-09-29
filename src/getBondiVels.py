"""
Adam M. Bauer

This code returns a .csv file that has radius values and the radial velocity at said values 
for the Bondi accretion case. 

To run: python getVelocities.py

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
from getMiscFuncs import getBondiExpressions
from scipy import optimize
import matplotlib.pyplot as plt
import sympy as sy

def getU(u, args):

    r = args[0]

    exp1_fcn = args[1][0]
    exp2_fcn = args[1][1]
    exp1_fcn_eval = exp1_fcn(0,0)
    exp2_fcn_eval = exp2_fcn(r,u)

    print(exp1_fcn_eval, exp2_fcn_eval)

    return np.float64(exp1_fcn_eval) - exp2_fcn_eval

mdot = 1

# zetas = [0, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175] # sgb zetas 
zetas = [-0.20, -0.15, -0.10, -0.05, 0.05, 0.10, 0.15, 0.20]
gammas = [float(13 * 9.**(-1)), float(5 * 3.**(-1))]
rcrits = [10., 30.]

for i in range(0, len(zetas)):
    for j in range(0, len(gammas)):
        for k in range(0, len(rcrits)):

            print("The is zeta = " + str(zetas[i]) + ", gamma = " + str(gammas[j]) + ", and rcrit = " + str(rcrits[k]))

            equation_list = getBondiExpressions("REZZBondi", zetas[i], mdot, gammas[j], rcrits[k])

            rh = 2.

            r_vals = np.linspace(rh-0.01, 25.1, 100)
            N = len(r_vals)

            vels = np.zeros(N)

            # FOR BONDI
            for h in range(0, N):
                #print(r_vals[h])

                if r_vals[h] > rcrits[k]:
                    guess = 10**(-6)
                
                else: 
                    guess = 0.4

                args = [r_vals[h], equation_list]
                solution = optimize.root(getU, x0=guess, args=args, method='lm', tol=10**(-8))
                vels[h] = solution.x

            written_data = np.vstack((r_vals, vels)).T 

            if gammas[j] == float(13 * 9.**(-1)):
                strgamma = 'g139'
            else: 
                strgamma = 'g53'
            
            PATH = "./data/vels/"

            if zetas[i] == 0:
                theory = "GR"
                filename = PATH + "vData_" + theory + "_" + strgamma + "_rcrit" + str(rcrits[k]) + ".csv"
            else: 
                theory = "REZZ"
                if zetas[i] < 0:
                    strzeta = "neg" + str(abs(zetas[i]))
                else: 
                    strzeta = str(zetas[i])
                
                filename = PATH + "vData_" + theory + "_z" + strzeta + "_" + strgamma + "_rcrit" + str(rcrits[k]) + ".csv"

            np.savetxt(filename, written_data, delimiter=",")

            plt.plot(r_vals, vels)
            plt.hlines(1,rh,26)
            plt.vlines(rh,0,1.1)
            plt.show()

#r1 = np.linspace(3, 20, 10)
#r2 = np.flip(np.linspace(3, 20, 10), axis=0)
#rtot = np.hstack((r2,r1))

#velinterp = np.interp(rtot, r_vals, vels)

#print(velinterp)

#plt.show()