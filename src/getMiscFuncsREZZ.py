"""
Adam M. Bauer
adammb4@illinois.ed

The file contains numerous functions that are used by MasterIntensityCalculator.py, as well
as other programs, such as Geodesic.py. Importantly, this file contains the symbolic algebra 
parsing that is used throughout the simulation to change between GR and sGB seamlessly. This uses
sympy. The 'checks' are all checks that are build into the 
numerical integration of the geodesic equation in Geodesic.py.

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
from scipy import optimize

from sympy.parsing import mathematica  
import sympy as sy 

def getEnergy(kt, g_tt):
    return - g_tt * kt # calc energy (from killing vector equation)

def getAngularMomentum(kp, g_phph):
    return g_phph * kp # calc angular momentum (from killing vector equation)

def getCriticalB(a1):
    return np.sqrt(27) - 4 * a1 *(np.sqrt(27))**(-1)

def getBRange(N, zeta):

    lowbound = getCriticalB(zeta) - 0.5
    upbound = getCriticalB(zeta) + 0.5

    finerange = np.linspace(lowbound, upbound, N/2)
    beginning = np.linspace(0,lowbound, N/4)
    end = np.linspace(upbound, 15, N/4)

    brange = np.hstack((beginning, finerange, end))

    return brange

def getExpressions(keyword, a1):

    if keyword == "GR":

        str_equations = []
        expressions = []

        r = sy.symbols('r') # make r a sympy symbol, i.e., a variable that sympy functions can "be a function of" later

        with open('./mathematica/GREquations.txt', 'r') as filehandle:
            filecontents = filehandle.readlines()

            for line in filecontents:
                tmp_eqn = line[:-1]
                str_equations.append(tmp_eqn)

        for i in range(0, len(str_equations)):
            tmp_eqn = mathematica.parse(str_equations[i]) # uses sympy parsing to take the string and form a function sympy can read
            tmp_sympy = mathematica.sympify(tmp_eqn) # makes interpretable string tmp_eqn into a sympy function, with variables defined as symbols above
            tmp_fcn = sy.lambdify(r, tmp_sympy) # makes the sympy function into a lambda function which is ***significantly faster*** than just straight subbing using sy.subs
            expressions.append(tmp_fcn) # append to expressions list! 
        
        return expressions

    if keyword == "REZZ":

        v = sy.symbols('v') # v = a1 from my mathematica code, sympy can't parse greek letters unfortunately 
        r = sy.symbols('r')

        str_equations = []
        expressions = []

        with open('./mathematica/REZZEquations.txt', 'r') as filehandle:
            filecontents = filehandle.readlines()

            for line in filecontents:
                tmp_eqn = line[:-1]
                str_equations.append(tmp_eqn)

        for i in range(0, len(str_equations)):
            tmp_eqn = mathematica.parse(str_equations[i])
            tmp_sympy = mathematica.sympify(tmp_eqn)
            tmp_sympy_zetaeval = tmp_sympy.subs(v, a1) # sub in x = zeta for the assigned value of zeta. this is slower than just lambdifying for both r and zeta, but it doesn't make a difference since we only sub once
            tmp_fcn = sy.lambdify(r, tmp_sympy_zetaeval)
            expressions.append(tmp_fcn)
        
        return expressions

# generates the initial conditions for the geodesics using adaptive mesh spacing 
def getIVsAndNames(N, a1):
    finalIVarray = []
    bs = getBRange(N, a1)

    tmpstr1 = ""
    finalNamearray = []
    
    for i in range(0, N):
        tmparray = [0., 1000., np.pi/2, 0., bs[i], -10.] # [t0, r0, theta0, phi0, b, k_t]
        finalIVarray.append(tmparray)

        tmpstr1 = "geo" + str(i+1)
        finalNamearray.append(tmpstr1)
    
    return finalNamearray, finalIVarray

def getCheck1(t,y):
    return y[5] - 2.01 # stops integration at the event horizon (y[5] = r, y[8] = eps; it has to be written in terms of t and y for the ODE integrator to recognize it)

def getCheck2(t,y):
    return y[5] - 25.05 # stop forward integration at r = 30 if it doesn't stop at the horizon (for all the geos that sling around the BH and don't fall in)

def getCheck3(t,y):
    return y[5] - 25. # stop backward intensity integration at r = 26, since there is no point in integrating after the intensity integration is over.

def getIntegration(integrand, xvals):
    tmp = 0

    for i in range(1, len(xvals)):
        tmp += 0.5 * (xvals[i] - xvals[i-1]) * (integrand[i] + integrand[i-1])
    
    return tmp