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

def getCriticalB(zeta):
    return np.sqrt(27) - 4397*zeta*(2430 * np.sqrt(3) )**(-1) + zeta**2 * (113374080 * np.sqrt(3))**(-1) *336780431 + zeta**3 * 114832399336847 * (np.sqrt(3) * 37192366944000)**(-1) + zeta**4 * 125183193573305833463 * (np.sqrt(3) * 52057412164177920000)**(-1)

def getBRange(N, zeta):

    lowbound = getCriticalB(zeta) - 0.5
    upbound = getCriticalB(zeta) + 0.5

    finerange = np.linspace(lowbound, upbound, N/2)
    beginning = np.linspace(0,lowbound, N/4)
    end = np.linspace(upbound, 15, N/4)

    brange = np.hstack((beginning, finerange, end))

    return brange

def getExpressions(keyword, zeta):

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

    if keyword == "SGB":

        x = sy.symbols('x') # x is zeta from my mathematica code, sympy can't parse greek letters unfortunately 
        r = sy.symbols('r')

        str_equations = []
        expressions = []

        with open('./mathematica/SGBEquations.txt', 'r') as filehandle:
            filecontents = filehandle.readlines()

            for line in filecontents:
                tmp_eqn = line[:-1]
                str_equations.append(tmp_eqn)

        for i in range(0, len(str_equations)):
            tmp_eqn = mathematica.parse(str_equations[i])
            tmp_sympy = mathematica.sympify(tmp_eqn)
            tmp_sympy_zetaeval = tmp_sympy.subs(x,zeta) # sub in x = zeta for the assigned value of zeta. this is slower than just lambdifying for both r and zeta, but it doesn't make a difference since we only sub once
            tmp_fcn = sy.lambdify(r, tmp_sympy_zetaeval)
            expressions.append(tmp_fcn)
        
        return expressions

# this function has the same energy as the previous function, just with significantly more variables
# since the accretion problem is more heavily parameterized than gravity theories
    
def getBondiExpressions(keyword, zeta, mdot, gamma, rcrit):
    if keyword == "SGBBondi":

        x = sy.symbols('x')
        m = sy.symbols('m')
        r = sy.symbols('r')
        u = sy.symbols('u')
        g = sy.symbols('g')
        b = sy.symbols('b')

        str_equations = []
        expressions = []

        with open('./mathematica/SGBBondi.txt', 'r') as filehandle:
            filecontents = filehandle.readlines()

            for line in filecontents:
                tmp_eqn = line[:-1]
                str_equations.append(tmp_eqn)

        for i in range(0, len(str_equations)):
            tmp_eqn = mathematica.parse(str_equations[i])
            tmp_sympy = mathematica.sympify(tmp_eqn)
            tmp_sympy_subs = tmp_sympy.subs([(x,zeta), (g, gamma), (b,rcrit) , (m, mdot)])
            tmp_fcn = sy.lambdify([r,u], tmp_sympy_subs)
            expressions.append(tmp_fcn)

        return expressions

    if keyword == "REZZBondi":

        v = sy.symbols('v')
        m = sy.symbols('m')
        r = sy.symbols('r')
        u = sy.symbols('u')
        g = sy.symbols('g')
        b = sy.symbols('b')

        str_equations = []
        expressions = []

        with open('./mathematica/REZZBondi.txt', 'r') as filehandle:
            filecontents = filehandle.readlines()

            for line in filecontents:
                tmp_eqn = line[:-1]
                str_equations.append(tmp_eqn)

        for i in range(0, len(str_equations)):
            tmp_eqn = mathematica.parse(str_equations[i])
            tmp_sympy = mathematica.sympify(tmp_eqn)
            tmp_sympy_subs = tmp_sympy.subs([(v,zeta), (g, gamma), (b,rcrit) , (m, mdot)])
            tmp_fcn = sy.lambdify([r,u], tmp_sympy_subs)
            expressions.append(tmp_fcn)

        return expressions

    if keyword == "GRBondi":

        str_equations = []
        expressions = []

        r = sy.symbols('r')
        m = sy.symbols('m')
        u = sy.symbols('u')
        g = sy.symbols('g')
        b = sy.symbols('b')

        with open('./mathematica/GRBondi.txt', 'r') as filehandle:
            filecontents = filehandle.readlines()

            for line in filecontents:
                tmp_eqn = line[:-1]
                str_equations.append(tmp_eqn)

        for i in range(0, len(str_equations)):
            tmp_eqn = mathematica.parse(str_equations[i])
            tmp_sympy = mathematica.sympify(tmp_eqn)
            tmp_sympy_subs = tmp_sympy.subs([(g, gamma), (b,rcrit) , (m, mdot)])
            tmp_fcn = sy.lambdify([r,u], tmp_sympy_subs)
            expressions.append(tmp_fcn)

        return expressions

# generates the initial conditions for the geodesics using adaptive mesh spacing 
def getIVsAndNames(N, zeta):
    brange = getBRange(N, zeta)
    finalIVarray = []

    tmpstr1 = ""
    finalNamearray = []
    
    for i in range(0, N):
        #print(tmpb)
        tmparray = [0., 1000., np.pi/2, 0., brange[i], -10.] # [t0, r0, theta0, phi0, b, k_t]
        finalIVarray.append(tmparray)

        tmpstr1 = "geo" + str(i+1)
        finalNamearray.append(tmpstr1)
    
    return finalNamearray, finalIVarray

def getIVsAndNames_emisfig(brange, zeta):
    finalIVarray = []

    tmpstr1 = ""
    finalNamearray = []

    N = len(brange)
    
    for i in range(0, N):
        #print(tmpb)
        tmparray = [0., 1000., np.pi/2, 0., brange[i], -10.] # [t0, r0, theta0, phi0, b, k_t]
        finalIVarray.append(tmparray)

        tmpstr1 = "geo" + str(i+1)
        finalNamearray.append(tmpstr1)
    
    return finalNamearray, finalIVarray

def getCheck1(t,y):
    return y[5] - (2.01 - 49 * y[8] * 40.**(-1)) # stops integration at the event horizon (y[5] = r, y[8] = zeta; it has to be written in terms of t and y for the ODE integrator to recognize it)

def getCheck2(t,y):
    return y[5] - 25.05# stop forward integration at r = 30 if it doesn't stop at the horizon (for all the geos that sling around the BH and don't fall in)

def getCheck3(t,y):
    return y[5] - 25. # stop backward intensity integration at r = 26, since there is no point in integrating after the intensity integration is over.

def getIntegration(integrand, xvals):
    tmp = 0

    for i in range(1, len(xvals)):
        tmp += 0.5 * (xvals[i] - xvals[i-1]) * (integrand[i] + integrand[i-1])
    
    return tmp