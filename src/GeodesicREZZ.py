"""
Adam M. Bauer
adammb4@illinois.edu

This code is a class object that generates a ray with some equation of motion, specified by
getGeodesicEvolution, called from a file of the same name. The object contains two objects: the ray
as traced from the camera outwards in spacetime, and the ray that was traced back to the camera.
The latter has an intensity profile associated with it as well. Other attributes of this object is the
energy and angular momentum of the ray, which can be checked as constant for the entire trajectory.

One can also get arrays of (x,y) values of the ray's trajectory if such a plot is useful. Functions with
'INT' attached to the end are associated with the ray that has an intensity profile, i.e., the ray that is 
traced back to the camera.

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
from scipy.integrate import ode, solve_ivp
from src.getGeodesicEvolutionREZZ import getGeodesicEvolution
from src.getMiscFuncsREZZ import getAngularMomentum, getEnergy, getCheck1, getCheck2, getCheck3
from src.getGeodesicEvolutionPlusREZZ import getGeodesicEvolutionPlus
import random 

import sympy as sy 

random.seed(1000) # for pretty pictures

getCheck1.terminal = True # make check in ray integration below terminal (no photon goes inside the event horizon and makes it out)

getCheck2.terminal = True # make check to not integrate past r = 30
getCheck2.direction = 1

getCheck3.terminal = True # make check to not backwards integrate after intensity is done being integrated 
getCheck3.direction = 1

class Geodesic: 

    def __init__(self, given_vals, expressions, a1, accword, power, rcrit, gamma):

        self.given_vals = given_vals # inital values 
        self.expressions = expressions # sympy functions
        self.accword = accword # accretion model
        self.power = power # power of emissivity 
        self.a1 = a1
        self.rcrit = rcrit 
        self.gamma = gamma 

        self.Earray = []
        self.Larray = []

        self.EarrayINT = []
        self.LarrayINT = []

        self.Xarray = []
        self.Yarray = []

        self.XarrayINT = []
        self.YarrayINT = []

        self.a_range = [0,10**6]
        self.atol = 10**(-8)
        self.rtol = 10**(-7)
        self.b = given_vals[4]

        self.color = (random.random(), random.random(), random.random())

        self.ray = solve_ivp(fun=lambda t, y: getGeodesicEvolution(t, y, self.expressions), t_span = self.a_range, y0 = self.getIC(), method="RK45", dense_output = True, events=(getCheck1, getCheck2), atol=self.atol, rtol=self.rtol) # make ray's trajectory an aspect of the class
        self.intensity = solve_ivp(fun=lambda t, y: getGeodesicEvolutionPlus(t, y, self.expressions, a1, accword, power, rcrit, gamma), t_span =self.a_range, y0 =self.getIntensityICs(), method="RK45", events=getCheck3, dense_output = True, atol=self.atol, rtol=self.rtol)


    def getIC(self):

        # given vals = [t, r, theta, phi, impact parameter, k_t0]
        t0 = self.given_vals[0]
        r0 = self.given_vals[1]
        th0 = self.given_vals[2]
        p0 = self.given_vals[3]

        # import metric functions from mathematica
        # all are evaluated at r = r0, and are 
        # ALREADY EVALUATED AT THETA = PI/2, M = 1
        # in the mathematica script. need to be converted into
        # numpy floats so numpy functions can be used on them 
        # down the line.

        g_tt_fcn = self.expressions[0]
        g_rr_fcn = self.expressions[1]
        g_thth_fcn = self.expressions[2]
        g_dpp_fcn = self.expressions[3]

        g_utt_fcn = self.expressions[4]
        g_urr_fcn = self.expressions[5]
        g_uthth_fcn = self.expressions[6]
        g_upp_fcn = self.expressions[7]

        g_tt = g_tt_fcn(r0)
        g_rr = g_rr_fcn(r0)
        g_thth = g_thth_fcn(r0)
        g_pp = g_dpp_fcn(r0)

        g_utt = g_utt_fcn(r0)
        g_urr = g_urr_fcn(r0)
        g_uthth = g_uthth_fcn(r0)
        g_upp = g_upp_fcn(r0)

        # raise index on k_t0

        kt0 = g_utt * self.given_vals[5]
        #kt0 = g_utt * k_t0

        # calculate k^phi0

        kp0 = g_upp * self.given_vals[4] * self.given_vals[5]
        #kp0 = g_upp * b0 * k_t0

        # normalizing four velocity gives k^r

        kr0 = - np.sqrt( -1 * (g_tt * (kt0)**2 + g_pp * (kp0)**2) * (g_rr)**(-1) )

        return [kt0, kr0, 0, kp0, t0, r0, th0, p0]

    def getIntensityICs(self):
        k = len(self.ray.t) 
        
        return [self.ray.y[0,k-1], self.ray.y[1,k-1], self.ray.y[2,k-1], self.ray.y[3,k-1], self.ray.y[4,k-1], self.ray.y[5,k-1], self.ray.y[6,k-1], self.ray.y[7,k-1], 0]

    def getEnergyArray(self):
        tmpE = 0

        g_tt_fcn = self.expressions[0]

        for i in range(0, len(self.ray.t)):
            kt0 = self.ray.y[0,i]
            r0 = self.ray.y[5, i]
            tmp_g_tt = g_tt_fcn(r0)
            tmpE = getEnergy(kt0, tmp_g_tt) # calculate energy (see MiscFuncs)
            self.Earray.append(tmpE)

    def getAngularMomentumArray(self):
        tmpL = 0

        g_dpp_fcn = self.expressions[3]

        for i in range(0, len(self.ray.t)):
            kph0 = self.ray.y[3,i]
            r0 = self.ray.y[5, i]
            tmp_g_phph = g_dpp_fcn(r0)
            tmpL = getAngularMomentum(kph0, tmp_g_phph) # calculate angular momentum (see MiscFuncs)
            self.Larray.append(tmpL)

    def getXArray(self):
        tmpx = 0

        for i in range(0, len(self.ray.t)):
            r0 = self.ray.y[5, i]
            phi0 = self.ray.y[7, i]
            tmpx = r0 * np.cos(phi0) # calculate x values for geodesic plot
            self.Xarray.append(tmpx)
        
    def getYArray(self):
        tmpy = 0

        for i in range(0, len(self.ray.t)):
            r0 = self.ray.y[5, i]
            phi0 = self.ray.y[7, i]
            tmpy = r0 * np.sin(phi0) # calculate x values for geodesic plot
            self.Yarray.append(tmpy)

    def getEnergyArrayINT(self):
        tmpE = 0

        g_tt_fcn = self.expressions[0]

        for i in range(0, len(self.intensity.t)):
            kt0 = self.intensity.y[0,i]
            r0 = self.intensity.y[5, i]
            tmp_g_tt = g_tt_fcn(r0)
            tmpE = getEnergy(kt0, tmp_g_tt) # calculate energy (see MiscFuncs) # calculate energy (see MiscFuncs)
            self.EarrayINT.append(tmpE)

    def getAngularMomentumArrayINT(self):
        tmpL = 0

        g_dpp_fcn = self.expressions[3]

        for i in range(0, len(self.intensity.t)):
            kph0 = self.intensity.y[3,i]
            r0 = self.intensity.y[5, i]
            tmp_g_phph = g_dpp_fcn(r0)
            tmpL = getAngularMomentum(kph0, tmp_g_phph) # calculate angular momentum (see MiscFuncs)
            self.LarrayINT.append(tmpL)

    def getXArrayINT(self):
        tmpx = 0

        for i in range(0, len(self.intensity.t)):
            r0 = self.intensity.y[5, i]
            phi0 = self.intensity.y[7, i]
            tmpx = r0 * np.cos(phi0) # calculate x values for geodesic plot
            self.XarrayINT.append(tmpx)
        
    def getYArrayINT(self):
        tmpy = 0

        for i in range(0, len(self.intensity.t)):
            r0 = self.intensity.y[5, i]
            phi0 = self.intensity.y[7, i]
            tmpy = r0 * np.sin(phi0) # calculate x values for geodesic plot
            self.YarrayINT.append(tmpy)