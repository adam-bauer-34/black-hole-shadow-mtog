# Generating synthetic black hole images
Written by: Adam Bauer

Affiliation: Department of Physics, University of Illinois at Urbana-Champaign, 1110 West Green St, Urbana, IL 61801, USA

To contact: adammb4 [at] illinois [dot] edu

## Contents:
1. [Code description](#codedesc)
2. [Publication](#pubs)
3. [Model assumptions](#ass)
   1. [Astrophysical model assumptions](#astroass)
   2. [Gravitational model assumptions](#gravass)
4. [Directory overview](#directover)
   1. [src](#src)
   2. [mathematica](#math)
   3. [data](#data)
   4. [analysis](#anal)
6. [Copyright statement](#dontgetgot)

## Code description <a name="codedesc"></a>
This code can be used to generate the intensity profile of plasma being accreted onto a black hole, which can then be transformed into an Event Horizon Telescope-esque image. Other science can be done on the intensity profile as well. 

## Publication <a name="pubs"></a>
We have submitted a publication which uses this code and its extensions; I only put here what I think will be of general use to the science community. Hence, *not every piece of code that was used in the paper is here.* If you’re interested in a specific calculation, don’t hesitate to reach out at the email above. I’m happy to help and share more of my science! A doi (or arxiv link) will be added when it is available. 

Paper title: *Spherical accretion in alternative theories of gravity*

Abstract: The groundbreaking image of the black hole at the center of the M87 galaxy has raised questions at the intersection of observational astronomy and black hole physics. How well the radius of the black hole shadow can be measured, and can this measurement be used to distinguish general relativity from other theories of gravity? We explore these questions using a simple spherical flow model in general relativity, scalar Gauss—Bonnet gravity, and a parameterized metric. We assume an optically thin plasma with power-law emissivity in radius. Along the way we provide a derivation of a generalized Bondi flow, as well as a piecewise-analytic model for the brightness profile of a cold inflow. We use the second moment of the image as a proxy for EHT observables, and compute the ratio of the second moment to the radius of the black hole shadow.   We show that corrections to this ratio from modifications to general relativity are subdominant compared to corrections to the critical impact parameter, and argue that this is generally true.  We find that astrophysical model parameters are the dominant source of uncertainty in this calculation, which highlights the importance of understanding the astrophysical model. Given a sufficiently accurate astrophysical model, however, it is possible using measurements of the black hole shadow to distinguish between general relativity and other theories of gravity.

Arxiv link: TBD

## Model assumptions <a name="ass"></a>
For a full overview of the model and its implementation, see the paper mentioned above, whenever it is available. :P

### Astrophysical model assumptions <a name="astroass"></a>
We make the following assumptions about the astrophysical flow:

#### Accretion models
This code supports the use of three accretion models, all of which are spherically symmetric and time independent: stationary emitting gas, radially infalling gas (i.e., Bondi accretion at zero temperature) and Bondi accretion. 

#### Emissivity model

This code only currently supports a frame-dependent emissivity that is described by a power law, i.e.,
j<sub>&nu;</sub> = j<sub>0</sub> (r<sub>0</sub>/r)<sup>&alpha;</sup>.
In our code,
r<sub>0</sub> = 1 and j<sub>0</sub> = (16 &pi;<sup>2</sup>)<sup>-1</sup>. This can be easily changed, however, in ```getGeodesicEvolutionPlus```, see below.

### Gravity model assumptions <a name="gravass"></a>
We assume there exists a black hole at the origin. This black hole has spherical symmetry and its metric is time independent. The code is, in principle, amendable with any theory of gravity whose black hole possesses this quality of spherical symmetry and time independence, given the correct input files are made. The below gives an overview of how the code is structured.

## Directory overview <a name="directover"></a>
In the top directory are the Master files. Those can be run to generate the intensity profile data. They rely on files in the ```mathematica``` and ```src``` folders.

### src <a name="src"></a>
Src has the dependencies fo the master file in it. Note SGB and REZZ have different dependencies. 
- Geodesic: A class that describes one ray of light. Can compute the energy and angular momentum of the ray as a sanity check.
- getGeodesicEvolution: Contains the equation of motion for a geodesic being forward integrated.
- getGeodesicEvolutionPlus: Contains equation of motion for backwards integration *plus* intensity calculation along the ray. 
- getMiscFuncs: Miscellaneous set of functions needed for integration.
- getBondiVels: Uses root-finding algorithm to calculate the bondi flow velocities in a generic theory of gravity.

### mathematica <a name="math"></a>
Contains numerous mathematica scripts for calculating the relevant equations for ray-tracing and intensity calculations. 

Files named get[theory]Equations.nb output a .txt file with relevant equations (metric components, Christoffel symbols, flow velocity in radial infall/stationary case) to be used by Master files.

Files named get[theory]Bondi.nb calculate the bondi flow equations that are used in getBondiVels to calculate the flow velocity.

These scripts should serve as a template for any other desirable theory of gravity that one would like to investigate.

### data <a name="data"></a>
Contains data for all accretion scenarios used in *Spherical accretion…* as well as Bondi velocity data. Any runs of Master files automatically save their output here in the relevant subdirectory.

### analysis <a name="anal"></a>
A few codes for model output analysis. 

#### computation
Contains code to calculate the ratio of the characteristic radius of emission to critical impact parameter, denoted as &vartheta;, in *Spherical accretion…*. Versions for both bondi and radial infall cases.

#### shadow
Contains code to make a grid of synthetic Event Horizon Telescope images.

## Copyright statement <a name="dontgetgot"></a>
Copyright (C) 2021, Adam M. Bauer

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
