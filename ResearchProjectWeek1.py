#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 22:31:11 2019

@author: RachelCampo
"""

#I am going to attempt to try and code a program to analyze my data
#without my data..... Here we go

import math
import fileinput
import numpy as np
import scipy as sp
import scipy.interpolate
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.misc import derivative
from scipy.integrate import odeint
from matplotlib import rcParams

#this code is going to be the part where it runs through my data and figures
#out what the wavelength is after the redshift, uses the dopler shift to find
#what the wavelength was before the redshift, and then determines what 
#part of the spectrum is was from.

 #wavelength = something from the data

for x in datalist:
    