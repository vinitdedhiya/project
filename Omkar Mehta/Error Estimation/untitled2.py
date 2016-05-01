# -*- coding: utf-8 -*-
"""
Created on Tue Feb 09 23:28:17 2016

@author: Omkar Mehta
"""

import numpy as np
from scipy.optimize import curve_fit
def func(x, a, b, c):
	return a * np.exp(-b * x) + c
xdata = np.linspace(0, 4, 50)
y = func(xdata, 2.5, 1.3, 0.5)
ydata = y + 0.2 * np.random.normal(size=len(xdata))
error=ydata-y
popt, pcov = curve_fit(func, xdata, ydata)
print popt
print error