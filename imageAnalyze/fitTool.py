## This file defines some image fitting tool functions
import numpy as np
from scipy.optimize import curve_fit
import operator
from polylog import *

def initialGauss(data):
	size = np.shape(data)

	xSlice = np.sum(data,0)    
	ySlice = np.sum(data,1)
	x0 = np.argmax(xSlice)
	y0 = np.argmax(ySlice)
	offset = np.nanmin(data)
	peak = np.nanmax(data)
	amplitude = peak - offset

	a = 0
	xOff = np.nanmin(xSlice)
	maxX = np.nanmax(xSlice)-xOff
	for i in range(len(xSlice)):
		if xSlice[i] - xOff > 0.5 * maxX:
			a += 1
	b = 0
	yOff = np.nanmin(ySlice)
	maxY = np.nanmax(ySlice)-yOff
	for i in range(len(ySlice)):
		if ySlice[i] - yOff > 0.5 * maxY:
			b += 1  

	return [x0, y0, a, b, amplitude, offset]

def gaussionDistribution(coordinates, x0, y0, a, b, amplitude, offset):
	"""gaussionParams = ((x0, y0, a, b, amplitude, offset)) """
	# print coordinates
	# print gaussionParams

	# x0, y0, a, b, theta0, amplitude, offset = gaussionParams

	dist = offset + amplitude * np.exp(- (coordinates[0] - x0) **2/a**2 - (coordinates[1] - y0)**2/b**2)
	return dist.ravel()


def fermionDistribution(coordinates, x0, y0, a, b, amplitude, offset, q):
	tmp = q - ((coordinates[0]-x0)**2/a**2 + (coordinates[1]-y0)**2/b**2)* f(np.exp(q))
	numerator = fermi_poly2(tmp.ravel())
	denumerator = fermi_poly2(q)
	dist = offset + amplitude * numerator/denumerator
	return dist


def bosonDistribution(x, y, bosonParams):
	"""BosonParams = ?"""

def fitData(data, distribution):
	size = np.shape(data)
	guess = initialGauss(data)
	if distribution == fermionDistribution:
		guess.append(1)
	elif distribution == bosonDistribution:
		guess.append(0.5, 0.5)

	coordinates = np.meshgrid(range(size[0]), range(size[1]))
	
	params, Cover = curve_fit(distribution, coordinates, data.ravel(), p0=guess)

	return params

def f(x):
    # if x == 0:
    #     return 0
    # if x < -1:
    #     return 0
    return (1+x)/x * np.log(1+x)


