## This file defines some image fitting tool functions
import numpy as np
from scipy.optimize import curve_fit
import operator
from polylog import *
import signal

def signal_handler(signum, frame):
	print "Timed Out!"
	raise Exception("Timed out!")



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

	dist = offset + amplitude * np.exp(- (coordinates[0] - x0) **2/a**2 - (coordinates[1] - y0)**2/b**2)
	return dist.ravel()


def fermionDistribution(coordinates, x0, y0, a, b, amplitude, offset, q):
	tmp = q - ((coordinates[0]-x0)**2/a**2 + (coordinates[1]-y0)**2/b**2)* f(np.exp(q))
	numerator = fermi_poly2(tmp.ravel())
	denumerator = fermi_poly2(q)
	dist = offset + amplitude * numerator/denumerator
	return dist


def bosonDistribution(coordinates, x0, y0, a, b, amplitudeC, offset, amplitudeT, Ca, Cb):
	"""BosonParams = ?"""
	thermalPart = amplitudeT * np.exp(- (coordinates[0] - x0) **2/a**2 - (coordinates[1] - y0)**2/b**2)
	condensatePart = amplitudeC *  np.maximum( (1 - Ca * (coordinates[0] - x0)**2 - Cb * (coordinates[1] - y0)**2), 0)
	dist = thermalPart + condensatePart + offset
	return dist.ravel()

def fitData(data, distribution, mode):
	if mode == "Single":
		signal.signal(signal.SIGALRM, signal_handler)
		signal.alarm(15)   # Ten seconds
		try:
			size = np.shape(data)
			guess = initialGauss(data)
			
			if distribution == fermionDistribution:
				guess.append(1)
			elif distribution == bosonDistribution:
				guess.append(1)
				guess.append(0.1)
				guess.append(0.1)
			
			coordinates = np.meshgrid(range(size[0]), range(size[1]))
			
			params, Cover = curve_fit(distribution, coordinates, data.ravel(), p0=guess)
			signal.alarm(0) 
			
			return params
		except Exception, msg:
			print "Timed out! This fitting looks much slower than usual. Please check your AOI!" 
	
	
	
	else:
		size = np.shape(data)
		guess = initialGauss(data)
		
		if distribution == fermionDistribution:
			guess.append(1)
		elif distribution == bosonDistribution:
			guess.append(1)
			guess.append(0.1)
			guess.append(0.1)
		
		coordinates = np.meshgrid(range(size[0]), range(size[1]))
		
		params, Cover = curve_fit(distribution, coordinates, data.ravel(), p0=guess)
		
		
		return params
	

def f(x):
    # if x == 0:
    #     return 0
    # if x < -1:
    #     return 0
    return (1+x)/x * np.log(1+x)


def radioDistribution(data, center, sigma):
	# print center
	# print sigma
	size = np.shape(data)
	# print size
	x1 = min(center[0], size[0]-center[0])/float(sigma[0])
	y1 = min(center[1], size[1]-center[1])/float(sigma[1])
	r0 = min(x1, y1)

	lr = int(0.95*r0)
	od_list = []

	for r in np.arange(0, lr, 0.01):
		od = 0
		for theta in range(0, 360, 5):
			x = center[0] + int(r*np.cos(theta) * sigma[0])
			y = center[1] + int(r*np.sin(theta) * sigma[1])
			od += data[y, x]
		od=od/360
		# print od
		od_list.append(od)

	return od_list



