from constant import *
import numpy as np
from scipy import stats
import os, sys, struct
from math import sqrt, log, atan
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpmath import mp
from polylog import *

# list square, for fitting
def square(list):
    return map(lambda x: x ** 2, list)

# read data, changed to absorption image
def readData(path):
    #Some option for setting path

    #__location__ = os.path.realpath(
    #os.path.join(os.getcwd(), os.path.dirname(__file__)))
    #
        #temperary solution 
    #location = "C:/ExperimentImages"
#print __location__
    b = bytearray(2*3*1024*1024+100)
    f = open(path, "rb")
    try:
	numBytesRead = f.readinto(b) 
#    print "Read ", numBytesRead, " bytes."
    finally:
	f.close()


#print "Filetype:" + chr(b[0])+chr(b[1])+chr(b[2])
    numBytesPerValue = struct.unpack('h',chr(b[3])+chr(b[4]))[0]
    rowTotal = struct.unpack('h',chr(b[5])+chr(b[6]))[0]
    colTotal = struct.unpack('h',chr(b[7])+chr(b[8]))[0]
    layerTotal = struct.unpack('h',chr(b[9])+chr(b[10]))[0]

    imageData = np.zeros((layerTotal,rowTotal,colTotal))
    byteIndex = 11
    for layer in range (layerTotal):
	for row in range(rowTotal):
            for col in range (colTotal):
                imageData[layer][row][col]= struct.unpack('h',chr(b[byteIndex])+chr(b[byteIndex+1]))[0]  
		byteIndex+=2      
###  Construct the transmittance map, with an np.maximum statement to avoid dividing by zero.
    absorbImg=(imageData[0]-imageData[2])/(np.maximum(imageData[1]-imageData[2],1))

###  Replace extremely low transmission pixels with a minimum meaningful transmission. 
    minT = 0.01
    temp = np.empty((rowTotal,colTotal))	
    temp.fill(minT)


    absorbImg = np.maximum(absorbImg,temp)
    return absorbImg


### atom number #######
def atomNumber(Img, offset):
    Img[:] = [x - offset for x in Img]
    return np.sum(np.sum(Img))  * (pixelToDistance**2)/crossSection

### fit image ###
def twoDGaussianFit(data):
    """Fits a two-dimensional Gaussian, A*exp(-0.5*(((x-x0)/sigmaX)^2+((y-y0)/sigmaY)**2))+offset, to given data 
    (a numpy array) and returns the fit parameters [[x0,y0],[sigmaX,sigmaY],offset,A]. """

### Get guess values by doing 1d fits to the data integrated in the X or Y direction.
    def oneDGaussian(x, centerX,sDevX,Amp,yOffset):
        return Amp*np.exp(-0.5*((x-centerX)/sDevX)**2)+yOffset
    ### Mask out only the area of interest then integrate the data long each axis
    
    xSlice = np.sum(data,0)    
    ySlice = np.sum(data,1)
        
    ### Initial guesses for 1d fits
    xOff = np.nanmin(xSlice)
    AmpX = np.nanmax(xSlice)-xOff
    x0 = np.argmax(xSlice)
    
    yOff = np.nanmin(ySlice)
    AmpY = np.nanmax(ySlice)-yOff
    y0 = np.argmax(ySlice)
    sigmaX =40
    sigmaY =40

    ### 1d fits
    xVals, yCovar = curve_fit(oneDGaussian,range(len(xSlice)),xSlice,p0=(x0,sigmaX,AmpX,xOff))
    x0 = xVals[0]
    sigmaX = xVals[1]

    
    yVals, yCovar = curve_fit(oneDGaussian,range(len(ySlice)),ySlice,p0=(y0,sigmaY,AmpY,yOff))
    y0 = yVals[0]
    sigmaY = yVals[1]
    
    ### 2d fit
#    def twoDGaussian(xVec, centerX,sDevX,Amp,yOffset):
#        return Amp*np.exp(-0.5*((x-centerX)/sDevX)^2)+yOffset
    offset = 0.5*(xVals[3]/(np.shape(data)[1]) + yVals[3]/(np.shape(data)[0]))
    A = 0.5 * ( xVals[2]/(sqrt(2.0*np.pi)*sigmaY) + yVals[2]/(sqrt(2.0*np.pi)*sigmaX) ) 
    
    return [[x0,y0],[sigmaX,sigmaY],A,offset]
    
# condensate fit
def twoDParbolicFit(data):
    """parabolic fit the function z = Amp * max(0, 1 - a*(x-x0)^2 - b*(y-y0)^2) +offset). 
    Since we only need center and width, this function returns center, width, amplitude and offset.
    """
    def oneDParabolic(x, centerX, Amp, a, offset):
        #print  - Amp *  ((x-centerX)**2) + C
        #print offset
        #return - Amp *  ((x-centerX)**2) + C * offset
        out = Amp *  np.maximum( (1 - a * (x - centerX)**2), 0)  ** 2+ offset

        return out
    #def oneDGaussian(x, centerX,sDevX,Amp,yOffset):
    #    return Amp*np.exp(-0.5*((x-centerX)/sDevX)**2)+yOffset
    ### Mask out only the area of interest then integrate the data long each axis
    



    xSlice = np.sum(data,0)    
    ySlice = np.sum(data,1)
        
    ### Initial guesses for 1d fits
    xOff = np.nanmin(xSlice)
    maxX = np.nanmax(xSlice)-xOff
    x0 = np.argmax(xSlice)
    
    yOff = np.nanmin(ySlice)
    maxY = np.nanmax(ySlice)-yOff
    y0 = np.argmax(ySlice)

    lengthX = 0
    for i in range(len(xSlice)):
        if xSlice[i] - xOff > 0.2 * maxX:
            lengthX += 1
 
    lengthY = 0
    for i in range(len(ySlice)):
        if ySlice[i] - yOff > 0.2 * maxY:
            lengthY += 1  
    
    aX = 4./lengthX**2
    aY = 4./lengthY**2
    AmpX = 4./3. * maxX/sqrt(aY)
    AmpY = 4./3. * maxY/sqrt(aX)

            
    ### 1d fits
    xVals, yCovar = curve_fit(oneDParabolic,range(len(xSlice)),xSlice,p0=(x0, AmpX, aX, xOff))
    x0 = xVals[0]
    widthX = sqrt(1/xVals[2])
    
    yVals, yCovar = curve_fit(oneDParabolic,range(len(ySlice)),ySlice,p0=(y0, AmpY, aY, yOff))
    y0 = yVals[0]
    widthY = sqrt(1/yVals[2])

    Amp = sqrt(xVals[1] * yVals[1] * 9./16. * sqrt(xVals[2] * yVals[2]) )
    
    offset = 0.5 * (xVals[3]/(np.shape(data)[1]) + yVals[3]/(np.shape(data)[0]))
    
    return [[x0, y0], [widthX, widthY], Amp, offset]


# partly condensate fit
def partlyCondensateFit(data):

    def oneDPartlyCondensate(x, centerX, AmpP, a, sDevX, AmpG, offset):
        return np.maximum(AmpG, 0)*np.exp(-0.5*((x-centerX)/sDevX)**2) + np.maximum(AmpP, 0) *  np.maximum( (1 - a * (x - centerX)**2), 0)  ** 2.5+ offset
       

    xSlice = np.sum(data,0)    
    ySlice = np.sum(data,1)
        
    ### Initial guesses for 1d fits
    xOff = np.nanmin(xSlice)
    AmpGX = np.nanmax(xSlice)-xOff
    x0 = np.argmax(xSlice)
    
    yOff = np.nanmin(ySlice)
    AmpGY = np.nanmax(ySlice)-yOff
    y0 = np.argmax(ySlice)
    sigmaX =40
    sigmaY =40

    lengthX = 0
    for i in range(len(xSlice)):
        if xSlice[i] - xOff > 0.5 * AmpGX:
            lengthX += 1
 
    lengthY = 0
    for i in range(len(ySlice)):
        if ySlice[i] - yOff > 0.5 * AmpGY:
            lengthY += 1  


    ### Standard deviation is a bit more work to estimate.
    check = AmpGX*np.exp(-0.5)
    for index in range(x0,len(xSlice)):
        if xSlice[index] < check:
            sigmaX = index-x0
            break
    check = AmpGY*np.exp(-0.5)
    for index in range(y0,len(ySlice)):
        if ySlice[index] < check:
            sigmaY = index-y0
            break
    
    aX = 4./lengthX**2
    aY = 4./lengthY**2
    # AmpPX = 4./3. * AmpGX/sqrt(aY)
    # AmpPY = 4./3. * AmpGY/sqrt(aX)
    AmpPX = AmpGX
    AmpPY = AmpGY

            
    ### 1d fits
    xVals, yCovar = curve_fit(oneDPartlyCondensate,range(len(xSlice)),xSlice,p0=(x0, AmpPX/2, aX, sigmaX, AmpGX/2, xOff))
    x0 = xVals[0]
    widthX = sqrt(1/xVals[2])
    sigmaX = xVals[3]
    
    yVals, yCovar = curve_fit(oneDPartlyCondensate,range(len(ySlice)),ySlice,p0=(y0, AmpPY/2, aY, sigmaY, AmpGY/2, yOff))
    y0 = yVals[0]
    widthY = sqrt(1/yVals[2])
    sigmaY = yVals[3]

    # AmpP = np.maximum(3./8. * xVals[2] * xVals[1], 0) + np.maximum(3./8. * yVals[2] * yVals[1], 0)
    AmpP = sqrt(np.maximum(xVals[1], 0) * np.maximum(yVals[1], 0) * 9./16. * sqrt(xVals[2] * yVals[2]) )
    AmpG = 0.5 * ( np.maximum(xVals[4], 0)/(sqrt(2.0*np.pi)*sigmaY) + np.maximum(yVals[4], 0)/(sqrt(2.0*np.pi)*sigmaX) ) 
    offset = 0.5 * (xVals[5]/(np.shape(data)[1]) + yVals[5]/(np.shape(data)[0]))

    return [[x0,y0], [widthX, widthY], [sigmaX,sigmaY], AmpP, AmpG, offset]

# for fermionic fit
def f(x):
    # if x == 0:
    #     return 0
    # if x < -1:
    #     return 0
    return (1+x)/x * np.log(1+x)

# fermionic fit
def fermionFit(data):

    def oneDPolylog(x, centerX, Rx, Amp , q, yOffset):
        x = np.array(x)
        
        temp = np.exp(q)
        numerator = polylog5half(-np.exp(q - (x-centerX)**2/Rx**2 * f(np.exp(q))))
        denuminator = polylog5half(-np.exp(q))
        out = numerator/denuminator * Amp + yOffset

        return out

        
        
    ### Mask out only the area of interest then integrate the data long each axis


    xSlice = np.sum(data,0)    
    ySlice = np.sum(data,1)
        
    ### Initial guesses for 1d fits
    xOff = np.nanmin(xSlice)
    AmpX = np.nanmax(xSlice)-xOff
    x0 = np.argmax(xSlice)
    
    yOff = np.nanmin(ySlice)
    AmpY = np.nanmax(ySlice)-yOff
    y0 = np.argmax(ySlice)
    # sigmaX = sXguess
    # sigmaY = sYguess
    sigmaX = 40
    sigmaY = 40
    q0 = 4

            
    ### 1d fits
    xVals, yCovar = curve_fit(oneDPolylog,range(len(xSlice)),xSlice,p0=(x0,sigmaX,AmpX, q0 ,xOff))
    yVals, yCovar = curve_fit(oneDPolylog,range(len(ySlice)),ySlice,p0=(y0,sigmaY,AmpY, q0 ,yOff))
    
    x0 = float(xVals[0])
    RX = float(xVals[1])

    y0 = float(yVals[0])
    RY = float(yVals[1])
    
    qx=float(xVals[3])
    qy=float(yVals[3])
    offset = float(0.5*(xVals[4]/(np.shape(data)[1]) + yVals[4]/(np.shape(data)[0])))
    A = np.array(data).max()
    
   
    return [[x0,y0],[RX,RY],A,[qx, qy],offset]
    
### extract static data #####

def temperatureSingleGaussianFit(ToF, gSigmaX, gSigmaY, OmegaAxial, OmegaRadial, atom):
    if atom == 'Li':
        m = mLi
    elif atom == 'Na':
        m = mNa
    else:
        return false
    
    Tempx = (m * gSigmaX**2 * OmegaAxial ** 2)/(2*kB*(1 + OmegaAxial ** 2 * ToF ** 2))
    Tempy = (m * gSigmaY**2 * OmegaRadial ** 2)/(2*kB*(1 + OmegaRadial ** 2 * ToF ** 2))

    # averageTemp = 0.5*(Tempx + Tempy)
    # deltaTemp = 0.5*abs(Tempx - Tempy)
    return [Tempx, Tempy]

# single shot -> trap frequency radial
def trapFrequencyRadial(ToF, rho, rho0):
    freq = (sqrt((rho/rho0)**2 - 1))/ToF
    return freq

# single shot -> trap frequency Axial
def trapFrequencyAxial(ToF, z, rho0, omegaRadial):
    tau = ToF * omegaRadial
    temp = tau * atan(tau) - log(sqrt(1+ tau**2))
    a = temp
    b = - z/rho0

    e = (sqrt(b**2-4*a) - b)/(2*a)

    return e*omegaRadial
    
# single shot -> chamical potential
def chemicalPotential(ToF, omegaRadial, rho0):
    m = mLi
    mu = 0.5 * m * ((omegaRadial**2/(1+(omegaRadial**2) * (ToF**2))) * rho0**2)
    return mu

#def scatteringLength(omegaRadial, omegaAxial):
#    m = mNa
#    omegaBar = omegaRadial**(2/3) * omegaAxial ** (1/3)
#    return sqrt(hbar/(omegaBar*m))
    
# used to get atom number    
def effectiveInteraction(a):
    m = mLi
    return 4*np.pi*(hbar**2)*a/m
    
# get atom number from chemical potential
def atomNumberFit(mu, omegaRadial, omegaAxial, U0):
    m= mLi
    omegaBar = (omegaRadial**2 * omegaAxial) ** (1./3.)
    #print 'omegaBar'
    #print omegaBar
    #print '(2*mu)/(m*omegaBar**2))**(3/2)'
    #print ((2*mu)/(m*omegaBar**2))**(3/2)
    N = (8 *np.pi/15) * ((2*mu)/(m*omegaBar**2))**(3./2.) * (mu/U0)
    return N
        
    

###### multiple shots functions ##########

#fit temperature from thermal gas
def fitTemperature(arr1, arr2, arr3, arr4, arr5):
    m = mLi

    timeOfFlight = arr1
    RX = arr2
    RY = arr3
    qX = arr4
    qY = arr5
    templist1 = np.zeros(len(arr1))
    templist2 = np.zeros(len(arr1))
    for i in range(len(arr1)):
        templist1[i] = RX[i]**2/f(np.exp(qX[i]))
        templist2[i] = RY[i]**2/f(np.exp(qY[i]))

    linregressX = stats.linregress(templist1, np.array(square(timeOfFlight)))
    # linregressX = stats.linregress(np.array(square(timeOfFlight)), templist1)
    slopeX = linregressX[0]


    Tempx = m  / (2 * kB * slopeX) 
 
    # print("temperature calculated from X is %f nK" %(temperatureX * 1E9))
    linregressY = stats.linregress(templist2, np.array(square(timeOfFlight)))
    # linregressY = stats.linregress(np.array(square(timeOfFlight)), templist2)
    slopeY = linregressY[0]
    Tempy = m / (2 * kB * slopeY)
    # print("temperature calculated from Y is %f nK" %(temperatureY * 1E9))

   
    return [Tempx, Tempy]


def fitInTrapRadialHalfWidth(tofList, pHalfWidthRadialList, TrapFreqRadial):
    
    return 0

def fitTrapFrequency(arr1, arr2, arr3, arr4, arr5):
    m = mLi

    timeOfFlight = arr1
    RX = arr2
    RY = arr3
    qX = arr4
    qY = arr5
    templist1 = np.zeros(len(arr1))
    templist2 = np.zeros(len(arr1))
    for i in range(len(arr1)):
        templist1[i] = RX[i]**2/f(np.exp(qX[i]))
        templist2[i] = RY[i]**2/f(np.exp(qY[i]))


    linregressX = stats.linregress(templist1, np.array(square(timeOfFlight)))
    bX = linregressX[1]
    
    omegaX = np.sqrt(1/bX)

    linregressY = stats.linregress(templist2, np.array(square(timeOfFlight)))
    bY = linregressY[1]
    
    omegaY = np.sqrt(1/bY)

    return [omegaX, omegaY]
# def fitTrapFrequency(tof, rho, z, rho_init, z_init):
    # n = len(tof)
    # rhoOverZ = []
    # for i in range(n):
    #     rhoOverZ.append(rho[i]/z[i])
    # # print rho_init
    # # print z_init
    # def rhoZRatio(t, omegaRadial, omegaAxial):
    #     e = omegaAxial/omegaRadial
    #     tau = omegaRadial * t
    #     temp1 = np.sqrt(1 + tau **2)
    #     temp2 = 1 + e**2 * (tau * np.arctan(tau) - np.log(np.sqrt(1 + tau**2)))
    #     return e * temp1/temp2 


    # trapVals, trapCov= curve_fit(rhoZRatio, tof, rhoOverZ , p0=(rho_init, z_init))
    # # print trapVals

    # templist = [None] * n
    # for i in range(n):
    #     templist[i] = 1 + trapVals[0] ** 2 * tof[i] ** 2
    # lineVals = stats.linregress(templist, rho)
    # rho0 = lineVals[0]
    # a = trapVals.tolist()
    # a.append(rho0)
    # return a

    

    
def Trap_frequency_fit_from_radius(rho_array, ToF_array):
    rhoregress = stats.linregress(square(ToF_array), square(rho_array))
    #print rhoregress
    
    ##y = rhoregress[1] + rhoregress[0] * (square(ToF_array))
    #t = [-0.02,0,0.05,0.1]
    #y = [rhoregress[1] + rhoregress[0] * x for x in t]
    #plt.plot(square(ToF_array), square(rho_array), 'ro')
    ##plt.plot(square(ToF_array), y)
    #line, = plt.plot(t, y, lw=2)
    #plt.show()
    return [sqrt(rhoregress[0]/rhoregress[1]), sqrt(rhoregress[1])]

# multiple fit for fermion
def fermionTemperature(tof, omegaAxial, omegaRadial, rx, ry, qx, qy):
    m = mLi

    p1 = (m*rx**2)/(2*kB)
    p2 = omegaAxial**2/(1+omegaAxial**2 * tof**2)
    p3 = 1/f(np.exp(qx))
    Tx = p1*p2*p3

    p1 = (m*ry**2)/(2*kB)
    p2 = omegaRadial**2/(1+omegaRadial**2 * tof**2)
    p3 = 1/f(np.exp(qy))
    Ty = p1*p2*p3

    return [Tx, Ty]


# temperature divided by fermionic temperature
def TOverTF(q):
    return (-6*polylog3(-np.exp(q)))**(-1./3.)

