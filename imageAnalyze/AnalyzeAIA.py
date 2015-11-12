import numpy as np
import os, sys, glob
import matplotlib as ml
import matplotlib.pyplot as plt
from imgFunc import *
from ImageUI import *
import matplotlib.pyplot as plt
from imagePlot import *



# def newFigure():
## some unit transformation

print "Reading Data..." + '\n'
# newest = max(glob.iglob('../*.aia'), key=os.path.getctime)
#print 'Latest image is ' + newest

# s = newest
## Override the image ###
#s = "C:\\ExperimentImages\\BEC_TOF_images\\in_trap.aia"
# s = "C:\\ExperimentImages\\BEC_TOF_images\\50ms.aia"


## imput time of flight in unit of ms
TimeOfFight = 110
s = '../BEC_TOF_images/' + str(TimeOfFight) + 'ms.aia'
s = '../BEC_TOF_images/in_trap.aia'
print 'Reading image' + s
ToF = TimeOfFight/1000

## read arguments ToF, Tmp freq FROM COMMEND LINE
#if len(sys.argv) == 5:
#    print 'Parameters are inputed!'
#    ToF = sys.argv[1]
#    OmegaX = sys.argv[2]
#    OmegaY = sys.argv[3]
#    OmegaZ = sys.argv[4]
#    print ("ToF is %s; Tmp of Freq is (%s, %s, %s)" %(ToF, OmegaX, OmegaY,OmegaZ)) 
#
#    f = open("parameters.txt", "w")
#    f.writelines(ToF+'\n')
#    f.writelines(OmegaX+'\n')
#    f.writelines(OmegaY+'\n')
#    f.writelines(OmegaZ)
#    f.close()
#
#else:
#    print 'No input parameters or Wrong format input. Read from file as last time'
#    f = open("parameters.txt", "r")
#    try:
#	ToF = f.readline()
#	OmegaX = f.readline()
#	OmegaY = f.readline()
#	OmegaZ = f.readline()
#    finally:
#	f.close()
#    print ("ToF is %s ms. Tmp of Freq is\n %s %s %s\n" %(ToF, OmegaX, OmegaY,OmegaZ)) 

#
#ToF = float(ToF) * 1E-3
#OmegaX = float(OmegaX) * 2 * np.pi
#OmegaY = float(OmegaY) * 2 * np.pi
#OmegaZ = float(OmegaZ) * 2 * np.pi





# Define Area of Interest in the form [upper left point, lower right point].
#AOI = [(200,300), (800,900)] 
AOI = [(0,0),(1024,1024)]
plotMin = 0.0
plotMax = 0.3


xLeft = AOI[0][0]
xRight = AOI[1][0]
yTop = AOI[0][1]
yBottom = AOI[1][1]

# Read Image
absorbImg = readData(s)

### Mapping from absorption to atom column density
atomImage = -np.log(absorbImg)

### Restrict to the are of interest
AOIImage = atomImage[yTop:yBottom,xLeft:xRight]



## Fit Image
print '\n'
print "Fitting Image..."

gVals = twoDGaussianFit(AOIImage)
gVals[0][0] += xLeft
gVals[0][1] += yTop
print '\n'
print 'Gaussian center at ('+ '%.1f'%(gVals[0][0]) + ', ' + '%.1f'%(gVals[0][1]) + ')'
print 'Gaussian sigmas are ('+ '%.1f'%(gVals[1][0]) + ', ' + '%.1f'%(gVals[1][1]) + ')'


pVals = twoDParbolicFit(AOIImage)
pVals[0][0] += xLeft
pVals[0][1] += yTop
print '\n'
print 'Parabolic center at ('+ '%.1f'%(pVals[0][0]) + ', ' + '%.1f'%(pVals[0][1]) + ')'
print 'Parabolic widths are ('+ '%.1f'%(pVals[1][0]*2) + ', ' + '%.1f'%(pVals[1][1]*2) + ')'


## view fit result
x = np.arange(xLeft, xRight, 1)
    
y = np.arange(yTop, yBottom, 1)
X,Y = np.meshgrid(x, y)
gaussianFitImage = gVals[2]*np.exp(-0.5*(((X-gVals[0][0])/gVals[1][0])**2+((Y-gVals[0][1])/gVals[1][1])**2))+gVals[3]
parabolicFitImage = pVals[2] * np.maximum(np.zeros((1024,1024)),  ( - (X - pVals[0][0])**2 + pVals[1][0] ** 2) ) * np.maximum(np.zeros((1024,1024)),   ( - (Y - pVals[0][1])**2 + pVals[1][1] ** 2) )+ pVals[3]


atomImagePlot([atomImage, gaussianFitImage, parabolicFitImage], ['original image','gaussian fit','parabolic fit'] )
print atomImage.max()
print gaussianFitImage.max()
print parabolicFitImage.max()
print atomImage.min()
print gaussianFitImage.min()
print parabolicFitImage.min()
print pVals
#atomImagePlot()
#atomImagePlot(parabolicFitImage)
# show UI Window

#app = wx.App()
#ui = ImageUI(None, title='Li/Na Image Analyze')
#
## ui.recalculate(AOIImage, ToF, sigmaX, sigmaY, OmegaX, OmegaY, crossSection)
#ui.TX.SetLabel('%f nK' %(Tmp[0] *1E9))
#ui.TY.SetLabel('%f nK' %(Tmp[1] *1E9))
#ui.atomN.SetLabel('%f' %AtomNumber)
#
#ui.LiOriginalFigure.drawFigure(AOIImage, plotMin, plotMax)
##ui.LiFittedFigure.drawFigure(fitImage, plotMin, plotMax)
#app.MainLoop()




#Science Calculation

#
### atom number
#print '\n'
#print 'extract quantities...'
#AtomNumber = atomNumber(AOIImage)
#print ('Atom Number is %.0f' %AtomNumber)
#print '\n'
#
### transfer units
#gSigmaX = gVals[1][0] * pixelToDistance
#gSigmaY = gVals[1][1] * pixelToDistance
#
### single image results ###
##ToF = 0.07
#omegaRadial = 250 * 2 *np.pi
#omegaAxial = 15 * 2 * np.pi
#
#
#
#rho0 = 0.00443E-3
#rho = pVals[1][1] * pixelToDistance
#z = pVals[1][0] * pixelToDistance
#
#
### Atom number using chemical potential
#mu = chemicalPotential(ToF, omegaRadial, rho)
##print 'mu%f'%(mu*1E30)
#U0 = effectiveInteraction(scatteringLength)
##print 'U0%f'%(U0*1E50)
#N = atomNumberFit(mu, omegaRadial, omegaAxial, U0)
#
#print 'Using parabolic fitting result and chemical potential, atom number is %.0f'%N
#
#
#
##print "in trap rho!!!!"
##print rho
#print 'Using gaussian fit results, calculate temperature for latset image,' +'\n' \
#+'with parameter time of flight = ' +'%.3f'%ToF+'s, trap frequency(radial) = '+'%.0f'%(omegaRadial/(2*np.pi))+'Hz and trap frequency(axial) = '+'%.0f'%(omegaAxial/(2*np.pi))+'Hz'
#
#Tmp = temperatureSingleGaussianFit(ToF, gSigmaX, gSigmaY, omegaAxial, omegaRadial, 'Na') 
#print '\n'
#print 'Temperature is (' + '%.1f'%(Tmp[0][0]* 1E9)  + '+-' + '%.1fnK)'%(Tmp[0][1]* 1E9)
#print ('Temperature from X is %.1f nK' %(Tmp[1][0] * 1E9) )
#print ('Temperature from Y is %.1f nK' %(Tmp[1][1] * 1E9) )
#
#
#print '\n'
#print 'Using parabolic fit results, calculate trapping frequency' + '\n'\
#+ 'with parameter time of flight = ' +'%.3f'%ToF+'s, rho0 = %.5fmm'%(rho0 * 1E3)
#
##print sqrt((rho/rho0)**2 - 1)/ToF
#omegaRadial = trapFrequencyRadial(ToF,  rho, rho0)
#print 'Trap frequency(radial) is %.1fHz'%(omegaRadial/(2*np.pi)) 
#
#omegaAxial = trapFrequencyAxial(ToF, z, rho0, omegaRadial)
#print 'Trap frequency(axial) is %.1fHz'%(omegaAxial/(2*np.pi))
#
#
#

# #PLOT FIGURES
# fig = plt.figure(1)
#
# plt.subplot(121)
# plt.imshow(AOIImage, vmin = plotMin, vmax = plotMax)
# plt.title("Li Original Img")
#
#
# plt.subplot(122)
# plt.imshow(fitImage, vmin = plotMin, vmax = plotMax)
# plt.title("Li Fitted Img")
#
#
# plt.draw() #need this to refresh image
# #plt.colorbar()

# plt.show()



#print(AOIImage)

