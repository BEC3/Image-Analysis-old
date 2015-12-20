#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys, glob
import wx, numpy
import matplotlib
matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from imagePlot import *
from imgFunc import *
from watchforchange import *


class ImageUI(wx.Frame):
  
    def __init__(self, parent, title):
        super(ImageUI, self).__init__(parent, title=title, 
            size=(1200, 700))
            
        self.InitUI()
        self.Centre()
        self.Show()
        # self.AOI = [(None,None),(None,None)]
        self.gVals = None
        self.pVals = None
        self.fVals = None
        self.AOIImage = None
        # self.Bind(wx.EVT_PAINT, self.OnPaint)

        self.qX = None
        self.qY = None

        self.filename = None
        self.Tmp = None
        self.data = None
        self.observer = Observer()
        self.observer.schedule(MyHandler(self.autoRun, self), path = self.imageFolderPath.GetValue())

    def InitUI(self):

    	panel = wx.Panel(self)
        font = wx.Font(18, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    	# font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
    	# font.SetPointSize(9)
######### file ############
        # set data path
        block1 = wx.StaticText(panel, pos=(10,20), size=(100, 25), label='Set Image Path')
        block1.SetFont(font)
        wx.StaticText(panel, pos=(20,55), size=(100, 25), label='Image Folder Path')
        # self.imageFolderPath = wx.TextCtrl(panel, pos=(20,75), size=(200,25), value="../BEC_TOF_images/")
        self.imageFolderPath = wx.TextCtrl(panel, pos=(20,75), size=(200,25), value="../data/fermion2/")
        wx.StaticText(panel, pos=(20,110), size=(100, 25), label='Image File Name')
        chooseFileButton = wx.Button(panel, pos=(20,130), size=(200,25), label = 'Choose File')
        chooseFileButton.Bind(wx.EVT_BUTTON, self.chooseFile)
        
        # self.imageFileName = wx.TextCtrl(panel, pos=(20,130), size=(200,25))
        self.filenameText = wx.TextCtrl(panel, pos=(20,170), size=(200, 25))
        #Fit Img button
        
        
        showImgButton = wx.Button(panel, pos=(10, 200), size = (90, 40), label = 'Fit Image')
        showImgButton.Bind(wx.EVT_BUTTON, self.showImg)
        autoButton = wx.Button(panel, pos=(10,520), size=(250,120), label = 'Run automatically')
        autoButton.Bind(wx.EVT_BUTTON, self.startAutoRun)
        autoButton.SetFont(font)

        #Fitting result
        block2 = wx.StaticText(panel, pos=(10,250), size=(100, 25), label='Fitting Result')
        block2.SetFont(font)
        wx.StaticText(panel, pos=(20,290), size=(100, 25), label='Gaussian Fit')
        wx.StaticText(panel, pos=(20,320), size=(100, 25), label='Center at')
        self.gCenter = wx.TextCtrl(panel, pos=(140,320), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(20,350), size=(100, 25), label='Sigma is')
        self.gSigma = wx.TextCtrl(panel, pos=(140,350), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(20,380), size=(100, 25), label='Paraboic Fit')
        wx.StaticText(panel, pos=(20,410), size=(100, 25), label='Center at')
        self.pCenter = wx.TextCtrl(panel, pos=(140,410), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(20,440), size=(100, 25), label='Half Width is')
        self.pWidth = wx.TextCtrl(panel, pos=(140,440), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(20,470), size=(100, 25), label='Fermion Fit')
        wx.StaticText(panel, pos=(20,500), size=(100, 25), label='(Rx, Ry)')
        self.fWidth = wx.TextCtrl(panel, pos=(140,500), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(20,530), size=(150, 25), label='q(=mu*beta)')
        self.fq = wx.TextCtrl(panel, pos=(140,530), size=(150, 25), value='', style=wx.TE_READONLY)


######### single shoot functions ############
        block3 = wx.StaticText(panel, pos=(380,20), size=(100, 25), label='Single Shot')
        block3.SetFont(font)
        block4 = wx.StaticText(panel, pos=(290,55), size=(100, 25), label='Input Values')
        block4.SetFont(font)

        #set input
        wx.StaticText(panel, pos=(300,90) , size=(50, 25), label = 'Time of Flight')
        self.tof = wx.TextCtrl(panel, pos=(420,90) , size=(50, 25), value='3')
        wx.StaticText(panel, pos=(470,90) , size=(30, 25), label = 'ms')

        wx.StaticText(panel, pos=(300,130) , size=(200, 30), label = 'Trapping Frequency')
        wx.StaticText(panel, pos=(320,160) , size=(50, 30), label = 'Axial')
        self.omegaAxial = wx.TextCtrl(panel, pos=(420,160) , size=(50, 30), value='200')
        wx.StaticText(panel, pos=(470,160) , size=(30, 25), label = '*2pi Hz')
        wx.StaticText(panel, pos=(320,190) , size=(50, 30), label = 'Radial')
        self.omegaRadial = wx.TextCtrl(panel, pos=(420,190) , size=(50, 30), value='250')
        wx.StaticText(panel, pos=(470,190) , size=(30, 25), label = '*2pi Hz')
        wx.StaticText(panel, pos=(300,230) , size=(50, 25), label = 'AOI: (x,y)->(x,y)')
        self.AOI1 = wx.TextCtrl(panel, pos=(410,230) , size=(30, 25), value='0')
        self.AOI2 = wx.TextCtrl(panel, pos=(440,230) , size=(30, 25), value='0')
        self.AOI3 = wx.TextCtrl(panel, pos=(470,230) , size=(40, 25), value='1024')
        self.AOI4 = wx.TextCtrl(panel, pos=(510,230) , size=(40, 25), value='1024')

        
        
        singleShootButton = wx.Button(panel, pos=(300, 260), size = (200, 25), label = 'Extract Scientific Properties')
        singleShootButton.Bind(wx.EVT_BUTTON, self.singleShoot)
        
        #atom number
        block5 = wx.StaticText(panel, pos=(290,295), size=(100, 25), label='Output Values')
        block5.SetFont(font)
        wx.StaticText(panel, pos=(300,320), size=(100, 25), label='Atom Number')
        wx.StaticText(panel, pos=(300,350), size=(100, 25), label='By integration')
        self.atomNumberInt = wx.TextCtrl(panel, pos=(450,350), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(300,380), size=(100, 25), label='By chemical potential')
        self.atomNumberChem = wx.TextCtrl(panel, pos=(450,380), size=(100, 25), value='', style=wx.TE_READONLY)
        #temperature
        wx.StaticText(panel, pos=(300,420), size=(100, 25), label='Temperature')
        self.gTemperature = wx.TextCtrl(panel, pos=(300,450), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(400,450), size=(30, 25), label='nK (from Gaussian fit)')
        self.fTemperature = wx.TextCtrl(panel, pos=(300,480), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(400,480), size=(30, 25), label='nK (from Fermion fit)')
        wx.StaticText(panel, pos=(300,510), size=(120, 25), label='T/T_F')
        self.tOverTF = wx.TextCtrl(panel, pos=(400,510), size=(150, 25), value='', style=wx.TE_READONLY)

######## save fitting result ##############
        saveButton = wx.Button(panel, pos=(300, 550), size = (200, 25), label = 'Save above results')
        saveButton.Bind(wx.EVT_BUTTON, self.saveResult)
        cleanButton = wx.Button(panel, pos=(300, 590), size = (200, 25), label = 'Remove all saved data')
        cleanButton.Bind(wx.EVT_BUTTON, self.cleanData)

######## multiple shoots functions ############
        block6 = wx.StaticText(panel, pos=(750,20), size=(100, 25), label='Multiple Shots')
    	block6.SetFont(font)

        block7 = wx.StaticText(panel, pos=(700,60), size=(100, 25), label='Input Values')
        block7.SetFont(font)
        readButton = wx.Button(panel, pos=(700, 100), size = (140, 25), label = 'Read saved data')
        readButton.Bind(wx.EVT_BUTTON, self.readData)
        self.dataReadedText = wx.TextCtrl(panel, pos=(850,100), size=(150, 25), value='input: 0', style=wx.TE_READONLY)


        block8 = wx.StaticText(panel, pos=(700,160), size=(100, 25), label='Output Values')
        block8.SetFont(font)
        fitTempButton = wx.Button(panel, pos=(700, 200), size = (200, 25), label = 'Fit temperature')
        fitTempButton.Bind(wx.EVT_BUTTON, self.fitTemp)
        self.fitTempText = wx.TextCtrl(panel, pos=(910,200), size=(110, 25), style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(1030,200), size=(30, 25), label='nK')
        wx.StaticText(panel, pos=(910,230) , size=(50, 30), label = 'rho0')
        self.fitrho0Text = wx.TextCtrl(panel, pos=(990,230), size=(50, 25), style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(1050,230) , size=(30, 25), label = 'um')
        fitTrapFreqButton = wx.Button(panel, pos=(700, 260), size = (200, 25), label = 'Fit Trap Frequency')
        fitTrapFreqButton.Bind(wx.EVT_BUTTON, self.fitTrapFreq)
        wx.StaticText(panel, pos=(910,260) , size=(50, 30), label = 'Axial')
        self.fitTrapAxialFreqText = wx.TextCtrl(panel, pos=(990,260), size=(50, 25), style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(1050,260) , size=(30, 25), label = '*2pi Hz')
        wx.StaticText(panel, pos=(910,290) , size=(50, 30), label = 'Radial')
        self.fitTrapRadialFreqText = wx.TextCtrl(panel, pos=(990,290), size=(50, 25), style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(1050,290) , size=(30, 25), label = '*2pi Hz')
        #

######## draw figures ############
        block9 = wx.StaticText(panel, pos=(750,350), size=(100, 25), label='Draw Figures')
        block9.SetFont(font)
        drawAtomNumberButton = wx.Button(panel, pos=(700, 390), size = (200, 25), label = 'Draw atom number figure')
        drawAtomNumberButton.Bind(wx.EVT_BUTTON, self.drawAtomNumber)

        block10 = wx.StaticText(panel, pos=(700,450), size=(100, 25), label='Recent Atom Number')
        block10.SetFont(font)
        self.recentAtomNumberInt = wx.TextCtrl(panel, pos=(700,490), size=(140, 180), style=wx.TE_READONLY | wx.TE_MULTILINE | wx.HSCROLL)
        self.recentAtomNumberChem = wx.TextCtrl(panel, pos=(850,490), size=(150, 180), style=wx.TE_READONLY | wx.TE_MULTILINE | wx.HSCROLL)
     

 #####################################################    

    def chooseFile(self, e):
        
        style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
        dialog = wx.FileDialog(None, 'Open', '', style=style)
        if dialog.ShowModal() == wx.ID_OK:
            self.filename = dialog.GetFilename()
            self.filenameText.SetValue(self.filename)
        else:
            self.filename = None

        dialog.Destroy()

    def test(self, e):
        print "Test for watchdog"

    def showImg(self, e):

        # Read Image
        path = self.imageFolderPath.GetValue()
        if not path:
            print "Wrong Folder!"
            return None
            # self.alert
        self.filename = self.filenameText.GetValue()

        
        if not self.filename:
            latest = max(glob.iglob(path + '*.aia'), key=os.path.getctime)
            s = latest
            s0 = s.split('/')
            self.filename = s0[-1]
        else:
            s = path + self.filename



        
        plotMin = 0.0
        plotMax = 0.3


        xLeft = int(self.AOI1.GetValue())
        # print xLeft
        xRight = int(self.AOI3.GetValue())
        # print xRight
        yTop = int(self.AOI2.GetValue())
        # print yTop
        yBottom = int(self.AOI4.GetValue())
        # print yBottom

        absorbImg = readData(s)


        ### Mapping from absorption to atom column density
        atomImage = -np.log(absorbImg)
        # print sort(np.array(atomImage))

        ### Restrict to the are of interest
        self.AOIImage = atomImage[yTop:yBottom,xLeft:xRight]
        ## view fit result

        print "Gaussian Fit"
        self.gVals = twoDGaussianFit(self.AOIImage)
        self.gVals[0][0] += xLeft
        self.gVals[0][1] += yTop
        self.offset = self.gVals[3]

        # print "parabolic Fit"
        # self.pVals = twoDParbolicFit(self.AOIImage)
        # self.pVals[0][0] += xLeft
        # self.pVals[0][1] += yTop
        # self.offset = self.pVals[3]

        # print "partly Condensate"
        # self.doubleVals = partlyCondensateFit(self.AOIImage)
        # self.doubleVals[0][0] += xLeft
        # self.doubleVals[0][1] += yTop
        # self.doubleoffset = self.doubleVals[5]

        print "fermion Fit"
        self.fVals = fermionFit(self.AOIImage)
        self.fVals[0][0] += xLeft
        self.fVals[0][1] += yTop
        self.foffset = self.fVals[4]
        q = 0.5* (self.fVals[3][0] + self.fVals[3][1])


        self.gCenter.SetValue('( %.0f'%self.gVals[0][0] + ' , %.0f )'%self.gVals[0][1])
        self.gSigma.SetValue('( %.0f'%self.gVals[1][0] + ' , %.0f )'%self.gVals[1][1])
        # self.pCenter.SetValue('( %.0f'%self.pVals[0][0] + ' , %.0f )'%self.pVals[0][1])
        # self.pWidth.SetValue('( %.0f'%(self.pVals[1][0]) + ' , %.0f )'%(self.pVals[1][1]))
        self.fWidth.SetValue('( %.0f'%(self.fVals[1][0]) + ' , %.0f )'%(self.fVals[1][1]))
        self.fq.SetValue('( %.2f'%(self.fVals[3][0]) + ' , ' + '%.2f )'%(self.fVals[3][1]))


        x = np.arange(xLeft, xRight, 1)
        y = np.arange(yTop, yBottom, 1)
        X,Y = np.meshgrid(x, y)

        print "redraw gaussian image"
        gaussianFitImage = self.gVals[2]*np.exp(-0.5*(((X-self.gVals[0][0])/self.gVals[1][0])**2+((Y-self.gVals[0][1])/self.gVals[1][1])**2))+self.gVals[3]
        
        # print "redraw parabolic image"
        # parabolicFitImage = np.maximum(np.zeros((1024,1024)), \
        # self.pVals[2] * (1 - ((X - self.pVals[0][0]) / self.pVals[1][0]) ** 2 - ((Y - self.pVals[0][1]) / self.pVals[1][1]) ** 2)  )+\
        # self.pVals[3]

        # print "redraw partly condensate image"
        # partlyFitImage = self.doubleVals[4] * \
        # np.exp(-0.5*(((X-self.doubleVals[0][0])/self.doubleVals[2][0])**2+((Y-self.doubleVals[0][1])/self.doubleVals[2][1])**2)) + \
        # self.doubleVals[3] * \
        # np.maximum(np.zeros((1024,1024)),  (1 - ((X - self.doubleVals[0][0]) / self.doubleVals[1][0]) ** 2 - ((Y - self.doubleVals[0][1]) / self.doubleVals[1][1])  ** 2 )  ) ** 2+ self.doubleoffset


        print "redraw fermion image"
        numerator_inside = - np.exp(q-((X-self.fVals[0][0])**2/self.fVals[1][0]**2+(Y-self.fVals[0][1])**2/self.fVals[1][1]**2)*f(np.exp(q)))
        l = len(numerator_inside)
        s = []
        for i in range(l):
            x = self.fVals[2] * polylog2(np.array(numerator_inside[i]))/polylog2(- np.exp(q))  + self.fVals[4]
            s.append(list(x))
        
        fermionFitImage = s
       

        

        # atomImagePlot([atomImage, gaussianFitImage,parabolicFitImage], ['original image','gaussian fit','parabolic fit'] )
        print "plot images"
        atomImagePlot([atomImage, gaussianFitImage, fermionFitImage], ['original image', 'gaussianFitImage', 'fermion fit'] )
        # atomImagePlot([atomImage, gaussianFitImage, parabolicFitImage, partlyFitImage], ['original image','gaussian fit','parabolic fit', 'partly BEC fit'] )

    def startAutoRun(self, e):
        print "Begin to watch new files"
        self.observer.start()


    def autoRun(self, e):
        self.showImg(e)
        self.singleShoot(e)
        # self.saveResult(e)
        self.recentAtomNumberInt.AppendText(self.atomNumberInt.GetValue() + '\n')
        # self.recentAtomNumberChem.AppendText(self.atomNumberChem.GetValue() + '\n')



    def singleShoot(self, e):
            ## Atom number using chemical potential
        ToF = float(self.tof.GetValue()) / 1000
        omegaAxial = float(self.omegaAxial.GetValue())
        omegaRadial = float(self.omegaRadial.GetValue())

        # mu = chemicalPotential(ToF, omegaRadial * 2 * np.pi, self.pVals[1][1]* pixelToDistance)
        #print 'mu%f'%(mu*1E30)
        # U0 = effectiveInteraction(scatteringLength)
        #print 'U0%f'%(U0*1E50)
        # N_chem = atomNumberFit(mu, omegaRadial * 2 * np.pi, omegaAxial * 2 * np.pi, U0)
        
        # self.atomNumberChem.SetValue(str("%.0f" % N_chem))

        N_int = atomNumber(self.AOIImage, self.offset)
        self.atomNumberInt.SetValue(str("%.0f" % N_int))

        Rx = self.fVals[1][0] * pixelToDistance
        Ry = self.fVals[1][1] * pixelToDistance
        self.qX = self.fVals[3][0]
        self.qY = self.fVals[3][1]

        tovertfx = TOverTF(self.qX)
        tovertfy = TOverTF(self.qY)

        gSigmaX = self.gVals[1][0] * pixelToDistance
        gSigmaY = self.gVals[1][1] * pixelToDistance

        self.gTmp = temperatureSingleGaussianFit(ToF, gSigmaX, gSigmaY, omegaAxial, omegaRadial, 'Li') 
        self.fTmp = fermionTemperature(ToF, omegaAxial, omegaRadial, Rx, Ry, self.qX, self.qY)
        self.gTemperature.SetValue(str('( %.0f' % (self.gTmp[0]*1E9)) + ' , ' + '%.0f )' % (self.gTmp[1]*1E9))
        self.fTemperature.SetValue(str('( %.0f' % (self.fTmp[0]*1E9)) + ' , ' + '%.0f )' % (self.fTmp[1]*1E9))
        self.tOverTF.SetValue(str('( %.3f' % tovertfx + ' , ' + '%.3f )'%tovertfy))

    def saveResult(self, e):
        f = open("../data.txt", "a")
        f.writelines(self.filename + ' , ' + self.tof.GetValue() + ' , ' + self.omegaAxial.GetValue() + ' , ' + self.omegaRadial.GetValue() + ' , '+ str(self.pVals[0][0]) + ' , ' + str(self.pVals[0][1]) + ' , ' + str(self.gVals[1][0]) + ' , ' + str(self.gVals[1][1]) + ' , ' + str(self.pVals[1][0]) + ' , ' + str(self.pVals[1][1]) + ' , '  + self.atomNumberInt.GetValue() + ' , ' + self.atomNumberChem.GetValue() + ' , ' + str(self.gTmp[0]*1E9) + ' , '  +  str(self.gTmp[1]*1E9) + ' , '  + str(self.fVals[1][0]) + ' , '+ str(self.fVals[1][1]) + ' , ' + str(self.fVals[3][0]) + ' , '+ str(self.fVals[3][1])+ '\n') 
        
        f.close()

    def cleanData(self, e):
        f = open("../data.txt", "w")
        f.close()

    def readData(self, e):
        f = open("../data.txt", "r")

        self.data = f.readlines()
        
        self.dataReadedText.SetValue("input: %i"%len(self.data))
        for i in range(len(self.data)):
            self.data[i] = self.data[i].split(' , ')
        
    
        f.close()

    def fitTemp(self, e):
        tofList = []
        fRXList = []
        fRYList = []
        qXList = []
        qYList = []
        n = len(self.data)
        for i in self.data:
            tofList.append(float(i[1])/1000.)
            fRXList.append(float(i[14]) * pixelToDistance)
            fRYList.append(float(i[15]) * pixelToDistance)
            qXList.append(float(i[16]))
            qYList.append(float(i[17]))

        t = fitTemperature(tofList, fRXList, fRYList, qXList, qYList)
        self.fitTempText.SetValue('(%.1f' %(t[0]*1E9) + ' , ' + '%.1f )' %(t[1]*1E9))


    def fitTrapFreq(self, e):
        tofList = []
        fRXList = []
        fRYList = []
        qXList = []
        qYList = []
        n = len(self.data)
        for i in self.data:
            tofList.append(float(i[1])/1000.)
            fRXList.append(float(i[14]) * pixelToDistance)
            fRYList.append(float(i[15]) * pixelToDistance)
            qXList.append(float(i[16]))
            qYList.append(float(i[17]))

        t = fitTrapFrequency(tofList, fRYList, fRYList, qXList, qYList)
        # tofList = []
        # pHalfWidthRadialList = []
        # pHalfWidthAxialList = []
        # initTrapFreqRadial = float(self.omegaRadial.GetValue()) * 2 * np.pi
        # initTrapFreqAxial = float(self.omegaAxial.GetValue()) * 2 * np.pi
        # for i in self.data:
        #     tofList.append(float(i[1])/1000.)
        #     pHalfWidthRadialList.append(float(i[9]) * pixelToDistance)
        #     pHalfWidthAxialList.append(float(i[8]) * pixelToDistance)
        # t = fitTrapFrequency(tofList, pHalfWidthRadialList, pHalfWidthAxialList, initTrapFreqRadial, initTrapFreqAxial)
        self.fitTrapAxialFreqText.SetValue(str('%.1f' % (t[1]/(2*np.pi))))
        self.fitTrapRadialFreqText.SetValue(str('%.1f' % (t[0]/(2*np.pi))))
        # self.fitrho0Text.SetValue(str('%.1f' % (t[2]*1E6)))

    def drawAtomNumber(self, e):
        atomNumberI = []
        atomNumberC = []
        n = len(self.data)
        for i in self.data:
            atomNumberI.append(int(i[10]))
            atomNumberC.append(int(i[11]))
        atomNumberPlot(n, atomNumberI, atomNumberC)

            
	
class FigurePanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent,  size=(600, 600))
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)

        # self.axes = self.figure.add_subplot(111)
        # self.canvas = FigureCanvas(self, -1, self.figure)
        # self.sizer = wx.BoxSizer(wx.VERTICAL)
        # self.sizer.Add(self.canvas, 1, wx.ALL)
        # print self.sizer
        # self.SetSizer(self.sizer)
        # self.Fit()
        self.canvas = FigureCanvas(self, -1, self.figure)

    def drawFigure(self, data, plotMin, plotMax):
        # t = arange(0.0, 3.0, 0.01)
        # s = sin(2*pi*t)
        # self.axes.plot(t, s)
        self.axes.imshow(data, vmin = plotMin, vmax = plotMax)



        # self.draw()

# plt.subplot(121)
# plt.imshow(AOIImage, vmin = plotMin, vmax = plotMax)
# plt.title("Li Original Img")





if __name__ == '__main__':
  
    app = wx.App()
    ui = ImageUI(None, title='Image Analyze')
    # ui.LiOriginalFigure.drawFigure()
    # ui.LiOriginalFigure.draw()
    app.MainLoop()