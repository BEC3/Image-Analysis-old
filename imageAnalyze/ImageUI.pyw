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
from localPath import *
from fitTool import *


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
        self.AOI = None
        # self.Bind(wx.EVT_PAINT, self.OnPaint)

        self.q = None
        

        self.filename = None
        self.Tmp = None
        self.data = None
        self.observer = Observer()
        self.observer.schedule(MyHandler(self.autoRun, self), path = self.imageFolderPath.GetValue())

        self.fitMethodFermion.SetValue(True)
        self.FermionFitChosen(e)

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
        self.imageFolderPath = wx.TextCtrl(panel, pos=(20,75), size=(200,25), value=LOCAL_PATH)
        wx.StaticText(panel, pos=(20,110), size=(100, 25), label='Image File Name')
        chooseFileButton = wx.Button(panel, pos=(20,135), size=(200,25), label = 'Choose File')
        chooseFileButton.Bind(wx.EVT_BUTTON, self.chooseFile)
        
        # self.imageFileName = wx.TextCtrl(panel, pos=(20,130), size=(200,25))
        self.filenameText = wx.TextCtrl(panel, pos=(20,170), size=(200, 25))
        #Fit Img button
        
        self.fitMethodFermion = wx.RadioButton(panel, label="Fermion", pos=(10,240), size = (80,20))
        self.fitMethodBoson = wx.RadioButton(panel, label="Boson", pos=(10,265), size=(80,20))
        self.Bind(wx.EVT_RADIOBUTTON, self.FermionFitChosen, id=self.fitMethodFermion.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.BosonFitChosen, id=self.fitMethodBoson.GetId())

        wx.StaticText(panel, pos=(20,210) , size=(100, 25), label = 'AOI: (x,y)->(x,y)')
        self.AOI1 = wx.TextCtrl(panel, pos=(130,210) , size=(35, 25), value='0')
        self.AOI2 = wx.TextCtrl(panel, pos=(165,210) , size=(35, 25), value='0')
        self.AOI3 = wx.TextCtrl(panel, pos=(200,210) , size=(40, 25), value='1024')
        self.AOI4 = wx.TextCtrl(panel, pos=(240,210) , size=(40, 25), value='1024')

        

        showImgButton = wx.Button(panel, pos=(140, 240), size = (90, 40), label = 'Fit Image')
        showImgButton.Bind(wx.EVT_BUTTON, self.showImg)
        autoButton = wx.Button(panel, pos=(10,560), size=(250,60), label = 'Run automatically')
        autoButton.Bind(wx.EVT_BUTTON, self.startAutoRun)
        autoButton.SetFont(font)

        #Fitting result
        block2 = wx.StaticText(panel, pos=(10,300), size=(100, 25), label='Fitting Result')
        block2.SetFont(font)
        wx.StaticText(panel, pos=(20,340), size=(100, 25), label='Gaussian Fit')
        wx.StaticText(panel, pos=(20,370), size=(100, 25), label='Center at')
        self.gCenter = wx.TextCtrl(panel, pos=(140,370), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(20,400), size=(100, 25), label='Sigma is')
        self.gSigma = wx.TextCtrl(panel, pos=(140,400), size=(100, 25), value='', style=wx.TE_READONLY)
        
        self.pLabel = wx.StaticText(panel, pos=(20,430), size=(100, 25), label='Paraboic Fit')
        self.pText1 = wx.StaticText(panel, pos=(20,460), size=(100, 25), label='Center at')
        self.pCenter = wx.TextCtrl(panel, pos=(140,460), size=(100, 25), value='', style=wx.TE_READONLY)
        self.pText2 = wx.StaticText(panel, pos=(20,490), size=(100, 25), label='Half Width is')
        self.pWidth = wx.TextCtrl(panel, pos=(140,490), size=(100, 25), value='', style=wx.TE_READONLY)
        
        self.fLabel = wx.StaticText(panel, pos=(20,430), size=(100, 25), label='Fermion Fit')
        self.fText1 = wx.StaticText(panel, pos=(20,460), size=(100, 25), label='(Rx, Ry)')
        self.fWidth = wx.TextCtrl(panel, pos=(140,460), size=(100, 25), value='', style=wx.TE_READONLY)
        self.fText2 = wx.StaticText(panel, pos=(20,490), size=(150, 25), label='q(=mu*beta)')
        self.fq = wx.TextCtrl(panel, pos=(140,490), size=(150, 25), value='', style=wx.TE_READONLY)


######### single shoot functions ############
        block3 = wx.StaticText(panel, pos=(420,20), size=(100, 25), label='Single Shot')
        block3.SetFont(font)
        block4 = wx.StaticText(panel, pos=(420,55), size=(100, 25), label='Set parameters')
        block4.SetFont(font)

        #set input
        wx.StaticText(panel, pos=(420,90) , size=(50, 25), label = 'Time of Flight')
        self.tof = wx.TextCtrl(panel, pos=(540,90) , size=(50, 25), value='3')
        wx.StaticText(panel, pos=(590,90) , size=(30, 25), label = 'ms')

        wx.StaticText(panel, pos=(420,130) , size=(200, 30), label = 'Trapping Frequency')
        wx.StaticText(panel, pos=(440,160) , size=(50, 30), label = 'Axial')
        self.omegaAxial = wx.TextCtrl(panel, pos=(540,160) , size=(50, 30), value='200')
        wx.StaticText(panel, pos=(590,160) , size=(60, 25), label = '*2pi Hz')
        wx.StaticText(panel, pos=(440,190) , size=(50, 30), label = 'Radial')
        self.omegaRadial = wx.TextCtrl(panel, pos=(540,190) , size=(50, 30), value='250')
        wx.StaticText(panel, pos=(590,190) , size=(60, 25), label = '*2pi Hz')
        
        
        singleShootButton = wx.Button(panel, pos=(420, 230), size = (200, 25), label = 'Extract Scientific Properties')
        singleShootButton.Bind(wx.EVT_BUTTON, self.singleShoot)
        
        #atom number
        block5 = wx.StaticText(panel, pos=(420,295), size=(100, 25), label='Output Values')
        block5.SetFont(font)
        wx.StaticText(panel, pos=(420,320), size=(100, 25), label='Atom Number')
        wx.StaticText(panel, pos=(420,350), size=(100, 25), label='By integration')
        self.atomNumberInt = wx.TextCtrl(panel, pos=(570,350), size=(100, 25), value='', style=wx.TE_READONLY)
        self.atomNumberChemLabel = wx.StaticText(panel, pos=(420,380), size=(100, 25), label='By chemical potential')
        self.atomNumberChem = wx.TextCtrl(panel, pos=(570,380), size=(100, 25), value='', style=wx.TE_READONLY)
        #temperature
        wx.StaticText(panel, pos=(420,420), size=(100, 25), label='Temperature')
        self.gTemperature = wx.TextCtrl(panel, pos=(420,450), size=(100, 25), value='', style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(520,450), size=(100, 25), label='nK (from Gaussian fit)')
        # self.fTemperature = wx.TextCtrl(panel, pos=(420,480), size=(100, 25), value='', style=wx.TE_READONLY)
        # self.fTempLabel = wx.StaticText(panel, pos=(520,480), size=(100, 25), label='nK (from Fermion fit)')
        self.tOverTFLabel = wx.StaticText(panel, pos=(420,510), size=(120, 25), label='T/T_F')
        self.tOverTF = wx.TextCtrl(panel, pos=(520,510), size=(150, 25), value='', style=wx.TE_READONLY)

######## save fitting result ##############
        self.saveFermionButton = wx.Button(panel, pos=(420, 580), size = (160, 25), label = 'Save fermion results')
        self.saveFermionButton.Bind(wx.EVT_BUTTON, self.saveFermionResult)
        self.saveBosonButton = wx.Button(panel, pos=(420, 580), size = (160, 25), label = 'Save boson results')
        self.saveBosonButton.Bind(wx.EVT_BUTTON, self.saveBosonResult)
        cleanButton = wx.Button(panel, pos=(420, 610), size = (200, 25), label = 'Remove all saved data')
        cleanButton.Bind(wx.EVT_BUTTON, self.cleanData)

######## multiple shoots functions ############
        block6 = wx.StaticText(panel, pos=(820,20), size=(100, 25), label='Multiple Shots')
    	block6.SetFont(font)

        block7 = wx.StaticText(panel, pos=(820,60), size=(100, 25), label='Read Data')
        block7.SetFont(font)
        readButton = wx.Button(panel, pos=(820, 100), size = (140, 25), label = 'Read saved data')
        readButton.Bind(wx.EVT_BUTTON, self.readData)
        self.dataReadedText = wx.TextCtrl(panel, pos=(970,100), size=(150, 25), value='input: 0', style=wx.TE_READONLY)


        block8 = wx.StaticText(panel, pos=(820,160), size=(100, 25), label='Output Values')
        block8.SetFont(font)
        fitListButton = wx.Button(panel, pos=(820, 200), size = (200, 25), label = 'List Data Fit')
        fitListButton.Bind(wx.EVT_BUTTON, self.fitListData)
        wx.StaticText(panel, pos=(850,230) , size=(50, 30), label = 'Temperature')
        self.fitTempText = wx.TextCtrl(panel, pos=(940,230), size=(110, 25), style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(1060,230), size=(30, 25), label='nK')
        # wx.StaticText(panel, pos=(910,230) , size=(50, 30), label = 'rho0')
        # self.fitrho0Text = wx.TextCtrl(panel, pos=(990,230), size=(50, 25), style=wx.TE_READONLY)
        # wx.StaticText(panel, pos=(1050,230) , size=(30, 25), label = 'um')
        # fitTrapFreqButton = wx.Button(panel, pos=(700, 260), size = (200, 25), label = 'List Data Fit')
        # fitTrapFreqButton.Bind(wx.EVT_BUTTON, self.fitListData)
        wx.StaticText(panel, pos=(850,260) , size=(50, 30), label = 'Axial')
        self.fitTrapAxialFreqText = wx.TextCtrl(panel, pos=(940,260), size=(50, 25), style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(1000,260) , size=(60, 25), label = '*2pi Hz')
        wx.StaticText(panel, pos=(850,290) , size=(50, 30), label = 'Radial')
        self.fitTrapRadialFreqText = wx.TextCtrl(panel, pos=(940,290), size=(50, 25), style=wx.TE_READONLY)
        wx.StaticText(panel, pos=(1000,290) , size=(60, 25), label = '*2pi Hz')
        #

######## draw figures ############
        block9 = wx.StaticText(panel, pos=(820,350), size=(100, 25), label='Draw Figures')
        block9.SetFont(font)
        drawAtomNumberButton = wx.Button(panel, pos=(820, 390), size = (200, 25), label = 'Draw atom number figure')
        drawAtomNumberButton.Bind(wx.EVT_BUTTON, self.drawAtomNumber)

        block10 = wx.StaticText(panel, pos=(820,450), size=(200, 25), label='Recent Atom Number List')
        block10.SetFont(font)
        self.recentAtomNumberInt = wx.TextCtrl(panel, pos=(820,490), size=(140, 180), style=wx.TE_READONLY | wx.TE_MULTILINE | wx.HSCROLL)
        # self.recentAtomNumberChem = wx.TextCtrl(panel, pos=(850,490), size=(150, 180), style=wx.TE_READONLY | wx.TE_MULTILINE | wx.HSCROLL)
     

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

    def FermionFitChosen(self, e):
        print "Fermion Fit"
        self.pLabel.Hide()
        self.pText1.Hide()
        self.pCenter.Hide()
        self.pText2.Hide()
        self.pWidth.Hide()
        self.atomNumberChemLabel.Hide()
        self.atomNumberChem.Hide()
        self.fLabel.Show()
        self.fText1.Show()
        self.fWidth.Show()
        self.fText2.Show()
        self.fq.Show()
        # self.fTempLabel.Show()
        # self.fTemperature.Show()
        self.tOverTF.Show()
        self.tOverTFLabel.Show()
        self.saveFermionButton.Show()
        self.saveBosonButton.Hide()

    def BosonFitChosen(self, e):

        print "Boson Fit"
        self.pLabel.Show()
        self.pText1.Show()
        self.pCenter.Show()
        self.pText2.Show()
        self.pWidth.Show()
        self.atomNumberChemLabel.Show()
        self.atomNumberChem.Show()
        self.fLabel.Hide()
        self.fText1.Hide()
        self.fWidth.Hide()
        self.fText2.Hide()
        self.fq.Hide()
        # self.fTempLabel.Hide()
        # self.fTemperature.Hide()
        self.tOverTF.Hide()
        self.tOverTFLabel.Hide()
        self.saveFermionButton.Hide()
        self.saveBosonButton.Show()

    def showImg(self, e):

        print "Begin to fit..."

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
        self.AOI = [[xLeft,yTop],[xRight,yBottom]]

        absorbImg = readData(s)


        ### Mapping from absorption to atom column density
        atomImage = -np.log(absorbImg)
        # print sort(np.array(atomImage))

        ### Restrict to the are of interest
        self.AOIImage = atomImage[yTop:yBottom,xLeft:xRight]
        ## view fit result

        print "Gaussian Fit"
        self.gaussionParams = fitData(self.AOIImage, gaussionDistribution)
        """x0, y0, a, b, amplitude, offset"""
        self.gaussionParams[0] += xLeft
        self.gaussionParams[1] += yTop
        self.offset = self.gaussionParams[5]

       
        x = np.arange(xLeft, xRight, 1)
        y = np.arange(yTop, yBottom, 1)
        
        x0, y0, a, b, amplitude, offset = self.gaussionParams
        print "redraw gaussian image"
        # gaussianFitImage = self.gVals[2]*np.exp(-0.5*(((X-self.gVals[0][0])/self.gVals[1][0])**2+((Y-self.gVals[0][1])/self.gVals[1][1])**2))+self.gVals[3]
        size = np.shape(self.AOIImage)
        coordinates = np.meshgrid(x, y)
        gaussianFitImage = gaussionDistribution(coordinates, x0, y0, a, b, amplitude, offset).reshape(1024,1024)
        
        # atomImagePlot([atomImage, gaussianFitImage], ['original image', 'gaussianFitImage'] )
        self.gCenter.SetValue('( %.0f'%x0 + ' , %.0f )'%y0)
        self.gSigma.SetValue('( %.0f'%a + ' , %.0f )'%b)
        


        if self.fitMethodFermion.GetValue():
            print "fermion Fit"
            self.fermionParams = fitData(self.AOIImage, fermionDistribution)
            """x0, y0, a, b, amplitude, offset, q"""
            x0, y0, a, b, amplitude, offset, q = self.fermionParams
        
            self.fermionParams[0] += xLeft
            self.fermionParams[1] += yTop
            self.foffset = self.fermionParams[5]
            self.fWidth.SetValue('( %.0f'%(self.fermionParams[2]) + ' , %.0f )'%(self.fermionParams[3]))
            self.fq.SetValue('%.2f'%(self.fermionParams[6]))
        
            print "redraw fermion image"
            fermionFitImage = fermionDistribution(coordinates, x0, y0, a, b, amplitude, offset, q).reshape(1024,1024)
       
            atomImagePlot([atomImage, gaussianFitImage, fermionFitImage], ['original image', 'gaussianFitImage', 'fermion fit'] )
        
        elif self.fitMethodBoson.GetValue():
            print "boson Fit"
        #     self.pVals = twoDParbolicFit(self.AOIImage)
        #     self.pVals[0][0] += xLeft
        #     self.pVals[0][1] += yTop
        #     self.offset = self.pVals[3]

        #     # print "partly Condensate"
        #     # self.doubleVals = partlyCondensateFit(self.AOIImage)
        #     # self.doubleVals[0][0] += xLeft
        #     # self.doubleVals[0][1] += yTop
        #     # self.doubleoffset = self.doubleVals[5]

        #     print "redraw parabolic image"
        #     parabolicFitImage = np.maximum(np.zeros((1024,1024)), \
        #     self.pVals[2] * (1 - ((X - self.pVals[0][0]) / self.pVals[1][0]) ** 2 - ((Y - self.pVals[0][1]) / self.pVals[1][1]) ** 2)  )+\
        #     self.pVals[3]

            print "redraw partly condensate image"
        #     # partlyFitImage = self.doubleVals[4] * \
        #     # np.exp(-0.5*(((X-self.doubleVals[0][0])/self.doubleVals[2][0])**2+((Y-self.doubleVals[0][1])/self.doubleVals[2][1])**2)) + \
        #     # self.doubleVals[3] * \
        #     # np.maximum(np.zeros((1024,1024)),  (1 - ((X - self.doubleVals[0][0]) / self.doubleVals[1][0]) ** 2 - ((Y - self.doubleVals[0][1]) / self.doubleVals[1][1])  ** 2 )  ) ** 2+ self.doubleoffset
        #     atomImagePlot([atomImage, gaussianFitImage, parabolicFitImage], ['original image','gaussian fit', 'parabolic fit'])
        #     # atomImagePlot([atomImage, gaussianFitImage, parabolicFitImage, partlyFitImage], ['original image','gaussian fit','parabolic fit', 'partly BEC fit'] )

       


        # # self.pCenter.SetValue('( %.0f'%self.pVals[0][0] + ' , %.0f )'%self.pVals[0][1])
        # # self.pWidth.SetValue('( %.0f'%(self.pVals[1][0]) + ' , %.0f )'%(self.pVals[1][1]))
       


      

        

        # atomImagePlot([atomImage, gaussianFitImage,parabolicFitImage], ['original image','gaussian fit','parabolic fit'] )
        
        

    def startAutoRun(self, e):
        print "Begin to watch new files"
        self.observer.start()


    def autoRun(self, e):
        self.showImg(e)
        self.singleShoot(e)
        if self.fitMethodFermion.GetValue():
            self.saveFermionResult(e)
        elif self.fitMethodBoson.GetValue():
            self.saveBosonResult(e)
        self.recentAtomNumberInt.AppendText(self.atomNumberInt.GetValue() + '\n')
        # self.recentAtomNumberChem.AppendText(self.atomNumberChem.GetValue() + '\n')



    def singleShoot(self, e):
            ## Atom number using chemical potential

        ToF = float(self.tof.GetValue()) / 1000
        omegaAxial = float(self.omegaAxial.GetValue()) * np.pi * 2
        omegaRadial = float(self.omegaRadial.GetValue()) * np.pi * 2

        atom = ""
        if self.fitMethodBoson.GetValue():
            atom = "Na"
        elif self.fitMethodFermion.GetValue():
            atom = "Li"

        # mu = chemicalPotential(ToF, omegaRadial * 2 * np.pi, self.pVals[1][1]* pixelToDistance)
        #print 'mu%f'%(mu*1E30)
        # U0 = effectiveInteraction(scatteringLength)
        #print 'U0%f'%(U0*1E50)
        # N_chem = atomNumberFit(mu, omegaRadial * 2 * np.pi, omegaAxial * 2 * np.pi, U0)
        
        # self.atomNumberChem.SetValue(str("%.0f" % N_chem))
        N_int = 0

        N_int = atomNumber(self.AOIImage, self.offset)
        self.atomNumberInt.SetValue(str("%.0f" % N_int))

        if self.fitMethodFermion.GetValue():
            Rx = self.fermionParams[2] * pixelToDistance
            Ry = self.fermionParams[3] * pixelToDistance
            self.q = self.fermionParams[6]

            tovertf = TOverTF(self.q)
            

            # self.fTmp = fermionTemperature(ToF, omegaAxial, omegaRadial, Rx, Ry, self.qX, self.qY)
            # self.fTemperature.SetValue(str('( %.0f' % (self.fTmp[0]*1E9)) + ' , ' + '%.0f )' % (self.fTmp[1]*1E9))
            self.tOverTF.SetValue(str('%.3f' % tovertf ))


        gSigmaX = self.gaussionParams[2] * pixelToDistance
        gSigmaY = self.gaussionParams[3] * pixelToDistance

        self.gTmp = temperatureSingleGaussianFit(ToF, gSigmaX, gSigmaY, omegaAxial, omegaRadial, atom) 
        
        self.gTemperature.SetValue(str('( %.0f' % (self.gTmp[0]*1E9)) + ' , ' + '%.0f )' % (self.gTmp[1]*1E9))
        
    def saveBosonResult(self, e):
        f = open("../boson_data.txt", "a")
        f.writelines(self.filename + ' , ' + self.tof.GetValue() + ' , '\
         # + self.omegaAxial.GetValue() + ' , ' + self.omegaRadial.GetValue() + ' , '\
         # + str(self.gVals[0][0]) + ' , ' + str(self.gVals[0][1]) + ' , ' \
         + self.atomNumberInt.GetValue() + ' , ' \
         + str(self.gaussionParams[2]) + ' , ' + str(self.gaussionParams[3]) + ' , ' \
         # + str(self.pVals[1][0]) + ' , ' + str(self.pVals[1][1])  \
            + '\n') 
        
        f.close()

    def saveFermionResult(self, e):
        f = open("../fermion_data.txt", "a")
        
        f.writelines(self.filename + ' , ' + self.tof.GetValue() + ' , '\
         # + self.omegaAxial.GetValue() + ' , ' + self.omegaRadial.GetValue() + ' , '\
         # + str(self.gVals[0][0]) + ' , ' + str(self.gVals[0][1]) + ' , ' \
         + self.atomNumberInt.GetValue() + ' , ' \
         + str(self.fermionParams[2]) + ' , ' + str(self.fermionParams[3]) + ' , ' \
         + str(self.fermionParams[6]) + '\n') 
        
        f.close()

    def cleanData(self, e):
        if self.fitMethodFermion.GetValue():
            f = open("../fermion_data.txt", "w")
        elif self.fitMethodBoson.GetValue():
            f = open("../boson_data.txt", "w")
        f.close()

    def readData(self, e):
        if self.fitMethodFermion.GetValue():
            f = open("../fermion_data.txt", "r")
            self.data = f.readlines()
        elif self.fitMethodBoson.GetValue():
            f = open("../boson_data.txt", "r")
            self.data = f.readlines()
        # f = open("../data.txt", "r")

        
        
        self.dataReadedText.SetValue("input: %i"%len(self.data))
        
        for i in range(len(self.data)):
            self.data[i] = self.data[i].split(' , ')
        
    
        f.close()

    def fitListData(self, e):
        tofList = []
        RXList = []
        RYList = []
       
        atom = ""

        n = len(self.data)
        for i in self.data:
            tofList.append(float(i[1])/1000.)
            RXList.append(float(i[3]) * pixelToDistance)
            RYList.append(float(i[4]) * pixelToDistance)
           
        if self.fitMethodBoson.GetValue():
            atom = "Na"
        elif self.fitMethodFermion.GetValue():
            atom = "Li"

        tx, ty, wx, wy = dataFit(atom, tofList, RXList, RYList)
        self.fitTempText.SetValue('(%.1f' %(tx*1E9) + ' , ' + '%.1f )' %(ty*1E9))

        self.fitTrapAxialFreqText.SetValue(str('%.1f' % (wy/(2*np.pi))))
        self.fitTrapRadialFreqText.SetValue(str('%.1f' % (wx/(2*np.pi))))
        # self.fitrho0Text.SetValue(str('%.1f' % (t[2]*1E6)))

    def drawAtomNumber(self, e):
        atomNumberI = []
        atomNumberC = []
        n = len(self.data)
        for i in self.data:
            atomNumberI.append(int(i[2]))
            # atomNumberC.append(int(i[11]))
        atomNumberPlot(n, atomNumberI)

            
	
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