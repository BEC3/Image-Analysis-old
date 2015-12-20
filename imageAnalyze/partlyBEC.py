# partly condensate fit
def partlyCondensateFit(data):

    ### get slides
    xSlice = np.sum(data,0)    
    ySlice = np.sum(data,1)


    ## define 1D distribution function
    def oneDPartlyCondensate(x, centerX, AmpP, a, sDevX, AmpG, offset):
        return np.maximum(AmpG, 0)*np.exp(-0.5*((x-centerX)/sDevX)**2) + np.maximum(AmpP, 0) *  np.maximum( (1 - a * (x - centerX)**2), 0)  ** 2.5+ offset
       

        
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
    AmpPX = AmpGX
    AmpPY = AmpGY

            
    ### Do 1D fitting
    xVals, yCovar = curve_fit(oneDPartlyCondensate,range(len(xSlice)),xSlice,p0=(x0, AmpPX/2, aX, sigmaX, AmpGX/2, xOff))
    x0 = xVals[0]
    widthX = sqrt(1/xVals[2])
    sigmaX = xVals[3]
    
    yVals, yCovar = curve_fit(oneDPartlyCondensate,range(len(ySlice)),ySlice,p0=(y0, AmpPY/2, aY, sigmaY, AmpGY/2, yOff))
    y0 = yVals[0]
    widthY = sqrt(1/yVals[2])
    sigmaY = yVals[3]

    ###### Transfer output to 2D
    AmpP = sqrt(np.maximum(xVals[1], 0) * np.maximum(yVals[1], 0) * 9./16. * sqrt(xVals[2] * yVals[2]) )
    AmpG = 0.5 * ( np.maximum(xVals[4], 0)/(sqrt(2.0*np.pi)*sigmaY) + np.maximum(yVals[4], 0)/(sqrt(2.0*np.pi)*sigmaX) ) 
    offset = 0.5 * (xVals[5]/(np.shape(data)[1]) + yVals[5]/(np.shape(data)[0]))


    return [[x0,y0], [widthX, widthY], [sigmaX,sigmaY], AmpP, AmpG, offset]
