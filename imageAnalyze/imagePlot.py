from pylab import *

import numpy as np

def atomImagePlot(data, caption):
    l = len(data)
    fig = figure(figsize=(15, 8))
    
    for r in range(len(data)):
        positionString = '1' + str(l) + str(r+1)
        subplot(positionString)
        imshow(data[r].real, cmap='jet', aspect='auto', vmin=-1, vmax=5)
        grid(True, color='white')
        title(caption[r])

    #fig.set_cmap('winter')
    #fig.colorbar()
    #colorbar(orientation='horizontal')
    show()

def atomNumberPlot(n, data1, data2):
    x = range(n)
    p1, = plot(x, data1, 'ro')
    p2, = plot(x, data2, 'bo')
    legend([p1, p2], ["By integration", "By Chemical Potential"])
    
    maxNumber = max(max(data1), max(data2))
    axis([-1, n, 0, maxNumber*1.4])
    grid(True)
    legend()
    title('Atom Number')
    show()
