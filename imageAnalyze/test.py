from polylog import *
from imgFunc import *
import numpy as np
import matplotlib.pyplot as plt


y = np.array(range(-100, 100))/10.

for x in range(20):
	z = np.exp(x/10.)
	n = polylog5half(np.exp(z-y**2*f(np.exp(z))))/polylog5half(np.exp(z))
	plt.plot(y, n)
	plt.show()
