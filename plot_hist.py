import numpy as np
import matplotlib.pyplot as plt

infile = "randomnumbers.txt"
figurepng = "histogram.png"
figurepdf = "histogram.pdf"

with open(infile) as f:
    lines = f.readlines()
    ihstr = [ float( line.split()[1] ) for line in lines]

numelmnts = len( ihstr )
iharr = np.asarray(ihstr)
bins  = np.arange( 0.0, 1.0, 0.05 )

fig=plt.figure(figsize=(6,5))
n, bins, patches = plt.hist(iharr, bins = bins, facecolor='red', align='mid')
plt.show()

