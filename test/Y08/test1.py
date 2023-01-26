import matplotlib.pyplot as plt
import numpy as np

nkpt = 151
nmesh = 501
kmesh = np.zeros(nkpt, dtype = np.int64)
rmesh = np.zeros(nmesh, dtype = np.float64)
chi = np.zeros((nmesh,nkpt), dtype = np.float64)

g = open("chir.data", "r")
for k in range(nkpt):
    for m in range(nmesh):
        line = g.readline().split()
        kmesh[k] = int(line[0])
        rmesh[m] = float(line[1])
        chi[m,k] = -float(line[3]) / np.pi
g.close()

extent = np.min(kmesh), np.max(kmesh), np.min(rmesh), np.max(rmesh)
im1 = plt.imshow(chi,
                 cmap=plt.cm.hot,   #viridis,
                 interpolation='bicubic',
                 extent = extent,
                 origin = "lower",
                 aspect = "auto")
plt.show()
