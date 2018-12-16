import sys
import numpy as np
import warnings
from matplotlib import pyplot as plt
import imageio
warnings.filterwarnings("ignore")

for i in range(8):
    img = np.loadtxt("data/localData"+str(i+1)+".txt")
    imageio.imwrite("data/localData"+str(i+1)+".png", img)
    plt.imshow(img,cmap='gray',vmin=-1,vmax=1)
    plt.title(str(i+1))
    plt.show()