import numpy as np
import numpy
import scipy.signal as signal
import matplotlib.pyplot as plt
from scipy.ndimage import laplace


stencil=numpy.array([[0,1,0],[1,-4,1],[0,1,0]])

fig=plt.figure(figsize=(10,10),dpi=100)

for ii,jj in enumerate([10,100,1000,10000]):

    x=numpy.linspace(-5,5,jj)
    xx,yy=numpy.meshgrid(x,x)
    step=x[1]-x[0]

    image=numpy.exp(-xx**2-yy**2)

    lap1=laplace(image)/step**2
    #lap2=filters.convolve(image,stencil,mode='wrap')/step**2
    #lap3=signal.convolve2d(image,stencil,mode='same')/step**2
    lap4=4*image*(xx**2+yy**2)-4*image

    ax=fig.add_subplot(2,2,ii+1)
    img=ax.imshow(lap1-lap4, cmap='plasma')
    ax.set_title('stencil - analytical (dx=%.4f)' %step)
    plt.colorbar(img)

fig.tight_layout()
plt.show()
