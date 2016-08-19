
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.io import FortranFile

Nx=300
Ny=300
nend=3000
iscan=20
dt=1.

# tikei 
f='./tikei_hokan.dat'
f=FortranFile(f,'r')
data=f.read_reals(dtype=float)
lat=data[0:Ny-1]
lon=data[Ny:Ny+Nx-1]
h=data[Ny+Nx:].reshape((Nx,Ny))
f.close()

imgs=[]

fig=plt.figure()
ax =fig.add_subplot(1,1,1)
ax.set_xlabel("lon")
ax.set_ylabel("lat")

print lon.min(0),lon.max(0),lat.max(0),lat.min(0)

for num in range(0,nend,iscan):
    fnum=str(num).zfill(5)
    f="./data/ts_saka_"+fnum+".dat"
    f=FortranFile(f,'r')

    data=f.read_reals(dtype=float)
    
    f.close()

    p=data[1:].reshape((Nx,Ny))

    image =[plt.imshow(p,vmin=-0.5,vmax=0.5,extent=[lon.min(0),lon.max(0),lat.max(0),lat.min(0)]), plt.annotate("t = " + fnum+ "(s)" , xy=(Nx-100,Ny-50))]
    
    imgs.append(image)
    
    if (num == 0):
        plt.colorbar()


ani = animation.ArtistAnimation(fig, imgs , interval=50, blit=True ,repeat_delay=1000)

Writer=animation.writers['ffmpeg']
writer=Writer(fps=30, metadata=dict(artist='Me'),bitrate=1800)

ani.save('./figure/ts.mp4',writer=writer)

plt.show()
