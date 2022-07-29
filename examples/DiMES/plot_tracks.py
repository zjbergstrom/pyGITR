import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import io
import libconf # pip install libconf

def plotTracks(filename='output/history.nc',geomFile='input/gitrGeom.cfg',plot=1):
    
    ncFile = netCDF4.Dataset(filename)
    nT = len(ncFile.dimensions['nT']) # number of timesteps
    nP = len(ncFile.dimensions['nP']) # number of particles
    
    x = np.reshape(np.array(ncFile.variables['x']),(nP,nT))
    y = np.reshape(np.array(ncFile.variables['y']),(nP,nT))
    z = np.reshape(np.array(ncFile.variables['z']),(nP,nT))
    vx = np.reshape(np.array(ncFile.variables['vx']),(nP,nT))
    vy = np.reshape(np.array(ncFile.variables['vy']),(nP,nT))
    vz = np.reshape(np.array(ncFile.variables['vz']),(nP,nT))

    r = np.sqrt(np.multiply(x,x) + np.multiply(y,y))

    with io.open(geomFile) as f:
        config = libconf.load(f)
    
    x1_geom = config['geom']['x1']
    x2_geom = config['geom']['x2']
    x3_geom = config['geom']['x3']
    y1_geom = config['geom']['y1']
    y2_geom = config['geom']['y2']
    y3_geom = config['geom']['y3']
    z1_geom = config['geom']['z1']
    z2_geom = config['geom']['z2']
    z3_geom = config['geom']['z3']



    if plot:
        plt.close()
        
        fig = plt.figure()
        ax = plt.axes(projection="3d")

        #ax.scatter3D(np.append(x1_geom,x1_geom[0]),np.append(y1_geom,y1_geom[0]),np.append(z1_geom,z1_geom[0]),color='g',s=0.1)
        #ax.plot3D(np.append(x1_geom,x1_geom[0]),np.append(y1_geom,y1_geom[0]),np.append(z1_geom,z1_geom[0]),color='k',linewidth=0.1)
        #ax.plot_wireframe(np.append(x1_geom,x1_geom[0]),np.append(y1_geom,y1_geom[0]),np.reshape(np.append(z1_geom,z1_geom[0]),(-1,1)), color='g', linewidth=0.1)
        #ax.plot_surface(np.append(x1_geom,x1_geom[0]),np.append(y1_geom,y1_geom[0]),np.reshape(np.append(z1_geom,z1_geom[0]),(-1,1)), color='g', linewidth=0.1)

        for i in range(0,nP,1):
            ax.plot3D(x[i,:],y[i,:],z[i,:],linewidth=0.3)

        x = np.sin(np.linspace(0,2*np.pi,101))
        y = np.cos(np.linspace(0,2*np.pi,101))
        z = [0.0001]*len(x)

        xdimes = x*0.025
        ydimes = y*0.025
        xld = x*0.005
        yld = y*0.005
        xsd = (x*0.0005)-0.01
        ysd = (y*0.0005)

        ax.plot3D(xdimes,ydimes,z,color='gray')
        ax.plot3D(xld,yld,z,color='purple')
        ax.plot3D(xsd,ysd,z,color='m')

        ax.set_xlim([-0.02,0.02])
        ax.set_ylim([-0.01,0.03])
        ax.set_zlim([-0.0,0.04])

        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')
        plt.show()
        

if __name__ == "__main__":
    plotTracks(filename='output/history.nc',geomFile='input/gitrGeom.cfg',plot=1)
