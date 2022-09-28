import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import io
import libconf # pip install libconf
import netCDF4
from netCDF4 import Dataset

minR = 1.37875
maxR = 1.59125
minY =  -0.10625
maxY = 0.10625
minZ = -1.25
maxZ = -1.0375

def Plasma(axis,plasmaFile):

    ncFile = netCDF4.Dataset(plasmaFile)
    R = np.array(ncFile.variables["r"])
    Z = np.array(ncFile.variables["z"])
    # R = np.array(ncFile.variables["gridR"])
    # Z = np.array(ncFile.variables["gridZ"])
    Br = np.array(ncFile.variables["br"])
    Te = np.array(ncFile.variables["te"])
    ne = np.array(ncFile.variables["ne"])

    R,Z = np.meshgrid(R,Z)
    R = np.reshape(R,(-1,1))
    Z = np.reshape(Z,(-1,1))
    # c=axis.scatter3D(R, np.array([0])*len(R), Z, vmax=np.max(Te), vmin=np.min(Te),  c=Te, cmap='Greens', label="Te", marker='s')
    c=axis.scatter3D(R, np.array([0])*len(R), Z, vmax=np.max(ne), vmin=np.min(ne),  c=ne, cmap='Blues', label="ne", marker='s')
    plt.colorbar(c)


def plotTracks(filename='../output/history_0.nc',geomFile='../input/gitrGeom.cfg',plasmaFile='../input/profiles.nc',plot=1):
    
    ncFile = netCDF4.Dataset(filename,'r',format="NETCDF4")
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

        for i in range(0,nP,1):
            ax.plot3D(x[i,:],y[i,:],z[i,:],linewidth=1.0)

        r_shift = 1.485
        z_shift = -1.250

        x = np.sin(np.linspace(0,2*np.pi,101))
        y = np.cos(np.linspace(0,2*np.pi,101))
        z = [0.0001]*len(x)
        z = np.array(z)+z_shift

        xdimes = x*0.025+r_shift
        ydimes = y*0.025
        xld = x*0.005+r_shift
        yld = y*0.005
        xsd = (x*0.0005)+r_shift
        ysd = (y*0.0005)+0.01
        ymr = [-0.10625,0.10625,0.10625,-0.10625,-0.10625]
        xmr = np.array([-0.085,-0.085,-0.035,-0.035,-0.085])+r_shift
        zmr = [0.0001+z_shift]*len(xmr)

        ax.plot3D(xdimes,ydimes,z,color='gray')
        # ax.plot3D(xld,yld,z,color='purple')
        # ax.plot3D(xsd,ysd,z,color='m')
        ax.plot3D(xmr,ymr,zmr,color='g')

        ax.set_xlim([minR,maxR])
        ax.set_ylim([minY,maxY])
        ax.set_zlim([minZ,maxZ])

        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')

        Plasma(ax,plasmaFile)

        plt.show()
        

if __name__ == "__main__":
    plotTracks(plot=1)
