# Math support functions
import numpy as np
import math3d as m3d
import matplotlib.pyplot as plt

def RotateCoordinates(x:np.ndarray, y:np.ndarray, z:np.ndarray, AxisVector:np.ndarray, Angle: float, Degree=True):
    v=np.vstack((x,y,z)).transpose()
    return RotateVector(v, AxisVector, Angle, Degree)


def RotateVector(v:list or np.ndarray, AxisVector:np.ndarray, Angle: float, Degree=True):
    if type(v) == list:
        v = np.array(v)

    if type(AxisVector) == list:
        AxisVector = np.array(AxisVector)
    if Degree:
        Angle = np.radians(Angle)
    r = m3d.Orientation.new_axis_angle(AxisVector, Angle)
    # print('r=',r)
    # print('v=',v.shape)
    vp=np.copy(v)
    if len(v.shape)>1:
        for i in range(0,v.shape[0]):
            vp[i,:] = r*v[i,:]
    else:
        vp=r*v
    return vp

def Integrale(f:np.ndarray, x:np.ndarray, Array=True, Normalized = False):
    Integral = np.cumsum((f[1:]+f[:-1])/2.0*(x[1:]-x[:-1]))
    Integral = np.insert(Integral,0,0)
    if Normalized and Integral[-1] != 0:
        Integral = Integral/Integral[-1]

    if Array:
        return Integral
    else:
        return Integral[-1]


# Methods for computations of geometrical properties

def Norm(Vector):
        return np.sqrt(np.sum(Vector**2, 1))


def Cross(V1, V2):
        CrossVec = np.copy(V1)
        CrossVec[:, 0] = V1[:, 1]*V2[:, 2]-V1[:, 2]*V2[:, 1]
        CrossVec[:, 1] = V1[:, 2]*V2[:, 0]-V1[:, 0]*V2[:, 2]
        CrossVec[:, 2] = V1[:, 0]*V2[:, 1]-V1[:, 1]*V2[:, 0]
        return CrossVec


def Dot(V1, V2):
        return np.sum(V1*V2, 1)

def PlotVector(V, L=1.2):

    if type(V) == list:
        V = np.array(V)

    assert type(V) == np.ndarray, "V must be a numpy array"
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    if len(V.shape)<2:
        V = V[:,None].transpose()
    ax.quiver(0, 0, 0, V[:,0], V[:,1], V[:,2])
    ax.quiver(0, 0, 0, 1,0,0,pivot='tail',color='k')
    ax.quiver(0, 0, 0, 0,1,0,pivot='tail',color='k')
    ax.quiver(0, 0, 0, 0,0,1,pivot='tail',color='k')
    ax.set_xlim3d(-L, L)
    ax.set_ylim3d(-L, L)
    ax.set_zlim3d(-L, L)
    ax.text(0, 0, 1, "z", color='k')
    ax.text(1, 0, 0, "x", color='k')
    ax.text(0, 1, 0, "y", color='k')