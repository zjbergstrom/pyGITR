import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import io
import libconf
from matplotlib import cm
from mpl_toolkits import mplot3d
from numpy import genfromtxt

def printInfo(ncData):
    print("\nDimensions")
    for dim in ncData.dimensions.values():
       print(dim)

    print("\nVariables")
    for var in ncData.variables.values():
        print(var,"\n")

def getHistoryData(ncData):
    ''' dim_dict: nP, nT
        var_dict: x, y, z, vx, vy, vz '''
    nP = len(ncData.dimensions["nP"])
    nT = len(ncData.dimensions["nT"])

    x = np.array(ncData.variables["x"])
    y = np.array(ncData.variables["y"])
    z = np.array(ncData.variables["z"])
    vx = np.array(ncData.variables["vx"])
    vy = np.array(ncData.variables["vy"])
    vz = np.array(ncData.variables["vz"])

    return {"nP":nP,"nT":nT}, {"x":x,"y":y,"z":z,"vx":vx,"vy":vy,"vz":vz}

def getSurfaceData(ncData):
    ''' dim_dict: nSurfaces, nEnergies, nAngles
        var_dict: grossDeposition[nSurfaces], grossErosion[nSurfaces], 
            aveSpyl[nSurfaces], spylCounts[nSurfaces],
            sumParticlesStrike[nSurfaces], sumWeightStrike[nSurfaces],
            surfEDist[nSurfaces, nEnergies, nAngles],
            surfReflDist[nSurfaces, nEnergies, nAngles],
            surfSputtDist[nSurfaces, nEnergies, nAngles] '''
    nSurfaces = len(ncData.dimensions["nSurfaces"])
    nEnergies = len(ncData.dimensions["nEnergies"])
    nAngles = len(ncData.dimensions["nAngles"])

    grossDeposition = np.array(ncData.variables["grossDeposition"])
    grossErosion = np.array(ncData.variables["grossErosion"])
    aveSpyl = np.array(ncData.variables["aveSpyl"])
    spylCounts = np.array(ncData.variables["spylCounts"])
    sumParticlesStrike = np.array(ncData.variables["sumParticlesStrike"])
    sumWeightStrike = np.array(ncData.variables["sumWeightStrike"])

    surfEDist = np.array(ncData.variables["surfEDist"])
    surfReflDist = np.array(ncData.variables["surfReflDist"])
    surfSputtDist = np.array(ncData.variables["surfSputtDist"])

    return {"nSurfaces":nSurfaces,"nEnergies":nEnergies,"nAngles":nAngles}, \
            {"grossDeposition":grossDeposition,"grossErosion":grossErosion,"aveSpyl":aveSpyl,
                "spylCounts":spylCounts, "sumParticlesStrike":sumParticlesStrike,"sumWeightStrike":sumWeightStrike,
                    "surfEDist":surfEDist,"surfReflDist":surfReflDist,"surfSputtDist":surfSputtDist}

def getSpecData(ncData):
    ''' dim_dict: nBins, nR, nY, nZ, nr_2d, nz_2d
        var_dict: n[nBins, nZ, nY, nR], n_2d[nBins, nz_2d, nr_2d]
              gridR[nR], gridY[nY], gridZ[nZ]
              gridr_2d[nr_2d], gridz_2d[nz_2d] '''
    nBins = len(ncData.dimensions["nBins"])
    nR = len(ncData.dimensions["nR"])
    nY = len(ncData.dimensions["nY"])
    nZ = len(ncData.dimensions["nZ"])
    nr_2d = len(ncData.dimensions["nr_2d"])
    nz_2d = len(ncData.dimensions["nz_2d"])

    n = np.array(ncData.variables["n"])
    n_2d = np.array(ncData.variables["n_2d"])
    gridR = np.array(ncData.variables["gridR"])
    gridY = np.array(ncData.variables["gridY"])
    gridZ = np.array(ncData.variables["gridZ"])
    gridr_2d = np.array(ncData.variables["gridr_2d"])
    gridz_2d = np.array(ncData.variables["gridz_2d"])

    return {"nBins":nBins,"nR":nR,"nY":nY,"nZ":nZ,"nr_2d":nr_2d,"nz_2d":nz_2d}, \
            {"n":n,"n_2d":n_2d,"gridR":gridR,"gridY":gridY,"gridZ":gridZ,"gridr_2d":gridr_2d,"gridz_2d":gridz_2d}

def getParticlesData(ncData):
    ''' dim_dict: nP
        var_dict: x, y, z, vx, vy, vz '''
    nP = len(ncData.dimensions["nP"])

    x = np.array(ncData.variables["x"])
    y = np.array(ncData.variables["y"])
    z = np.array(ncData.variables["z"])
    vx = np.array(ncData.variables["vx"])
    vy = np.array(ncData.variables["vy"])
    vz = np.array(ncData.variables["vz"])

    return {"nP":nP}, {"x":x,"y":y,"z":z,"vx":vx,"vy":vy,"vz":vz}

def getPlasmaData(ncData):
    ''' dim_dict: nR, nZ
        var_dict: r, z, br, bt, bz, psi '''
    nR = len(ncData.dimensions["nR"])
    nZ = len(ncData.dimensions["nZ"])

    r = np.array(ncData.variables["r"])
    z = np.array(ncData.variables["z"])
    br = np.array(ncData.variables["br"])
    bt = np.array(ncData.variables["bt"])
    bz = np.array(ncData.variables["bz"])
    psi = np.array(ncData.variables["psi"])

    return {"nR":nR,"nZ":nZ}, {"r":r,"z":z,"br":br,"by":bt,"bz":bz,"psi":psi}

def getPositionsData(ncData):
    ''' dim_dict: nP
        var_dict: x, y, z, vx, vy, vz 
            transitTime, hitWall, surfaceHit
            weight, charge, hasLeaked,
            distTraveled, time, dt'''
    nP = len(ncData.dimensions["nP"])

    x = np.array(ncData.variables["x"])
    y = np.array(ncData.variables["y"])
    z = np.array(ncData.variables["z"])
    vx = np.array(ncData.variables["vx"])
    vy = np.array(ncData.variables["vy"])
    vz = np.array(ncData.variables["vz"])
    transitTime = np.array(ncData.variables["transitTime"])
    hitWall = np.array(ncData.variables["hitWall"])
    surfaceHit = np.array(ncData.variables["surfaceHit"])
    weight = np.array(ncData.variables["weight"])
    charge = np.array(ncData.variables["charge"])
    hasLeaked = np.array(ncData.variables["hasLeaked"])
    distTraveled = np.array(ncData.variables["distTraveled"])
    time = np.array(ncData.variables["time"])
    dt = np.array(ncData.variables["dt"])

    return {"nP":nP}, {"x":x,"y":y,"z":z,"vx":vx,"vy":vy,"vz":vz, \
            "transitTime":transitTime,"hitWall":hitWall,"surfaceHit":surfaceHit, \
                "weight":weight,"charge":charge,"hasLeaked":hasLeaked, \
                    "distTraveled":distTraveled,"time":time,"dt":dt}


def plotHistory2D(dim_dict,var_dict,skip=10):
    for i in range(0,dim_dict["nP"],skip):
        X = var_dict["x"][i]
        Y = var_dict["y"][i]
        Z = var_dict["z"][i]
        R = np.sqrt(np.multiply(X,X) + np.multiply(Y,Y))
        plt.plot(R,Z,linewidth=1.0)
    plt.axis('scaled')
    vessel_x, vessel_y = getVessel()
    # plt.plot(wall[0,:],wall[1,:],c='k')
    plt.plot(vessel_x,vessel_y,c='k')
    plt.title('2D tracks, skipping every {} particles'.format(skip))
    plt.show()

def plotHistory3D(dim_dict,var_dict,skip=100):
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    for i in range(0,dim_dict["nP"],skip):
        X = var_dict["x"][i]
        Y = var_dict["y"][i]
        Z = var_dict["z"][i]
        # ax.scatter3D(X,Y,Z,s=20)
        ax.plot3D(X,Y,Z,linewidth=1.0)
    plt.title('3D tracks, skipping every {} particles'.format(skip))
    plt.show()

def plotPlasma(dim_dict,var_dict):
    R = var_dict["r"]
    Z = var_dict["z"]
    # Br = var_dict["br"]
    # plt.scatter(R,Z,s=10)
    # plt.show()
    for key in var_dict.keys():
        # plt.scatter(R,Z,s=10)
        plt.scatter(R,Z,s=10,vmax=max(var_dict[key]), vmin=min(var_dict[key]),  c=var_dict[key], label="{}".format(key))
        plt.legend()
        plt.show()

def getVessel():
    wall = genfromtxt('wallcurve_nx_option2.txt')
    np.append(wall,np.reshape(wall[0],(1,2)),axis=0)
    return wall[:,0], wall[:,1]

def findHits(dim_dict,var_dict,skip=1,plotLongTracks=True):
    hit_dict = {}
    for i in range(0,dim_dict["nP"],1):
        X = var_dict["x"][i]
        Y = var_dict["y"][i]
        Z = var_dict["z"][i]
        R = np.sqrt(np.multiply(X,X) + np.multiply(Y,Y))
        for j in range(1,dim_dict["nT"]):
            delR = R[j]-R[j-1]
            if delR==0:
                # print("particle",i,"hit the wall at",j)
                hit_dict[i] = j
                break

    plt.hist(hit_dict.values())
    plt.yscale("log")
    plt.title('Histogram of hit times')
    plt.ylabel('Counts')
    plt.xlabel('nT')
    plt.show()

    if plotLongTracks:
        long_tracks = []
        for key in hit_dict.keys():
            if hit_dict[key] > 50:
                # print(key,hit_dict[key])
                long_tracks.append(int(key))
        
        # fig = plt.figure()
        # ax = plt.axes(projection="3d")
        for i in long_tracks:
            X = var_dict["x"][i]
            Y = var_dict["y"][i]
            Z = var_dict["z"][i]
            R = np.sqrt(np.multiply(X,X) + np.multiply(Y,Y))
            plt.plot(R,Z,linewidth=1.0)
            # ax.plot3D(X,Y,Z)

        plt.axis('scaled')
        vessel_x, vessel_y = getVessel()
        plt.plot(vessel_x,vessel_y,c='k')
        plt.title('2D tracks of long-lived hits ($\Delta$nT > 50)')
        plt.show()

    return hit_dict

def plotSurface(dim_dict,var_dict,variable=None):
    for var in var_dict.keys():
        if variable==None or var in variable:
            if len(var_dict[var].shape)==1:
                print(var, var_dict[var].shape, "\n", var_dict[var],"\n")
                plt.plot(np.arange(1,dim_dict["nSurfaces"]+1), var_dict[var], 'b-d', label="{}".format(var))
                plt.ylabel("{}".format(var))
                plt.xlabel("surface number")
                plt.legend()
                plt.show()
            else:
                for surface in range(0,1): #dim_dict["nSurfaces"]):
                    # print(var_dict["surfEDist"][surface])
                    x = np.arange(0,dim_dict["nEnergies"],1)
                    y = np.arange(0,dim_dict["nAngles"],1)
                    c = plt.pcolormesh(var_dict["surfEDist"][surface], cmap ='Greens', vmin = np.min(var_dict["surfEDist"][surface]), \
                        vmax = np.max(var_dict["surfEDist"][surface]), label="{}, surface #{}".format(var,surface))
                    plt.colorbar(c)
                    plt.legend()
                    plt.show()

def plotSpec(dim_dict,var_dict):

    c = plt.pcolormesh(var_dict["n_2d"][0], cmap ='Blues', vmin = np.min(var_dict["n_2d"][0]), \
        vmax = np.max(var_dict["n_2d"][0]), label="{}".format("n_2d"))
    plt.colorbar(c)
    plt.legend()
    plt.show()

    n3D = np.array(var_dict["n"][0])

    # # Across Z
    for i in range(50,56):
        c = plt.pcolormesh(n3D[i], cmap ='Reds', vmin = np.min(n3D[i]), \
            vmax = np.max(n3D[i]), label="{}".format("n"))
        plt.colorbar(c)
        plt.title("{}".format(i))
        plt.show()

    # Across Y
    for i in range(38,184,2):
        c = plt.pcolormesh(n3D[:,i,:], cmap ='Reds', vmin = np.min(n3D[:,i,:]), \
            vmax = np.max(n3D[:,i,:]), label="{}".format("n"))
        plt.colorbar(c)
        plt.title("{}".format(i))
        plt.show()

    # Across R
    for i in range(38,184,2):
        c = plt.pcolormesh(n3D[:,:,i], cmap ='Reds', vmin = np.min(n3D[:,:,i]), \
            vmax = np.max(n3D[:,:,i]), label="{}".format("n"))
        plt.colorbar(c)
        plt.title("{}".format(i))
        plt.show()