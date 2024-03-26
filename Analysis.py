import numpy as np
import random
import math
import sys
from Simulation import Simulation
from DataFile import DataSet
import matplotlib.pyplot as plt
from matplotlib import colors
import os.path



def LoadGroup(idMin = 1, idMax = 2, versionMin = 1, versionMax = 16, replicates=1000, s_length=1000000, lp=20, tp=100, tb=4000, path = "", name=""):
    Data1 = DataSet(replicates=replicates,s_length=s_length,lp=lp,tp=tp,tb=tb,path=path, id=idMin, version=versionMin, name= name)
    Data1.LoadData()
    for id in range(idMin,idMax+1):
        for version in range(versionMin,versionMax+1):
            if version==versionMin and id == idMin:
                continue
            else:
                Data2 = DataSet(replicates=replicates,s_length=s_length,lp=lp,tp=tp,tb=tb,path=path, id=id, version=version)
                Data2.LoadData()
                Data1.Average(Data2)
    return Data1


def PlotLogMSD(Data):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_yscale('log',base=10)
    ax1.set_xscale('log',base=10)
    ax2.set_ylabel("Peptides Remaining")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mean Squared Distance")
    ax1.plot((Data.timescale),(Data.MSD))
    ax1.plot((Data.timescale),(Data.timescale))
    ax2.plot((Data.timescale),Data.peptides,c = 'r')
    ax2.set_ylim(0,np.max(Data.peptides))
def PlotMultipleLogMSD(DataList):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_yscale('log',base=10)
    ax1.set_xscale('log',base=10)
    ax2.set_ylabel("Peptides Remaining")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mean Squared Distance")

    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 0.5, 1.0),
                     (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (0.5, 1.0, 1.0),
                       (1.0, 0.0, 1.0)),

             'blue': ((0.0, 0.0, 0.0),
                      (0.5, 0.0, 1.0),
                      (1.0, 0.0, 1.0))}
    cmap = colors.LinearSegmentedColormap('custom', cdict)
    pepMax = 0
    maxTime =0
    count  = 0
    countp =1
    diff = .4/len(DataList)
    count  = diff
    countp =.5+diff
    maxTime = 0
    for Data in DataList:
        ax1.plot((Data.timescale), (Data.MSD), label = Data.name,c=cmap(count))
        ax2.plot((Data.timescale), Data.peptides, label = Data.name+" peptides",c=cmap(countp))
        count = count+diff
        countp = countp+diff
        if pepMax < np.max(Data.peptides):
            pepMax = np.max(Data.peptides)
        if(maxTime<np.max(Data.timescale)):
            maxTime=np.max(Data.timescale)

    ax1.plot([0,maxTime], [0,maxTime],c='k')

    ax1.legend()
    ax2.legend(loc='center left')
    ax2.set_ylim(0, pepMax)
def PlotTotalMSD(Data):

    fig2, ax1_b = plt.subplots()
    ax2_b = ax1_b.twinx()

    ax1_b.plot(Data.timescale,(Data.MSD))
    ax2_b.plot(Data.timescale,Data.peptides,c = 'r')
    ax2_b.set_ylim(0, np.max(Data.peptides))

    ax2_b.set_ylabel("Peptides Remaining")
    ax1_b.set_xlabel("Time")
    ax1_b.set_ylabel("Mean Squared Distance")


def PlotRunningRMSD(Data):
    fig3, ax1_c = plt.subplots()
    ax1_c.plot(Data.timescale[-len(Data.RMSDw):],Data.RMSDw)
    ax2_c = ax1_c.twinx()
    ax2_c.plot(Data.timescale,Data.peptides,c = 'r')
    ax2_c.set_ylim(0, np.max(Data.peptides))

    ax2_c.set_ylabel("Peptides Remaining")
    ax1_c.set_xlabel("Time")
    ax1_c.set_ylabel("Running Root Mean Squared Distance ")


def plotTestFunction():
    fig, ax = plt.subplots(figsize=(6, 6))


    cdict = {'red': ((0.0, 0.22, 0.0),
                     (0.5, 1.0, 1.0),
                     (1.0, 0.89, 1.0)),

             'green': ((0.0, 0.49, 0.0),
                       (0.5, 1.0, 1.0),
                       (1.0, 0.12, 1.0)),

             'blue': ((0.0, 0.72, 0.0),
                      (0.5, 0.0, 0.0),
                      (1.0, 0.11, 1.0))}
    cmap = colors.LinearSegmentedColormap('custom', cdict)

    for i in np.linspace(0, 1):
        # Plot 50 lines, from y = 0 to y = 1, taking a corresponding value from the cmap
        ax.plot([-1, 1], [i, i], c=cmap(i))

#plotTestFunction()
#plt.show()

Data1  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="Directional", name="Normal")

Data3  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=1000,lp=20, tb=4000,path="Directional_2", name="Half Insertion Rate")

Data4  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=40, tb=4000,path="Directional_2", name="Douple Move Distance")



#DataSet(tp=500,lp=100, tb=1000,id=1,version=2,path="DataFolder")

#Data2  = DataSet(replicates=100,s_length=1000000, tp=500)

PlotMultipleLogMSD([Data1,Data4, Data3])



#Data2.LoadData()

#Data1.Average(Data2)


#PlotLogMSD(Data1)
#PlotRunningRMSD(Data1)
#PlotRunningLogMSD(Data1)
#PlotTotalMSD(Data1)

plt.show()