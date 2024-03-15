import numpy as np
import random
import math
import sys
from Simulation import Simulation
from DataFile import DataSet
import matplotlib.pyplot as plt
import os.path



def LoadGroup(idMin = 1, idMax = 2, versionMin = 1, versionMax = 16, replicates=1000, s_length=1000000, lp=20, tp=100, tb=2000, path = ""):
    Data1 = DataSet(replicates=replicates,s_length=s_length,lp=lp,tp=tp,tb=tb,path=path, id=idMin, version=versionMin)
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
    ax1.plot((Data.timescale),(Data.MSD))
    ax1.plot((Data.timescale),(Data.timescale))
    ax2.plot((Data.timescale),Data.peptides,c = 'r')
    ax2.set_ylim(0,600)

def PlotTotalMSD(Data):

    fig2, ax1_b = plt.subplots()
    ax2_b = ax1_b.twinx()

    ax1_b.plot(Data.timescale,(Data.MSD))
    ax2_b.plot(Data.timescale,Data.peptides,c = 'r')
    ax2_b.set_ylim(0, 600)


def PlotRunningMSD(Data):
    fig3, ax1_c = plt.subplots()
    ax1_c.plot(Data.timescale[-len(Data.RMSDw):],Data.RMSDw)
    ax2_c = ax1_c.twinx()
    ax2_c.plot(Data.timescale,Data.peptides,c = 'r')
    ax2_c.set_ylim(0, 600)


Data1  = LoadGroup(idMin=1,idMax=20,versionMin=1, versionMax=1,tp=500,lp=20, tb=4000,path="DataFolder3")
#DataSet(tp=500,lp=100, tb=1000,id=1,version=2,path="DataFolder")

#Data2  = DataSet(replicates=100,s_length=1000000, tp=500)





#Data2.LoadData()

#Data1.Average(Data2)


PlotLogMSD(Data1)
PlotRunningMSD(Data1)
PlotTotalMSD(Data1)

plt.show()