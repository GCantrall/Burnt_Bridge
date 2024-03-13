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

def PlotRMSD(Data):
    fig, ax1 = plt.subplots()
    fig2, ax1_b = plt.subplots()
    fig3, ax1_c = plt.subplots()

    ax2 = ax1.twinx()
    ax2_c = ax1_c.twinx()

    ax1_c.plot(timescale[MSDLen:],MSD)
    ax2_c.plot(timescale,peptide_remaining,c = 'r')

    ax1.set_yscale('log',base=10)
    ax1.set_xscale('log',base=10)
    ax1.plot((timescale),(distances))
    ax1.plot((timescale),(timescale))
    ax2.plot((timescale),peptide_remaining,c = 'r')
    ax2.set_ylim(0,600)

    ax1_b.plot(timescale,(distances))
    ax2_b = ax1_b.twinx()
    ax2_b.plot(timescale,peptide_remaining,c = 'r')






Data1  = LoadGroup(idMin=1,idMax=9,versionMin=1, versionMax=7,tp=500,lp=100, tb=1000,path="DataFolder")
#DataSet(tp=500,lp=100, tb=1000,id=1,version=2,path="DataFolder")

#Data2  = DataSet(replicates=100,s_length=1000000, tp=500)




Data1.LoadData()
#Data2.LoadData()

#Data1.Average(Data2)


PlotLogMSD(Data1)

plt.show()