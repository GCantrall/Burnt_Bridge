import numpy as np
import random
import math
import sys
from Simulation import Simulation
from DataFile import DataSet
import matplotlib.pyplot as plt
from matplotlib import colors
import os.path


"""Loads a group of Dataset objects and averages the values"""
def LoadGroup(idMin = 1, idMax = 2, versionMin = 1, versionMax = 1, replicates=1000, s_length=1000000, lp=20, tp=100, tb=4000, path = "", name=""):
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


"""Plots Log of the Mean Squared Distance"""
def PlotLogMSD(Data):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_yscale('log',base=10)
    ax1.set_xscale('log',base=10)
    ax2.set_ylabel("Peptides Remaining")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mean Squared Distance")
    ax1.plot((Data.timescale),(Data.MSD), c = 'g')
    ax1.plot((Data.timescale),(Data.timescale),c='k')
    ax2.plot((Data.timescale),Data.peptides,c = 'r')
    ax2.set_ylim(0,np.max(Data.peptides))

""" Plots multiple LogMSD graphs on the same plot"""
def PlotMultipleLogMSD(DataList):
    fig, ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    #ax2 = ax1.twinx()

    ax1.set_yscale('log',base=10)
    ax2.set_xscale('log',base=10)

    ax1.set_xscale('log',base=10)
    ax2.set_ylabel("Peptides Remaining")
    ax1.set_xlabel("Time")
    ax2.set_xlabel("Time")
    ax1.set_ylabel("Mean Squared Distance")


    pepMax = 0
    maxTime =0
    count  = 0
    countp =1
    diff = .4/len(DataList)
    count  = diff
    countp =.5+diff
    maxTime = 0
    for Data in DataList:
        ax1.plot((Data.timescale), (Data.MSD), label = Data.name)
        ax2.plot((Data.timescale), Data.peptides, label = Data.name)
        count = count+diff
        countp = countp+diff
        if pepMax < np.max(Data.peptides):
            pepMax = np.max(Data.peptides)
        if(maxTime<np.max(Data.timescale)):
            maxTime=np.max(Data.timescale)

    ax1.plot([0,maxTime], [0,maxTime],c='k',linestyle='--', label = "Normal Diffusion")
    ax1.plot(np.arange(0, maxTime), np.arange(0, maxTime)*4, c='grey', linestyle='--', label = "Fast Diffusion")
    ax1.plot(np.arange(0, maxTime), np.arange(0, maxTime) * 3.2, c='lightgrey', linestyle='--', label="Facilitated Motion")
    ax1.legend()
    ax2.legend()
    ax2.set_ylim(0, pepMax)
    #ax1.set_xlim(DataList[0].timescale[-1])
    #ax1.set_ylim(bottom=100)
    #ax2.set_xlim(100, DataList[0].timescale[-1])

"""Plots Total means squared distance vs time"""
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

def PlotAngleFrequency(Data):
    plt.plot(np.arange(len(Data.angles) - 1) / (10) - np.pi, Data.angles[1:] / sum(Data.angles[1:]))

def PlotMultipleAngleFrequency(Datas):
    fig, ax1 = plt.subplots()
    for Data in Datas:
        ax1.plot(np.arange(len(Data.angles) - 1) / (10) - np.pi, Data.angles[1:] / sum(Data.angles[1:]),label = Data.name)
    ax1.set_ylabel("Probability Distribution")
    ax1.legend()
    ax1.set_xlabel("Angle")
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

def PlotKuhn(Data):
    plt.boxplot(Data.Kuhn,0,'')
    print(np.min(Data.Kuhn))
#plotTestFunction()
#plt.show()

#Data1  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=5000000, versionMax=1,tp=500,lp=20, tb=4000,path="Augmented", name="Normal")
#Data2  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=5000000, versionMax=1,tp=1000,lp=20, tb=4000,path="Augmented", name="Half Insertion Rate")
#Data3  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=5000000, versionMax=1,tp=500,lp=40, tb=4000,path="Augmented", name="Douple Move Distance")



#Distince1  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="Distince", name="Normal")
#Distince2  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=1000,lp=20, tb=4000,path="Distince", name="Half Insertion Rate")
#Distince3  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=40, tb=4000,path="Distince", name="Douple Move Distance")

"""
Normal1  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="Normal", name="Normal")
Normal2  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=1000,lp=20, tb=4000,path="Normal", name="Half Insertion Rate")
Normal3  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=40, tb=4000,path="Normal", name="Douple Move Distance")
"""

#Data1 = LoadGroup(idMin=1,idMax=20,replicates= 10000, versionMin=1, s_length=10000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_Unblocked_Updated", name="Unblocked")
#Data2 = LoadGroup(idMin=1,idMax=20,replicates= 1000, versionMin=1, s_length=10000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_Unblocked", name="Unblocked2")
Data1  = LoadGroup(idMin=1,idMax=20,versionMin=1, replicates=10000, s_length=10000, versionMax=1,tp=500,lp=20, tb=4000, path="Analytics_Diffusion_Only", name="Diffusion")
#unblocked  = LoadGroup(idMin=1,idMax=20,versionMin=1,replicates=10000, s_length=10000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_Unblocked_Updated", name="Unblocked")
#directional  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="Distince", name="Directional")
#Data = DataSet(replicates=500,s_length=10000,version=2)
#Small = DataSet(replicates=100,s_length=2000, lp=20,sim="b")
#Small.LoadData()


#PlotAngleFrequency(Small)


PlotMultipleLogMSD([Data1])



plt.show()