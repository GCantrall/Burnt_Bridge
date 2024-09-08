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
def LoadGroup(idMin = 1, idMax = 2, versionMin = 1, versionMax = 16, replicates=1000, s_length=1000000, lp=20, tp=100, tb=4000, path = "", name="", simType = ""):
    Data1 = DataSet(replicates=replicates,s_length=s_length,lp=lp,tp=tp,tb=tb,path=path, id=idMin, version=versionMin, name= name,simType=simType)
    Data1.LoadData()
    for id in range(idMin,idMax+1):
        for version in range(versionMin,versionMax+1):
            if version==versionMin and id == idMin:
                continue
            else:
                Data2 = DataSet(replicates=replicates,s_length=s_length,lp=lp,tp=tp,tb=tb,path=path, id=id, version=version, simType=simType)
                Data2.LoadData()
                Data1.Average(Data2)
    return Data1



"""Loads a group of Dataset objects and averages the values"""
def LoadGroupPath(path = "", name=""):

    files = os.listdir(path)

    Data1 = DataSet()
    t = 0
    for file in files:
        if "Simulation_" in file:
            if t==0:
                Data1.LoadData(path+"/"+ file)
                Data1.name = name
                t = 1
            else:
                Data2 = DataSet()
                Data2.LoadData(path+"/"+ file)
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


def PlotSingleMSD(path):
    files = os.listdir(path)
    fig, ax1 = plt.subplots()
    Data1 = DataSet()
    t = 0
    maxTime = 0

    ax1.set_yscale('log',base=10)
    ax1.set_xscale('log',base=10)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mean Squared Distance")

    for file in files:
        if "Simulation_" in file:
            if t==0:
                Data1.LoadData(path+"/"+ file)
                t = 1
            else:
                Data2 = DataSet()
                Data2.LoadData(path+"/"+ file)
                Data1.Average(Data2)
                t = t+1
                if t >20:
                    continue

                ax1.plot((Data2.timescale[10:]), (Data2.MSD[10:]), c='gray', alpha = .5)
                if(maxTime<np.max(Data2.timescale)):
                    maxTime=np.max(Data2.timescale)
    ax1.plot((Data1.timescale[10:]), (Data1.MSD[10:]))
    ax1.plot([1,maxTime], [1,maxTime],c='k',linestyle='--', label = "Normal Diffusion")




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
        ax1.plot((Data.timescale[10:]), (Data.MSD[10:]), label = Data.name)
        ax2.plot((Data.timescale[10:]), Data.peptides[10:], label = Data.name)
        count = count+diff
        countp = countp+diff
        if pepMax < np.max(Data.peptides):
            pepMax = np.max(Data.peptides)
        if(maxTime<np.max(Data.timescale)):
            maxTime=np.max(Data.timescale)

    ax1.plot([1,maxTime], [1,maxTime],c='k',linestyle='--', label = "Normal Diffusion")
    #ax1.plot(np.arange(0, maxTime), np.arange(0, maxTime)*4.8, c='grey', linestyle='--', label = "Rolling Diffusion")
    ax1.legend()
    ax2.legend()
    ax2.set_ylim(0, pepMax)
    #ax1.legend(title="Rate of Clevage")
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

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def PlotDM(DataSets):
    fig, ax1 = plt.subplots()
    n = 2
    for Data in DataSets:
        ax1.plot((Data.timescale[10+n:]), (moving_average(Data.DM[11:],n)), label = Data.name)
    ax1.set_xscale('log', base=10)
    ax1.legend()


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
    fig, ax1 = plt.subplots()
    ax1.plot(np.arange(len(Data.angles) - 1) / (10) - np.pi, Data.angles[1:] / sum(Data.angles[1:]))
    ax1.set_ylabel("Probability Distribution")
    ax1.set_xlabel("Angle")

def PlotMultipleAngleFrequency(Datas):
    fig, ax1 = plt.subplots()
    for Data in Datas:
        ax1.plot(np.arange(len(Data.angles) - 1) / (10) - np.pi, Data.angles[1:] / sum(Data.angles[1:]),label = Data.name)
    ax1.set_ylabel("Probability Distribution")
    ax1.legend(title = "Strength of Directed Motion")
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

#Data1  = DataSet( replicates=500, version=6, s_length=50, tp=500,lp=20, tb=4000,path="",simType="d", name="No Peptides Cluster")
#Data1.LoadData()
#Data2  = DataSet( replicates=1000, s_length=1000, tp=500,lp=20, tb=4000,path="",simType="d", name="Diffusion Only")
#Data2.LoadData()
#Data3  = DataSet( replicates=1000, s_length=1000, version = 3,  tp=500,lp=20, tb=4000,path="",simType="d", name="No Directed Motion")
#Data3.LoadData()
#Data4  = DataSet( replicates=200, s_length=1000, tp=500,lp=20, tb=4000,path="",simType="d", name="Fast Peptide Diffusion")
#Data4.LoadData()
#Data5  = DataSet( replicates=201, s_length=1000, tp=500,lp=20, tb=4000,path="",simType="d", name="No Directed Motion")
#Data5.LoadData()
#Data4.Average(Data5)


#Data1 = LoadGroupPath("Analytics_AttDM_Final_2_tc_100",name="tc 100")
#Data2 = LoadGroupPath("Analytics_AttDM_Final_2_tc_200", name = "tc 200")
#Data3 = LoadGroupPath("Analytics_AttDM_Final_2_tc_400", name = "tc 400")
#Data1_msd = LoadGroupPath("AttDM_Final_tc_100",name="tc 100")
#Data2_msd = LoadGroupPath("AttDM_Final_tc_200", name = "tc 200")
#Data3_msd = LoadGroupPath("AttDM_Final_tc_400", name = "tc 400")
#Data4 = LoadGroupPath("Analytics_AttDM_tc_400", name = "tc 400")

#Data4 = LoadGroupPath("AttDM_NoDM", name = "No DM")
#Data5 = LoadGroupPath("AttDM_NoPeptide", name = "No Peptide Attraction")

#Data6 = LoadGroupPath("AttDM_Diff", name = "Diff")
"""
Data2 = LoadGroupPath("Analytics_AttDiff_tc_200",name="CR 0.005")
Data1 = LoadGroupPath("Analytics_AttDiff_tc_100", name = "CR 0.01")
Data3 = LoadGroupPath("Analytics_AttDiff_tc_400", name = "CR 0.0025")
"""


#Data1  = LoadGroup(idMin=1,idMax=50,versionMin=1, replicates=50, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="AttDiff_z", name= "Run1",simType="d")
#Normal  = LoadGroup(idMin=1,idMax=60,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="AttDif_tc_200", name="Normal", simType="d")
#Data3  = LoadGroup(idMin=1,idMax=60,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="AttDif_tc_400", name="0.0025", simType="d")


#PlotMultipleLogMSD([Data1_msd,Data2_msd,Data3_msd])
#PlotDM([Data1,Data2,Data3])
#PlotMultipleAngleFrequency([Data1,Data2, Data3])

PlotSingleMSD("Single_AttDM_Final")
PlotSingleMSD("Single_AttDM_Final_tc_100")
PlotSingleMSD("Single_AttDM_Final_tc_100_kp_3")
#Distince1  = LoadGroup(replicates=10000, idMin=1,idMax=20,versionMin=1, s_length=10000, versionMax=3,tp=500,lp=20, tb=4000,path="Analytics_Directional_Updated_Angle", name="Normal")
#Distince2  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=1000,lp=20, tb=4000,path="Distince", name="Half Insertion Rate")
#Distince3  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=40, tb=4000,path="Distince", name="Douple Move Distance")


"""
Normal1  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="Normal", name="Normal")
Normal2  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=1000,lp=20, tb=4000,path="Normal", name="Half Insertion Rate")
Normal3  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=40, tb=4000,path="Normal", name="Douple Move Distance")
"""
#Data = LoadGroup(idMin=1,idMax=20,versionMin=1,replicates=20000, s_length=10000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_Repulsion_2", name="Repulsive")
#Data1 = LoadGroup(idMin=1,idMax=20,replicates= 10000, versionMin=1, s_length=10000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_Unblocked_Updated", name="Unblocked")
#Data2 = LoadGroup(idMin=1,idMax=20,replicates= 1000, versionMin=1, s_length=10000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_Unblocked", name="Unblocked2")
#normal_2  = LoadGroup(idMin=1,idMax=40,versionMin=1, replicates=10000, s_length=30000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_AttDM", name="0.005",simType="a")
#fast_2  = LoadGroup(idMin=1,idMax=40,versionMin=2, replicates=10000, s_length=30000, versionMax=2,tp=500,lp=20, tb=4000,path="Analytics_AttDM", name="0.0025",simType="a")
#slow_2  = LoadGroup(idMin=1,idMax=40,versionMin=1, replicates=10000, s_length=30000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_AttDM_tc_100", name="0.01",simType="a")


#normal  = LoadGroup(idMin=1,idMax=40,versionMin=1, replicates=10000, s_length=30000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_AttDM_kp_1", name="1",simType="a")
#slow  = LoadGroup(idMin=1,idMax=40,versionMin=1, replicates=10000, s_length=30000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_AttDM_kpb_5", name="0.5",simType="d")
#fast  = LoadGroup(idMin=1,idMax=40,versionMin=1, replicates=10000, s_length=30000, versionMax=1,tp=500,lp=20, tb=4000,path="Analytics_AttDM_kpb_1", name="0.1",simType="d")


#Fast  = LoadGroup(idMin=1,idMax=60,versionMin=2,replicates=10000, s_length=10000, versionMax=2,tp=500,lp=20, tb=4000,path="Analytics_Diffuse_Peptide", name="Fast Cutting",simType="b")

#directional  = LoadGroup(idMin=1,idMax=20,versionMin=1, s_length=1000000, versionMax=1,tp=500,lp=20, tb=4000,path="Distince", name="Directional")
#Data = DataSet(replicates=500,s_length=10000,version=2)
#Data.LoadData()
#PlotMultipleLogMSD([Data,Data1])
#PlotKuhn(Data)
#PlotAngleFrequency(Data2)
#PlotMultipleAngleFrequency([normal,slow,fast])
#PlotMultipleAngleFrequency([fast_2,normal_2,slow_2])
#PlotAngleFrequency(Data)
#PlotAngleFrequency(Distince1)
#DataSet(tp=500,lp=100, tb=1000,id=1,version=2,path="DataFolder")

#Data2  = DataSet(replicates=100,s_length=1000000, tp=500)

#PlotMultipleLogMSD([Data])
#

#PlotMultipleLogMSD([Data1,Data2, Data3])

#PlotMultipleLogMSD([Distince1,Distince2, Distince3])

#PlotMultipleLogMSD([Normal1,Normal2, Normal3])

#Data2.LoadData()

#Data1.Average(Data2)


#PlotLogMSD(Data1)
#PlotRunningRMSD(Data1)
#PlotRunningLogMSD(Data1)
#PlotTotalMSD(Data1)

plt.show()