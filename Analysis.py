import numpy as np
import random
import math
import sys
from Simulation import Simulation
from DataFile import DataSet
import matplotlib.pyplot as plt
import os.path




def loadData(replicates=1000, s_length=1000000, lp = 20, tp = 100, tb = 1000, id=-1, version = -1):
    filename = "Simulation_" + str(replicates) + "r_" + str(s_length) + "s_" + str(lp) + "lp_" + str(tp) + "tp_" + str(tb) + "tb"
    if(id!=-1):
        filename = filename+"_"+id
    if version != -1:
        filename = filename +"("+str(version)+")"
    filename = filename+".npz"


    file = np.load(filename, allow_pickle=True)
    return file['peptides'] , file['MSD'], file['RMSDw'], file['timescale']





Data1  = DataSet(replicates=100,s_length=1000000, tp=500,version=2)
Data2  = DataSet(replicates=100,s_length=1000000, tp=500)




Data1.LoadData()
Data2.LoadData()

Data1.Average(Data2)

peptides, MSD, RMSDw, timescale = loadData(replicates=100,s_length=1000000, tp=500,version=2)

MSDLen = 1

def PlotLogMSD():

    RMSDs = np.power(np.array(RMSDs), 2)

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
