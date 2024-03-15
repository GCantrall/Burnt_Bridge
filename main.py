import numpy as np
import random
import math
import sys
from Simulation import Simulation
import matplotlib.pyplot as plt
import os.path




def GetRMSDW(x,y,inc):
    RMSD = []
    for i in range(inc,len(x)):
        RMSD.append(np.sqrt(pow(x[i]-x[i-inc],2)+pow(y[i]-y[i-inc],2)))

    return RMSD

def GetMSD(x,y):
    MSD = []
    for i in range(len(x)):
        MSD.append(pow(x[i],2)+pow(y[i],2))
    return MSD

def GetParams():
    s_length  = 1000000
    replicates = 1000
    tp = 500
    lp = 20
    id = -1
    tb = 4000
    arguments = sys.argv
    w = 10000
    plot = -1
    for i in range(int((len(arguments)-1)/2)):
        if arguments[(2*i+1)] == "-s":
            s_length = int(arguments[2*i+2])
        elif arguments[(2*i+1)] == "-r":
            replicates = int(arguments[2*i+2])
        elif arguments[(2*i+1)] == "-tp":
            tp = int(arguments[2*i+2])
        elif arguments[(2 * i + 1)] == "-lp":
            lp = int(arguments[2 * i + 2])
        elif arguments[(2*i+1)] == "-id":
            id = int(arguments[2*i+2])
        elif arguments[(2*i+1)] == "-tb":
            tb = int(arguments[2*i+2])
        elif arguments[(2*i+1)] =="-w":
            w = int(arguments[2*i+2])
        elif arguments[(2*i+1)] =="-p":
            plot = int(arguments[2*i+2])

    return [s_length, replicates, tp, lp, tb, id,w,plot]


if __name__ == '__main__':
    s_length, replicates, tp, lp, tb, id, RMSDw_l,plot = GetParams()

    S = Simulation(tp,lp,tb)
    times = []

    p = np.arange(0,1,.1)
    l = np.arange(1, s_length)
    timescale = np.concatenate((p,l))
    peptide_remaining = np.zeros(len(timescale))
    x = np.zeros(len(timescale))
    y = np.zeros(len(timescale))
    RMSDw = np.zeros(len(timescale)- RMSDw_l)
    MSD = np.zeros(len(timescale))
    for i in range(replicates):
        print(i)

        time, x_tracker, y_tracker = S.RunSimulation(s_length)
        if(plot!=-1):
            if(plot==1 or plot==3):
                S.PlotPath(x_tracker,y_tracker)
            if(plot==2 or plot==3):
                S.particle.PlotParticle()
            plt.show()
        peptide_unif = []
        x_unif = []
        y_unif = []

        k = 0

        for j in range(len(timescale)):
            while True:
                if time[k+1] < timescale[j]:
                    k = k+1
                else:
                    break


            peptide_unif.append(S.peptide_remain[k])
            x_unif.append(x_tracker[k])
            y_unif.append(y_tracker[k])
        MSD = (i * MSD + np.array(GetMSD(x_unif, y_unif))) / (i + 1)
        RMSDw = (i * RMSDw + np.array(GetRMSDW(x_unif, y_unif, RMSDw_l))) / (i + 1)
        peptide_remaining = (i*peptide_remaining+np.array(peptide_unif))/(i+1)

    filename  = "Simulation_"+str(replicates)+"r_"+str(s_length)+"s_"+str(lp)+"lp_"+str(tp)+"tp_"+str(tb)+"tb"
    if id != -1:
        filename = filename+"_"+str(id)


    if os.path.isfile(filename+".npz"):
        k=2
        while True:
            if os.path.isfile(filename+"("+str(k)+")"+".npz"):
                k = k+1
            else:
                filename = filename+"("+str(k)+")"
                break




    np.savez(filename, peptides = peptide_remaining, MSD =MSD, RMSDw = RMSDw, timescale = timescale)











# See PyCharm help at https://www.jetbrains.com/help/pycharm/
