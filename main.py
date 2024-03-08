import numpy as np
import random
import math
import sys
from Simulation import Simulation
import matplotlib.pyplot as plt
import os.path




def GetMSDR(x,y,inc):
    MSD = []
    for i in range(inc,len(x)):
        MSD.append(np.sqrt(pow(x[i]-x[i-inc],2)+pow(y[i]-y[i-inc],2)))

    return MSD

def GetParams():
    s_length  = 1000000
    replicates = 1000
    tp = 500
    lp = 20
    id = -1
    tb = 1000
    arguments = sys.argv
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

    return [s_length, replicates, tp, lp, tb, id]


if __name__ == '__main__':
    s_length, replicates, tp, lp, tb, id = GetParams()

    S = Simulation(tp,lp,tb)
    times = []

    p = np.arange(0,1,.1)
    l = np.arange(1, s_length)
    timescale = np.concatenate((p,l))
    peptide_remaining = np.zeros(len(timescale))
    x = np.zeros(len(timescale))
    y = np.zeros(len(timescale))

    for i in range(replicates):
        print(i)

        time, x_tracker, y_tracker = S.RunSimulation(s_length)

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




    np.savez(filename, peptides = peptide_remaining, x =x_unif, y=y_unif, timescale = timescale)


    #file = np.load("InfoFile.npz",allow_pickle=True)










# See PyCharm help at https://www.jetbrains.com/help/pycharm/
