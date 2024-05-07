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
    analytic = -1
    ptype = "n"
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
        elif arguments[(2*i+1)]== "-type":
            ptype = arguments[2*i+2]
        elif arguments[(2*i+1)]== "-a":
            analytic = arguments[2*i+2]


    return [s_length, replicates, tp, lp, tb, id,w,plot, ptype, analytic]


if __name__ == '__main__':
    s_length, replicates, tp, lp, tb, id, RMSDw_l,plot, ptype, analytic = GetParams()

    S = Simulation(tp,lp,tb,ptype, analytic)
    times = []

    p = np.arange(0,1,.1)
    l = np.arange(1, s_length)
    timescale = np.concatenate((p,l))
    peptide_remaining = np.zeros(len(timescale))
    x = np.zeros(len(timescale))
    y = np.zeros(len(timescale))
    RMSDw = np.zeros(len(timescale)- RMSDw_l)
    MSD = np.zeros(len(timescale))
    KuhnTotal = []

    angleDistTotal = np.zeros(round(2*np.pi*10))
    for i in range(replicates):
        print(i)


        if analytic!=-1:
            time, x_tracker, y_tracker, [Kuhn, angles] = S.RunSimulation(s_length)

            KuhnTotal.append(Kuhn)
            #angleDist = np.zeros(round(2 * np.pi * 100) + 1)
            #for angle in angles:
            #    angleDist[round((angle+np.pi)*100)] +=1
            #angleDist[0] +=angleDistTotal[-1]
            #angleDist = angleDist[:-1]
            #angleDistTotal = (i * angleDistTotal + angleDist)/(i + 1)


        else:
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
        if(plot!=-1):
            if(plot==1 or plot==3 or plot==5):
                S.PlotPath(x_tracker,y_tracker)
            if(plot==2 or plot==3):
                S.particle.PlotParticle()
            if(plot==4 or plot == 5):
                S.PlotPathRange(timescale,x_unif,y_unif)
            plt.show()

        if analytic != -1:
            a = [1,0]
            prev_coords = [0, 0]
            prev_dist = 1
            angles = []
            for j in range(1,len(timescale)):
                if timescale[j]%100==0:
                    b = [x_unif[j] - prev_coords[0], y_unif[j] - prev_coords[1]]
                    ncross = np.cross(a, b)
                    dist = np.sqrt(np.power(b[0],2)+np.power(b[1],2))

                    diff = np.dot(a, b) / (prev_dist*dist)
                    if diff > 1:
                        diff = .999999
                    elif (diff < -1):
                        diff = -.9999999

                    if (ncross == 0 and timescale[j]!=100):
                        angles.append(0)
                    elif(timescale[j]!= 100):
                        angles.append(np.arccos(diff) * ncross / np.sqrt(ncross.dot(ncross)))
                    a = b
                    prev_dist  = dist
                    prev_coords = [x_unif[j], y_unif[j]]

            angleDist = np.zeros(round(2 * np.pi * 10) + 1)
            for angle in angles:
                angleDist[int(math.floor((angle+np.pi)*10))] +=1
            #angleDist[0] +=angleDistTotal[-1]
            angleDist = angleDist[:-1]
            angleDistTotal = (i * angleDistTotal + angleDist)/(i + 1)


        MSD = (i * MSD + np.array(GetMSD(x_unif, y_unif))) / (i + 1)

        RMSDw = (i * RMSDw + np.array(GetRMSDW(x_unif, y_unif, RMSDw_l))) / (i + 1)

        peptide_remaining = (i*peptide_remaining+np.array(peptide_unif))/(i+1)
    plt.plot(angleDistTotal[1:-1])
    plt.show()
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



    if analytic!=-1:
        np.savez(filename, peptides = peptide_remaining, MSD =MSD, RMSDw = RMSDw, timescale = timescale, ptype=ptype, angle = angleDistTotal, Kuhn = KuhnTotal)
    else:
        np.savez(filename, peptides=peptide_remaining, MSD=MSD, RMSDw=RMSDw, timescale=timescale, ptype=ptype)










# See PyCharm help at https://www.jetbrains.com/help/pycharm/
