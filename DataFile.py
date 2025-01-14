import numpy as np
import random
import math
import sys
from Simulation import Simulation
import matplotlib.pyplot as plt
import os.path

class DataSet:
    def __init__(self, replicates=1000, s_length=1000000, lp=20, tp=500, tb=4000, id=-1, version=-1, path = "",name="", simType = ""):
        self.replicates = replicates
        self.s_length = s_length
        self.lp = lp
        self.tp = tp
        self.path = path
        self.tb = tb
        self.id = id
        self.version = version
        self.DM = []
        self.peptides = []
        self.MSD = []
        self.RMSDw = []
        self.timescale = []
        self.angles = []
        self.Kuhn = []
        self.simType = simType
        if (name==""):
            self.name = path
        else:
            self.name = name


    def LoadData(self, fileIn = ""):

        if self.simType != "":
            type = self.simType+"_"
        else:
            type = ""
        if fileIn =="":
            filename = "Simulation_" +type + str(self.replicates) + "r_" + str(self.s_length) + "s_" + str(self.lp) + "lp_" + str(
                self.tp) + "tp_" + str(self.tb) + "tb"
            if (self.id != -1):
                filename = filename + "_" + str(self.id)
            if self.version != -1 and self.version!=1:
                filename = filename + "(" + str(self.version) + ")"
            filename = filename + ".npz"
            if(self.path!=""):
                filename = self.path+"/"+ filename
        else:
            filename = fileIn
        print(filename)
        file = np.load(filename, allow_pickle=True)
        self.peptides = file['peptides']
        self.MSD = file['MSD']
        self.RMSDw = file['RMSDw']
        self.timescale = file['timescale']
        if 'angle' in file.keys():
            self.angles = file['angle']
        if 'Kuhn' in file.keys():
            self.Kuhn  = file['Kuhn']
        if 'DM' in file.keys():
            self.DM = file['DM']


    def Average(self, DataSet2):
        self.MSD = (self.MSD*self.replicates+DataSet2.MSD*DataSet2.replicates)/(self.replicates+DataSet2.replicates)
        self.peptides = (self.peptides * self.replicates + DataSet2.peptides * DataSet2.replicates) / (self.replicates + DataSet2.replicates)
        self.replicates = self.replicates +DataSet2.replicates
        for Ku in DataSet2.Kuhn:
            np.append(self.Kuhn,Ku)
        #if len(self.DM)>1:
        #    self.DM = (self.DM * self.replicates + DataSet2.DM * DataSet2.replicates) / (
        #                self.replicates + DataSet2.replicates)
        if(len(self.angles)>1 and len(DataSet2.angles)>1):
            self.angles = (self.angles * self.replicates + DataSet2.angles * DataSet2.replicates) / (self.replicates + DataSet2.replicates)
class Trajectory:
    def __init__(self,_filename):
        self.x = []
        self.time = []
        self.y = []
        self.x_pep = []
        self.y_pep = []
        self.s_pep = []
        self.filename = _filename

    def LoadData(self):
        with open(self.filename, "r") as file:
            line = file.readline()
            while line:
                self.time.append(int(line[:-3]))

                line = file.readline()
                line = line[1:-2].split(',')
                self.x.append([float(z) for z in line])

                line = file.readline()
                line = line[1:-2].split(',')
                self.y.append([float(z) for z in line])

                line = file.readline()
                if line == "[]\n":
                    self.x_pep.append([])
                else:
                    line = line[1:-2].split(',')
                    self.x_pep.append([float(z) for z in line])

                line = file.readline()
                if line == "[]\n":
                    self.y_pep.append([])
                else:
                    line = line[1:-2].split(',')
                    self.y_pep.append([float(z) for z in line])

                line = file.readline()
                if line == "[]\n":
                    self.s_pep.append([])
                else:
                    line = line[1:-2].split(',')
                    self.s_pep.append([int(z) for z in line])

                line = file.readline()

    def PlotTrj(self,distance=-1,start= 0):
        numFigs = 5

        fig, ax = plt.subplots(1,numFigs)
        plt.subplots_adjust(wspace = 0.05)
        if start<0:
            start = len(self.x)+start-4
        x_min = np.min(np.array(self.x)[start:start+5])-5
        x_max = np.max(np.array(self.x)[start:start+5])+5
        x_dist = x_max-x_min
        y_min = np.min(np.array(self.y)[start:start+5])-5
        y_max = np.max(np.array(self.y)[start:start+5])+5
        y_dist = y_max-y_min
        print(x_dist)
        print(y_dist)
        if distance==-1:
            if y_dist>x_dist:
                x_min -= (y_dist-x_dist)/2
                x_max += (y_dist-x_dist)/2
            if x_dist>y_dist:
                y_min -= (x_dist-y_dist)/2
                y_max += (x_dist-y_dist)/2
        else:
            if distance>x_dist:
                x_min -= (distance-x_dist)/2
                x_max += (distance-x_dist)/2
            if distance>y_dist:
                y_min -= (distance-y_dist)/2
                y_max += (distance-y_dist)/2


        for i in range(start, start+5):


            ax[i-start].get_xaxis().set_visible(False)
            ax[i-start].get_yaxis().set_visible(False)
            ax[i-start].spines[:].set_linewidth(1)
            ax[i-start].set_xlim( x_min, x_max)
            ax[i-start].set_ylim(y_min, y_max)
            ax[i-start].set_aspect('equal')


            cax_2 = ax[i-start].scatter(self.x_pep[i], self.y_pep[i], zorder=1, c="red",
                                s=3 * np.array(self.s_pep[i]))
            for z in range(start,i):
                ax[i-start].plot(self.x[z], self.y[z], c='lightsteelblue', zorder=0)
            ax[i-start].plot(self.x[i], self.y[i], c= 'royalblue', zorder=0)
            #ax1.plot(x_end, y_end, c='blue', zorder=0)
            ax[i-start].scatter(self.x[i][-1], self.y[i][-1], s=30, c="black", zorder=3)
