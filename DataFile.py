import numpy as np
import random
import math
import sys
from Simulation import Simulation
import matplotlib.pyplot as plt
import os.path

class DataSet:
    def __init__(self, replicates=1000, s_length=1000000, lp=20, tp=500, tb=4000, id=-1, version=-1, sim="n", path = "",name=""):
        self.replicates = replicates
        self.s_length = s_length
        self.lp = lp
        self.tp = tp
        self.path = path
        self.tb = tb
        self.sim = sim
        self.id = id
        self.version = version
        self.peptides = []
        self.MSD = []
        self.RMSDw = []
        self.timescale = []
        self.angles = []
        self.Kuhn = []
        if (name==""):
            self.name = path
        else:
            self.name = name


    def LoadData(self):
        filename = "Simulation_"
        if self.sim !="n":
            filename = filename+self.sim+"_"

        filename = filename+ str(self.replicates) + "r_" + str(self.s_length) + "s_" + str(self.lp) + "lp_" + str(
            self.tp) + "tp_" + str(self.tb) + "tb"
        if (self.id != -1):
            filename = filename + "_" + str(self.id)
        if self.version != -1 and self.version!=1:
            filename = filename + "(" + str(self.version) + ")"
        filename = filename + ".npz"
        if(self.path!=""):
            filename = self.path+"/"+ filename
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


    def Average(self, DataSet2):
        self.MSD = (self.MSD*self.replicates+DataSet2.MSD*DataSet2.replicates)/(self.replicates+DataSet2.replicates)
        self.peptides = (self.peptides * self.replicates + DataSet2.peptides * DataSet2.replicates) / (self.replicates + DataSet2.replicates)
        self.replicates = self.replicates +DataSet2.replicates
        for Ku in DataSet2.Kuhn:
            np.append(self.Kuhn,Ku)
        if(len(self.angles)>1 and len(DataSet2.angles)>1):
            self.angles = (self.angles * self.replicates + DataSet2.angles * DataSet2.replicates) / (self.replicates + DataSet2.replicates)
