import numpy as np
import random
import math
import sys
from Simulation import Simulation
import matplotlib.pyplot as plt
import os.path

class DataSet:
    def __init__(self, replicates=1000, s_length=1000000, lp=20, tp=100, tb=1000, id=-1, version=-1):
        self.replicates = replicates
        self.s_length = s_length
        self.lp = lp
        self.tp = tp
        self.tb = tb
        self.id = id
        self.version = version
        self.peptides = []
        self.MSD = []
        self.RMSDw = []
        self.timescale = []


    def LoadData(self):
        filename = "Simulation_" + str(self.replicates) + "r_" + str(self.s_length) + "s_" + str(self.lp) + "lp_" + str(
            self.tp) + "tp_" + str(self.tb) + "tb"
        if (self.id != -1):
            filename = filename + "_" + self.id
        if self.version != -1:
            filename = filename + "(" + str(self.version) + ")"
        filename = filename + ".npz"

        file = np.load(filename, allow_pickle=True)
        self.peptides = file['peptides']
        self.MSD = file['MSD']
        self.RMSDw = file['RMSDw']
        self.timescale = file['timescale']

    def Average(self, DataSet2):
        self.MSD = (self.MSD*self.replicates+DataSet2.MSD*DataSet2.replicates)/(self.replicates+DataSet2.replicates)
