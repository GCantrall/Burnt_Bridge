import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Particle import Particle
from Particle import Peptide
import numpy as np
import math
import subprocess
import os

import random

class Simulation:
    def __init__(self, _tp, _lp, _tb,_tc,_pType="", _analytic = -1,_kb= 1, _plot=-1):
        self.particle = Particle()
        self.x = 0.
        self.y = 0.
        self.tp = _tp
        self.td = 1
        self.lp = _lp
        self.ld = 1.
        self.peptide_size = 40
        self.tb = _tb
        self.peptide_remain = [600]
        self.neighbors = []
        self.pType = _pType
        self.analytic = _analytic
        self.energy = 0
        self.kb= _kb
        self.tc =_tc
        self.plot = _plot
        self.y_peptide_tracker =[]
        self.x_peptide_tracker = []
        self.strength_peptide_tracker = []


    # Check if Peptides are close enough to merge
    def CheckPeptides(self):

        for i in range(len(self.x_peptide)):
            for j in range(i+1,len(self.x_peptide)):
                if(j>len(self.x_peptide)-1):
                        break
                if (self.CalculatePeptideDistance(i,j)<np.sqrt(self.time_peptide[i]+self.time_peptide[j])) and self.strength_peptide[i]+self.strength_peptide[j]<11:

                    self.x_peptide[i] = round((self.x_peptide[i]+self.x_peptide[j])/2)
                    self.y_peptide[i] = round((self.y_peptide[i] + self.y_peptide[j]) / 2)
                    self.strength_peptide[i] = self.strength_peptide[i]+self.strength_peptide[j]
                    self.time_peptide[i] = 1
                    del self.x_peptide[j]
                    del self.y_peptide[j]
                    del self.strength_peptide[j]
                    del self.time_peptide[j]

        self.SetNeighbors()
        self.energy = self.CalculateEnergy(self.x, self.y)


    # Calculate the Distance between two peptides
    def CalculatePeptideDistance(self,i,j):
        return np.sqrt(pow(self.x_peptide[i]-self.x_peptide[j],2)+pow(self.y_peptide[i]-self.y_peptide[j],2))

    # Plot Final Trajectory
    def PlotPath(self,x_tracker, y_tracker):

        fig, ax1 = plt.subplots()
        ax1.set_aspect('equal')
        cax = ax1.scatter(self.x_peptide,self.y_peptide, c =self.strength_peptide,zorder=1, cmap='viridis')
        ax1.plot(x_tracker, y_tracker,c = 'r', zorder=0)
        cbar = fig.colorbar(cax, ticks=range(1,10))
        cbar.ax.set_yticklabels(range(1,10))

    # Plot Final Trajectory but only every 25 steps
    def PlotPathRange(self,time,x_tracker,y_tracker):
        x = []
        y = []
        fig, ax1 = plt.subplots()
        for i in range(len(time)):
            if time[i]%25 ==0:
                x.append(x_tracker[i])
                y.append(y_tracker[i])
        ax1.set_aspect('equal')
        ax1.scatter(self.x_peptide,self.y_peptide, c ='b',zorder=1)
        ax1.plot(x, y,c = 'r', zorder=0)


    # Creates imgages of the trajectorie
    def PlotPathVideo(self,time,x_tracker,y_tracker,x_peptide_unif,y_peptide_unif,strength_peptide_unif):
        x = []
        y = []
        folder= "testImages"
        fig, ax1 = plt.subplots()
        k = 1
        for i in range(len(time[-30000:])):
            j = i-30000
            if time[j]%25 ==0:
                x.append(x_tracker[j])
                y.append(y_tracker[j])
                fig, ax1 = plt.subplots()
                ax1.set_xlim(min(x_tracker[-30000:]),max(x_tracker[-30000:]))
                ax1.set_ylim(min(y_tracker[-30000:]), max(y_tracker[-30000:]))
                ax1.set_aspect('equal')
                cax  = ax1.scatter(x_peptide_unif[j], y_peptide_unif[j], c=strength_peptide_unif[j], zorder=1, vmin=1, vmax=10)
                ax1.plot(x, y, c='r', zorder=0)
                ax1.scatter(x[-1],y[-1],c="magenta")
                cbar = fig.colorbar(cax, ticks= range(1, np.amax(strength_peptide_unif[-1])))
                cbar.ax.set_yticklabels(range(1, np.amax(strength_peptide_unif[-1])))
                plt.savefig(folder + "/%03d.png" % k)
                k +=1
                plt.close(fig)

    # Creates imgages of the trajectorie
    def PlotPathVideoBeginEnd(self,time,x_tracker,y_tracker,x_peptide_unif,y_peptide_unif,strength_peptide_unif):
        x_ini = []
        y_ini = []
        x_end = []
        y_end = []
        folder= "testImages"
        fig, ax1 = plt.subplots()
        k = 1
        for i in range(30000):
            j = i-30000
            if time[j]%25 ==0:
                x_end.append(x_tracker[j]-x_tracker[-30000])
                y_end.append(y_tracker[j] - y_tracker[-30000])
                y_ini.append(y_tracker[i])
                x_ini.append([x_tracker[i]])

        fig, ax1 = plt.subplots()
        ax1.set_xlim(min([min(x_tracker[-30000:])-x_tracker[-30000],min(x_tracker[:30000])]), max([max(x_tracker[-30000:])-x_tracker[-30000],max(x_tracker[:30000])]))
        ax1.set_ylim(min([min(y_tracker[-30000:])-y_tracker[-30000],min(y_tracker[:30000])]), max([max(y_tracker[-30000:])-y_tracker[-30000],max(y_tracker[:30000])]))
        ax1.set_aspect('equal')

        #cax  = ax1.scatter(np.array(x_peptide_unif[j])-x_tracker[-30000], np.array(y_peptide_unif[j])-y_tracker[-30000], c=strength_peptide_unif[j], zorder=1, vmin=1, vmax=10,edgecolors="darkblue",cmap ="viridis")
        #cax_2  = ax1.scatter(x_peptide_unif[i], y_peptide_unif[i], c=strength_peptide_unif[i], zorder=1, vmin=1, vmax=10,edgecolors="red",cmap="viridis")

        ax1.plot(x_ini, y_ini, c='orange', zorder=0)
        ax1.plot(x_end, y_end, c='blue', zorder=0)
        #ax1.scatter(x_ini[-1],y_ini[-1], s=8,c="orange",edgecolors="black")
        #ax1.scatter(x_end[-1],y_end[-1], s=8, c="blue",edgecolors="black")


        #cbar = fig.colorbar(cax, ticks= range(1, 10))
        #cbar.ax.set_yticklabels(range(1, 10))
        plt.savefig(folder + "/%03d.png" % k)
        k +=1
        plt.close(fig)


    # Set Neighbors of Peptides
    def SetNeighbors(self):
        self.neighbors = []
        for k in range(len(self.x_peptide)):
            if(np.sqrt(math.pow(self.x_peptide[k]- self.x,2) + math.pow(self.y_peptide[k]- self.y,2))<(self.lp+2*self.peptide_size)):
                self.neighbors.append(k)

    def RunSimulation(self, totalTime):
        print("Error")

class BurntBridge(Simulation):
    def RunSimulation(self, totalTime):
        x_tracker = [0]
        y_tracker = [0]
        time  = [0]
        self.x = 0
        self.y = 0
        self.x_peptide = []
        self.y_peptide = []
        self.neighbors = []
        numDiff = 0
        numRoll = 0
        angle = []


        if(self.pType=="c"):
            self.particle.CreateFake()
        else:
            self.particle.CreateParticle()
        self.peptide_remain = [np.sum(self.particle.peptide)]
        current_location = 0
        vector = [1.,0.]


        while (time[-1]<totalTime):
            if(len(time)%30==0):
                self.SetNeighbors()
            withPeptide = []
            withoutPeptide = []
            options  = self.particle.GetEdges(current_location)
            for option in options:
                if (self.particle.peptide[option]==1):
                    withPeptide.append(option)
                else:
                    withoutPeptide.append(option)
            time.append(time[-1] + -math.log(random.random())/(len(withPeptide)/self.tp+1/self.td+len(withoutPeptide)/self.tb))

            rand = random.random()
            if(rand<(len(withPeptide)/(self.tp))/(len(withPeptide)/(self.tp)+1/self.td+len(withoutPeptide)/(self.tb))):
                self.SetNeighbors()
                choice  = random.random()*len(withPeptide)
                chosen = -1
                for i in range(len(withPeptide)):
                    if choice< (i+1):
                        chosen = i
                        break


                degree = self.particle.GetDirection(withPeptide[chosen])
                x2_d = math.cos(degree)*vector[0]-math.sin(degree)*vector[1]
                y2_d = math.sin(degree)*vector[0]+math.cos(degree)*vector[1]

                x2 = self.x+x2_d*self.lp
                y2 =self.y+y2_d*self.lp

                toClose = False
                for k in self.neighbors:
                    for j in range(self.lp):
                        if(np.sqrt(math.pow(self.x_peptide[k]- (self.x+x2_d*(j+1)),2) + math.pow(self.y_peptide[k]- (self.y+y2_d*(j+1)),2))<self.peptide_size):
                            toClose = True
                            break
                if not toClose:
                    degree, peptide = self.particle.MoveParticle(withPeptide[chosen])
                    vector = [x2_d,y2_d]

                    self.x_peptide.append(self.x)
                    self.y_peptide.append(self.y)
                    self.neighbors.append(len(self.x_peptide)-1)
                    self.x = x2
                    self.y = y2
                    numRoll = numRoll+1
                    current_location = withPeptide[chosen]

            elif(rand<(len(withPeptide)/(self.tp)+len(withoutPeptide)/(self.tb))/(len(withPeptide)/(self.tp)+1/self.td+len(withoutPeptide)/(self.tb))):
                choice  = random.random()*len(withoutPeptide)
                chosen = -1
                for i in range(len(withoutPeptide)):
                    if choice< (i+1):
                        chosen = i
                        break
                degree, peptide = self.particle.MoveParticle(withoutPeptide[chosen])
                x2 = math.cos(degree)*vector[0]-math.sin(degree)*vector[1]
                y2 = math.sin(degree)*vector[0]+math.cos(degree)*vector[1]
                vector = [x2,y2]
                current_location = withoutPeptide[chosen]
            else:
                degree = random.random()*2*math.pi
                x2 = self.x+self.ld*math.cos(degree)
                y2 = self.y+self.ld*math.sin(degree)
                toClose = False

                for k in self.neighbors:
                    if(np.sqrt(math.pow(self.x_peptide[k]- x2,2) + math.pow(self.y_peptide[k]- y2,2))<self.peptide_size):
                        toClose = True
                        break

                if not toClose:
                    self.x = x2
                    self.y = y2
                    numDiff = numDiff+1
                    previousLength = 1

            x_tracker.append(self.x)
            y_tracker.append(self.y)
            self.peptide_remain.append(np.sum(self.particle.peptide))

        #if self.analytic != -1:
        #    lContour = self.ld*numDiff+self.lp*numRoll
        #    Kuhn = (np.power(self.x,2)+np.power(self.y,2))/lContour

        #    return time, x_tracker, y_tracker, [Kuhn, angle]
        #else:
        return time, x_tracker, y_tracker

class RepulseToAttractPeptide(Simulation):

    def DiffusePeptide(self):
        for i in range(len(self.x_peptide)):
            degree = random.random() * 2 * math.pi
            self.x_peptide[i] = self.x_peptide[i] + 10 * math.cos(degree)/self.strength_peptide[i]
            self.y_peptide[i] = self.y_peptide[i] + 10 * math.sin(degree)/self.strength_peptide[i]


    def CalculateEnergy(self,x2,y2):
        energy = 0
        self.SetNeighbors()
        for k in self.neighbors:
            dist = np.sqrt(math.pow(self.x_peptide[k] - x2, 2) + math.pow(self.y_peptide[k] - y2, 2))
            energydiff = 0
            if (dist < self.peptide_size):
                if dist<.1:
                    energydiff = 10
                else:
                    energydiff = 1/dist
            energy += energydiff*(5-self.strength_peptide[k])
            energy -= energydiff
        return energy
    def RunSimulation(self, totalTime):
        x_tracker = [0]
        y_tracker = [0]
        time  = [0]
        self.x = 0
        self.y = 0
        self.x_peptide = []
        self.y_peptide = []
        self.neighbors = []
        previousLength = 1
        numDiff = 0
        numRoll = 0
        angle = []


        if(self.pType=="c"):
            self.particle.CreateFake()
        else:
            self.particle.CreateParticle()
        self.peptide_remain = [np.sum(self.particle.peptide)]
        current_location = 0
        vector = [1.,0.]



        while (time[-1]<totalTime):
            if(len(time)%30==0):
                self.SetNeighbors()
            withPeptide = []
            withoutPeptide = []

            options  = self.particle.GetEdges(current_location)
            for option in options:
                if (self.particle.peptide[option]==1):
                    withPeptide.append(option)
                else:
                    withoutPeptide.append(option)
            time.append(time[-1] + -math.log(random.random())/(len(withPeptide)/self.tp+1/self.td+len(withoutPeptide)/self.tb))

            rand = random.random()
            totalChance = (len(withPeptide) / (self.tp) +
                           1 / self.td +
                           len(withoutPeptide) / (self.tb)+
                           self.particle.peptide[current_location]/self.tc)
            if(rand<(self.particle.peptide[current_location]/self.tc)/totalChance):

                self.SetNeighbors()
                self.x_peptide.append(self.x)
                self.y_peptide.append(self.y)
                self.neighbors.append(len(self.x_peptide)-1)
                self.energy = self.energy+10
                self.particle.peptide[current_location] = 0

            elif (rand < (len(withPeptide) / (self.tp) +self.particle.peptide[current_location]/self.tc) /totalChance):
                choice  = random.random()*len(withPeptide)
                chosen = -1
                for i in range(len(withPeptide)):
                    if choice< (i+1):
                        chosen = i
                        break
                degree, peptide = self.particle.MoveParticle(withPeptide[chosen])
                x2 = math.cos(degree)*vector[0]-math.sin(degree)*vector[1]
                y2 = math.sin(degree)*vector[0]+math.cos(degree)*vector[1]
                vector = [x2,y2]
                current_location = withPeptide[chosen]

            elif(rand<(len(withPeptide)/(self.tp)+len(withoutPeptide)/(self.tb) +self.particle.peptide[current_location]/self.tc)/(totalChance)):
                choice  = random.random()*len(withoutPeptide)
                chosen = -1
                for i in range(len(withoutPeptide)):
                    if choice< (i+1):
                        chosen = i
                        break
                degree, peptide = self.particle.MoveParticle(withoutPeptide[chosen])
                x2 = math.cos(degree)*vector[0]-math.sin(degree)*vector[1]
                y2 = math.sin(degree)*vector[0]+math.cos(degree)*vector[1]
                vector = [x2,y2]
                current_location = withoutPeptide[chosen]
            else:
                degree = random.random()*2*math.pi

                x2 = self.x+(self.ld+3*(1-self.particle.peptide[current_location]))*math.cos(degree)
                y2 = self.y+(self.ld+3*(1-self.particle.peptide[current_location]))*math.sin(degree)

                energy_n = self.CalculateEnergy(x2,y2)

                if (energy_n-self.energy<=0 or random.random()<np.exp(-self.kb*(energy_n-self.energy))):
                    self.energy= energy_n
                    self.x = x2
                    self.y = y2
                    numDiff = numDiff+1
                    previousLength = 1

            x_tracker.append(self.x)
            y_tracker.append(self.y)
            self.peptide_remain.append(np.sum(self.particle.peptide))

        #if self.analytic != -1:
        #    lContour = self.ld*numDiff+self.lp*numRoll
        #    Kuhn = (np.power(self.x,2)+np.power(self.y,2))/lContour

        #    return time, x_tracker, y_tracker, [Kuhn, angle]
        #else:
        return time, x_tracker, y_tracker

class RepulsivePeptidesDirectedMotion(Simulation):

    def DiffusePeptide(self):
        for i in range(len(self.x_peptide)):
            degree = random.random() * 2 * math.pi
            self.x_peptide[i] = self.x_peptide[i] + 10 * math.cos(degree)/self.strength_peptide[i]
            self.y_peptide[i] = self.y_peptide[i] + 10 * math.sin(degree)/self.strength_peptide[i]

    def RunSimulation(self, totalTime):
            x_tracker = [0]
            y_tracker = [0]
            self.x_peptide_tracker = []
            self.y_peptide_tracker = []
            self.strength_peptide_tracker = []
            time  = [0]
            self.x = 0
            self.y = 0
            self.x_peptide = []
            self.y_peptide = []
            self.strength_peptide = []
            self.time_peptide = []
            self.neighbors = []
            self.kbd = 1
            numDiff = 0
            numRoll = 0
            angle = []


            if(self.pType=="c"):
                self.particle.CreateFake()
            else:
                self.particle.CreateParticle()
            self.peptide_remain = [np.sum(self.particle.peptide)]
            current_location = 0
            vector = [1.,0.]
            DM = -1

            roll_tracker = np.zeros(7)
            while (time[-1]<totalTime):
                if(len(time)%100==0):
                    self.DiffusePeptide()
                    self.CheckPeptides()
                elif(len(time)%30==0):
                    self.SetNeighbors()

                withPeptide = []
                withoutPeptide = []

                options  = self.particle.GetEdges(current_location)
                for option in options:
                    if (self.particle.peptide[option]==1):
                        withPeptide.append(option)
                    else:
                        withoutPeptide.append(option)
                deltaT =  -math.log(random.random())/(len(withPeptide)/self.tp+1/self.td+len(withoutPeptide)/self.tb)
                time.append(time[-1] + deltaT)
                #for i in range(len(self.time_peptide)):
                #    self.time_peptide[i] +=deltaT
                rand = random.random()
                totalChance = (len(withPeptide) / (self.tp) +
                               1 / self.td +
                               len(withoutPeptide) / (self.tb)+
                               self.particle.peptide[current_location]/self.tc)
                # Cleave Peptide
                if(rand<(self.particle.peptide[current_location]/self.tc)/totalChance):

                    self.SetNeighbors()
                    self.x_peptide.append(self.x)
                    self.y_peptide.append(self.y)
                    self.strength_peptide.append(1)
                    self.time_peptide.append(100)
                    self.neighbors.append(len(self.x_peptide)-1)
                    self.energy = self.CalculateEnergy(self.x, self.y)
                    self.particle.peptide[current_location] = 0
                    DM = -1
                    self.CheckPeptides()


                # Insert Peptide
                elif (rand < (len(withPeptide) / (self.tp) +self.particle.peptide[current_location]/self.tc) /totalChance):
                    choice  = random.random()*len(withPeptide)
                    chosen = -1
                    for i in range(len(withPeptide)):
                        if choice< (i+1):
                            chosen = i
                            break
                    degree, peptide = self.particle.MoveParticle(withPeptide[chosen])
                    x2 = math.cos(degree)*vector[0]-math.sin(degree)*vector[1]
                    y2 = math.sin(degree)*vector[0]+math.cos(degree)*vector[1]
                    if DM== -1:
                        DM =  random.random()*2*math.pi
                    vector = [x2,y2]
                    current_location = withPeptide[chosen]

                # Rotate To Other site
                elif(rand<(len(withPeptide)/(self.tp)+len(withoutPeptide)/(self.tb) +self.particle.peptide[current_location]/self.tc)/(totalChance)):
                    choice  = random.random()*len(withoutPeptide)
                    chosen = -1
                    for i in range(len(withoutPeptide)):
                        if choice< (i+1):
                            chosen = i
                            break
                    degree, peptide = self.particle.MoveParticle(withoutPeptide[chosen])
                    x2 = math.cos(degree)*vector[0]-math.sin(degree)*vector[1]
                    y2 = math.sin(degree)*vector[0]+math.cos(degree)*vector[1]
                    vector = [x2,y2]
                    DM = -1
                    current_location = withoutPeptide[chosen]
                # Diffuse
                else:
                    degree = random.random()*2*math.pi


                    roll_tracker[math.floor(degree)] +=1
                    #print(degree)


                    diff = np.abs(degree-DM)
                    if(diff> math.pi):
                        diff = math.pi*2 - diff
                    diff = diff/math.pi


                    if DM ==-1:
                        diff = 0
                    #if diff!= 0:
                    #    print("hi")
                    #if random.random()>np.exp(-self.kbd*diff):
                    #    continue


                    x2 = self.x+(self.ld)*math.cos(degree)
                    y2 = self.y+(self.ld)*math.sin(degree)

                    energy_n = self.CalculateEnergy(x2,y2)




                    if (energy_n-self.energy<=0 or random.random()<np.exp(-self.kb*(energy_n-self.energy)) ) and random.random()<np.exp(-self.kbd*diff):
                        self.energy= energy_n
                        self.x = x2
                        self.y = y2
                        numDiff = numDiff+1

                x_tracker.append(self.x)
                y_tracker.append(self.y)
                if self.plot==6:
                    self.y_peptide_tracker.append(self.y_peptide.copy())
                    self.x_peptide_tracker.append(self.x_peptide.copy())
                    self.strength_peptide_tracker.append(self.strength_peptide.copy())
                self.peptide_remain.append(np.sum(self.particle.peptide))

            #if self.analytic != -1:
            #    lContour = self.ld*numDiff+self.lp*numRoll
            #    Kuhn = (np.power(self.x,2)+np.power(self.y,2))/lContour

            #    return time, x_tracker, y_tracker, [Kuhn, angle]

            return time, x_tracker, y_tracker


    # Calculate The energy of peptides
    def CalculateEnergy(self,x2,y2):
        energy = 0
        self.SetNeighbors()
        for k in self.neighbors:
            dist = np.sqrt(math.pow(self.x_peptide[k] - x2, 2) + math.pow(self.y_peptide[k] - y2, 2))
            energydiff = 0
            if (dist < self.peptide_size):
                if dist<.1:
                    energydiff = 10
                else:
                    energydiff = 1/dist

            """
            if self.strength_peptide[k]>3:
                energy += -energydiff
            else:
                energy +=energydiff"""
            energydiff = energydiff*((self.strength_peptide[k]/10)*9+1)
            #energy +=0
            #energy += energydiff*(5-self.strength_peptide[k])/
            energy -= energydiff
        return energy


def GetMSDR(x,y,inc):
    MSD = []
    for i in range(inc,len(x)):
        MSD.append(np.sqrt(pow(x[i]-x[i-inc],2)+pow(y[i]-y[i-inc],2)))

    return MSD