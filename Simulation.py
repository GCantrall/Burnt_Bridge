import matplotlib.pyplot as plt
from Particle import Particle
import numpy as np
import math

import random

class Simulation:
    def __init__(self, _tp, _lp, _tb,_tc,_pType="", _analytic = -1,_kb= 5):
        self.particle = Particle()
        self.x = 0.
        self.y = 0.
        self.tp = _tp
        self.td = 1
        self.lp = _lp
        self.ld = 1.
        self.peptide_size = 20
        self.tb = _tb
        self.peptide_remain = [600]
        self.neighbors = []
        self.pType = _pType
        self.analytic = _analytic
        self.energy = 0
        self.kb= _kb
        self.tc =_tc

    def CalculateEnergy(self,x2,y2):
        energy = 0
        for k in self.neighbors:
            dist = np.sqrt(math.pow(self.x_peptide[k] - x2, 2) + math.pow(self.y_peptide[k] - y2, 2))

            if (dist < self.peptide_size):
                if dist<.1:
                    energy = 10+energy
                else:
                    energy = 1/dist+energy
        return energy

    def PlotPath(self,x_tracker, y_tracker):

        fig, ax1 = plt.subplots()
        ax1.set_aspect('equal')
        ax1.scatter(self.x_peptide,self.y_peptide, c ='b',zorder=1)
        ax1.plot(x_tracker, y_tracker,c = 'r', zorder=0)

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


    def SetNeighbors(self):
        self.neighbors = []
        for k in range(len(self.x_peptide)):
            if(np.sqrt(math.pow(self.x_peptide[k]- self.x,2) + math.pow(self.y_peptide[k]- self.y,2))<(self.lp+2*self.peptide_size)):
                self.neighbors.append(k)


    def RunSimulationAttr(self, totalTime):
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

        if self.analytic != -1:
            lContour = self.ld*numDiff+self.lp*numRoll
            Kuhn = (np.power(self.x,2)+np.power(self.y,2))/lContour

            return time, x_tracker, y_tracker, [Kuhn, angle]
        else:
            return time, x_tracker, y_tracker

    def RunSimulationRep(self, totalTime):
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

        if self.analytic != -1:
            lContour = self.ld*numDiff+self.lp*numRoll
            Kuhn = (np.power(self.x,2)+np.power(self.y,2))/lContour

            return time, x_tracker, y_tracker, [Kuhn, angle]
        else:
            return time, x_tracker, y_tracker

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
                    previousLength = self.lp
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

        if self.analytic != -1:
            lContour = self.ld*numDiff+self.lp*numRoll
            Kuhn = (np.power(self.x,2)+np.power(self.y,2))/lContour

            return time, x_tracker, y_tracker, [Kuhn, angle]
        else:
            return time, x_tracker, y_tracker






def GetMSDR(x,y,inc):
    MSD = []
    for i in range(inc,len(x)):
        MSD.append(np.sqrt(pow(x[i]-x[i-inc],2)+pow(y[i]-y[i-inc],2)))

    return MSD