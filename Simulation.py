from Particle import Particle
import numpy as np
import math

import random

class Simulation:
    def __init__(self, _tp, _lp, _tb):
        self.particle = Particle()
        self.x = 0.
        self.y = 0.
        self.tp = _tp
        self.td = 1
        self.lp = _lp
        self.ld = 1.
        self.tb = _tb
        self.peptide_remain = [600]

    def RunSimulation(self, totalTime):
        x_tracker = [0]
        y_tracker = [0]
        time  = [0]
        self.x = 0
        self.y = 0

        self.peptide_remain = [600]
        self.particle.CreateParticle()
        current_location = 0
        vector = [1.,0.]
        #self.particle.PlotParticle()


        while (time[-1]<totalTime):

            #if(len(time)==10000):
                #self.particle.PlotParticle()
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

                self.x = self.x+x2*self.lp
                self.y = self.y+ y2*self.lp
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
                self.x = self.x+self.ld*math.cos(degree)
                self.y = self.y+self.ld*math.sin(degree)
            x_tracker.append(self.x)
            y_tracker.append(self.y)
            self.peptide_remain.append(np.sum(self.particle.peptide))


        #plt.plot(x_tracker[94000:],y_tracker[94000:])
        #print(time[6000])
        #print(time[94000])
        #self.particle.PlotParticle()

        """
        plt.gca().set_aspect('equal')
        plt.plot(x_tracker[:6000],y_tracker[:6000])
        print(str(time[0])+" - "+str(time[6000]))
        plt.plot(np.array(x_tracker[-6000:])-x_tracker[-6000],np.array(y_tracker[-6000:])- y_tracker[-6000])
        print(str(time[-6000])+" - "+str(time[-1]))
        plt.show()
        """

        distance = []
        for j in range(len(x_tracker)):
            distance.append(math.sqrt(math.pow(x_tracker[j],2)+math.pow(y_tracker[j],2)))
        #self.PlotMSD(time, distance)
        return time, x_tracker, y_tracker







def GetMSDR(x,y,inc):
    MSD = []
    for i in range(inc,len(x)):
        MSD.append(np.sqrt(pow(x[i]-x[i-inc],2)+pow(y[i]-y[i-inc],2)))

    return MSD