import numpy as np
import random
import math
from Particle import Particle
import matplotlib.pyplot as plt






class Simulation:
    def __init__(self):
        self.particle = Particle()
        self.x = 0.
        self.y = 0.
        self.tp = 150
        self.td = 1.
        self.lp = 30.
        self.ld = 1.
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
            time.append(time[-1] + -math.log(random.random())/(1/self.tp+1/self.td))
            if(random.random()<((1/self.tp)/(1/(self.tp)+1/self.td))):
                #print("large")
                options  = self.particle.GetEdges(current_location)
                choice  = random.random()*len(options)
                chosen = -1
                for i in range(len(options)):
                    if choice< (i+1):
                        chosen = i
                        break
                peptide = self.particle.peptide[options[chosen]]
                if(peptide==1):
                    degree, peptide = self.particle.MoveParticle(options[chosen])
                    x2 = math.cos(degree)*vector[0]-math.sin(degree)*vector[1]
                    y2 = math.sin(degree)*vector[0]+math.cos(degree)*vector[1]
                    vector = [x2,y2]
                    if peptide==1:
                        self.x = self.x+x2*self.lp
                        self.y = self.y+ y2*self.lp
                    current_location = options[chosen]
            else:
                #print("small")
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


    def PlotMSD(self, times, distances):
        msd = []
        for distance in distances:
            msd.append(math.pow(distance,2))

        plt.plot(times,msd)
        plt.show()




def GetMSDR(x,y,inc):
    MSD = []
    for i in range(inc,len(x)):
        MSD.append(np.sqrt(pow(x[i]-x[i-inc],2)+pow(y[i]-y[i-inc],2)))

    return MSD
#x = particle.GetDirection(8)

#particle.CreateFake()
#particle.PlotParticle()

#particle.Plot2D()
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    S = Simulation()
    simLength = 1000000
    MSDLen=1
    #x = S.particle.GetDistance([0,math.pi/2], [0,0])
    times = []

    p = np.arange(0,1,.1)
    l = np.arange(1,simLength)
    timescale = np.concatenate((p,l))
    distances = np.zeros(len(timescale))
    peptide_remaining = np.zeros(len(timescale))
    x = np.zeros(len(timescale))
    y = np.zeros(len(timescale))
    MSD = np.zeros(len(timescale)- MSDLen)

    for i in range(100):
        print(i)

        time, x_tracker, y_tracker = S.RunSimulation(simLength)

        distance = []
        for j in range(len(x_tracker)):
            distance.append(math.sqrt(math.pow(x_tracker[j],2)+math.pow(y_tracker[j],2)))


        distance_unif = []
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


            distance_unif.append(distance[k])
            peptide_unif.append(S.peptide_remain[k])
            x_unif.append(x_tracker[k])
            y_unif.append(y_tracker[k])

        peptide_remaining = (i*peptide_remaining+np.array(peptide_unif))/(i+1)
        distances = (i*distances+np.array(distance_unif))/(i+1)
        MSD = (i*MSD+np.array(GetMSDR(x_unif,y_unif,MSDLen)))/(i+1)


    distances = np.power(np.array(distances),2)


    np.savez("InfoFile_1000_1", peptides = peptide_remaining, distance = distances, WMSD = MSD, timescale = timescale)


    #file = np.load("InfoFile.npz",allow_pickle=True)

    #test = file.files
    """
    print("hi")
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


    plt.show()"""









# See PyCharm help at https://www.jetbrains.com/help/pycharm/
