import numpy as np
import random
import matplotlib.pyplot as plt
import math

class Particle:

    def __init__(self):
        self.edges = []
        self.coordinates = []
        self.peptide = np.ones(600)
        self.current = 0
        self.previous = 1

    def GetEdges(self,i):
        return self.edges[i]

    def MoveParticle(self, new):

        degree = self.GetDirection(new)
        peptide = self.peptide[new]
        self.peptide[self.current] = 0
        self.previous = self.current
        self.current = new

        return (degree, peptide)

    def GetDirection(self, new):

        co = self.coordinates[self.previous]

        cc = self.coordinates[self.current]
        cn = self.coordinates[new]
        if self.previous == new:
            return math.pi


        xo = math.sin(co[1]) * math.cos(co[0])
        yo = math.sin(co[1]) * math.sin(co[0])
        zo = math.cos(co[1])

        xc = math.sin(cc[1]) * math.cos(cc[0])
        yc = math.sin(cc[1]) * math.sin(cc[0])
        zc = math.cos(cc[1])


        xn = math.sin(cn[1]) * math.cos(cn[0])
        yn = math.sin(cn[1]) * math.sin(cn[0])
        zn = math.cos(cn[1])




        ko = -(xo*xc+yo*yc+zo*zc)/(pow(zc,2)+pow(xc,2)+pow(yc,2))
        kn = -(xn*xc+yn*yc+zn*zc)/(pow(zc,2)+pow(xc,2)+pow(yc,2))

        xo_p = xo+ko*xc
        yo_p = yo+ko*yc
        zo_p = zo+ko*zc

        xn_p = xn+kn*xc
        yn_p = yn+kn*yc
        zn_p = zn+kn*zc


        """
        fig = plt.figure()

        # add axes
        ax = fig.add_subplot(111,projection='3d')

        xx, yy = np.meshgrid(np.arange(-1,1,.1), np.arange(-1,1,.1))
        z = -(xc*(xx-xc)+yc*(yy-yc))/zc+zc
        ax.plot_surface(xx, yy, z, alpha=0.5)

        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x2 = np.cos(u)*np.sin(v)
        y2 = np.sin(u)*np.sin(v)
        z2 = np.cos(v)
        ax.plot_wireframe(x2, y2, z2, color="r")
        ax.set_zlim((-1.5, 1.5))
        # plot the plane
        ax.scatter(xo,yo,zo,c='r')
        ax.scatter(xn,yn,zn,c='r')
        ax.scatter(xc,yc,zc,c='g')
        ax.scatter(xo_p,yo_p,zo_p,c='b')
        ax.scatter(xn_p,yn_p,zn_p,c='b')
        

        plt.show()"""


        degree = math.acos((xo_p*xn_p+yo_p*yn_p)/(math.sqrt(math.pow(xo_p,2)+math.pow(yo_p,2))*math.sqrt(math.pow(xn_p,2)+math.pow(yn_p,2))))+np.pi
        return degree

    def MaxCoords(self):
        maxx= 0
        maxy = 0
        for i in range(len(self.coordinates)):
            if self.coordinates[i][0] > maxx:
                maxx = self.coordinates[i][0]
            if self.coordinates[i][1] > maxy:
                maxy = self.coordinates[i][0]
        if( maxx >6.5):
            print("hi")
        print("x: "+str(maxx)+ "   y: "+str(maxy))

    def GetDirection_f(self, new):
        co = [math.pi+.07,math.pi+1]# self.coordinates[self.previous]
        cc = [0, 0]#self.coordinates[self.current]
        cn = [math.pi-1, math.pi/2-1]#self.coordinates[new]


        xo = math.sin(co[1]) * math.cos(co[0])
        yo = math.sin(co[1]) * math.sin(co[0])
        zo = math.cos(co[1])

        xc = math.sin(cc[1]) * math.cos(cc[0])
        yc = math.sin(cc[1]) * math.sin(cc[0])
        zc = math.cos(cc[1])

        xn = math.sin(cn[1]) * math.cos(cn[0])
        yn = math.sin(cn[1]) * math.sin(cn[0])
        zn = math.cos(cn[1])

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')


        ax.scatter(xo, yo, math.cos(co[1]),c='r')
        ax.scatter(xn, yn, math.cos(cn[1]),c='r')
        ax.scatter(xc, yc, math.cos(cc[1]), c='b')

        co[0] = co[0] - cc[0]
        co[1] = co[1] - cc[1]
        cn[0] = cn[0] - cc[0]
        cn[1] = cn[1] - cc[1]

        xo = math.sin(co[1]) * math.cos(co[0])
        yo = math.sin(co[1]) * math.sin(co[0])

        xn = math.sin(cn[1]) * math.cos(cn[0])
        yn = math.sin(cn[1]) * math.sin(cn[0])
        ax.scatter(0, 0, 0, c='m')
        ax.scatter(xo, yo, 0,c='g')
        ax.scatter(xn, yn, 0,c='g')


        u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        x = np.cos(u) * np.sin(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(v)
        ax.plot_wireframe(x, y, z, color="black")

        plt.show()

    def CreateParticle(self):

        self.peptide = np.ones(600)
        self.current = 0
        self.previous = 1
        self.coordinates = []
        self.coordinates.append([2*math.pi*random.random(),2*math.pi*random.random()])
        self.edges = []
        i = 0
        maxDelta = .1
        while i < 599:
            theta = 2*math.pi*random.random()
            phi = 2*math.pi*random.random()
            toClose = False
            for coord in self.coordinates:
                #print(self.GetDistance(coord,coord))
                if (self.GetDistance([theta,phi], coord)<maxDelta):
                    g = self.GetDistance([theta,phi],coord)
                    toClose = True
                    break
            if toClose:
                continue
            else:
                self.coordinates.append([theta,phi])
                i = i+1
        for i in range(len(self.coordinates)):
            edgeList = []
            for j in range(len(self.coordinates)):
                if i==j:
                    continue
                c1 = self.coordinates[i]
                c2 = self.coordinates[j]
                if ((self.GetDistance(c1,c2 )))<.2:
                   edgeList.append(j)
            self.edges.append(edgeList)

    def CreateFake(self):
        count  = 1
        self.edges = []
        self.coordinates = []
        self.current = 0
        self.previous = 1
        for i in np.arange(0,2*math.pi,.1):
            #self.coordinates.append([0,i])
            self.coordinates.append([i, math.pi/2])
            edgelist = [count]
            count = count + 1

            self.edges.append(edgelist)
        self.peptide = np.ones(count-1)
        self.edges[-1] = [0]



    def GetDistance(self, coord1, coord2):

        x1= (math.sin(coord1[1])*math.cos(coord1[0]))
        y1 = (math.sin(coord1[1])*math.sin(coord1[0]))
        z1 = (math.cos(coord1[1]))

        x2 = math.sin(coord2[1])*math.cos(coord2[0])
        y2 = math.sin(coord2[1])*math.sin(coord2[0])
        z2 = (math.cos(coord2[1]))

        return math.sqrt(math.pow(x2-x1,2)+math.pow(y2-y1,2)+math.pow(z2-z1,2))

    def GetPeriodicDistance(self,c1,c2):
        diff = abs(c1-c2)
        if diff>math.pi:
            diff = 2*math.pi-diff
        return diff


    def PlotParticle(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        x_r = []
        x_u = []
        y_r = []
        y_u = []
        z_r = []
        z_u = []
        for i in range(len(self.coordinates)):
            coord = self.coordinates[i]
            if self.peptide[i] == 1:
                x_r.append(math.sin(coord[1])*math.cos(coord[0]))
                y_r.append(math.sin(coord[1])*math.sin(coord[0]))
                z_r.append(math.cos(coord[1]))
            else:
                x_u.append(math.sin(coord[1])*math.cos(coord[0]))
                y_u.append(math.sin(coord[1])*math.sin(coord[0]))
                z_u.append(math.cos(coord[1]))

        ax.scatter(x_r,y_r,z_r,c='blue')
        ax.scatter(x_u,y_u, z_u, c='r')
        ax.set_axis_off()
        plt.show()

    def Plot2D(self):
        x = []
        y = []
        for coord in self.coordinates:
            x.append(coord[0])
            y.append(coord[1])

        plt.scatter(x,y)
        plt.show()
