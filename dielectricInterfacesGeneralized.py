import matplotlib.pyplot as plt
import numpy as np
from math import inf
from math import isinf

def testFunctionality():

    d = buildDomain([1000000,1.2,1.3],[0,2],[1],[1],10)
    
    d.checkSolutionAccuracy(.00001)


    #Plot all charges in the set
    [q,x,y,z] = d.returnImageCharges()
    plt.plot(z[1],q[1],'bo')
    for region in d.regions:
        if(not isinf(region.zleft)):
            plt.axvline(x=region.zleft,linestyle=':')
    plt.xlabel('Position')
    plt.ylabel('Charge')
    plt.show()

    #Plot the potential
    [z,v] = d.getPotentialFull(30)
    plt.plot(z,v,'b')
    plt.xlabel('Position')
    plt.ylabel('Potential')
    plt.show()
    

#Computes total charge induced at the surfaces of one or more regions
def getNetCharge(domain,regionIndices):
    qtotal = 0
    for ind in regionIndices:
        qtotal = qtotal + domain.regions[ind].leftSurfaceCharge
        if(ind+1 not in regionIndices):
            qtotal = qtotal + domain.regions[ind].rightSurfaceCharge
    return(qtotal)
   

def buildDomain(eps,boundaries,qinit,zinit,numIterations):
    if(len(eps)-1 != len(boundaries)):
        print("ERROR: INPROPER SPECFICATION OF REGIONS")
    d = domain(eps,boundaries)
    for i in range(len(qinit)):
        d.addPointCharge(qinit[i],0,0,zinit[i])

    for i in range(numIterations):
        d.iterate()

    return(d)
    

class pointCharge():
    #destLeft/destRight indicate whether this charge is to be reflected across the left/right boundaries in the next iteration
    def __init__(self,q,x,y,z,destLeft,destRight,freeCharge):
        self.q = q
        self.x = x
        self.y = y
        self.z = z
        self.destLeft = destLeft
        self.destRight = destRight
        self.freeCharge = freeCharge

    def dist(self,x1,y1,z1,x2,y2,z2):
        return(np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2))
    
    def potentialPoint(self,epsilon,q,x1,y1,z1,x2,y2,z2):
        return((1/(4*np.pi*epsilon))*q/self.dist(x1,y1,z1,x2,y2,z2))

    def potential(self,epsilon,xt,yt,zt):
        return(self.potentialPoint(epsilon,self.q,self.x,self.y,self.z,xt,yt,zt))


    


class region():
    def __init__(self,pointCharges,epsilon,zleft,zright):
        self.pointCharges = pointCharges
        self.epsilon = epsilon
        self.zleft = zleft
        self.zright = zright
        self.leftSurfaceCharge = 0
        self.rightSurfaceCharge = 0

    def potential(self,xt,yt,zt):
        v = 0
        for p in self.pointCharges:
            v = v + p.potential(self.epsilon,xt,yt,zt)

        return(v)

class domain():
    def __init__(self,eps,boundaries):
        self.regions = [region(list(),eps[0],-inf,boundaries[0])]
        for i in range(len(eps))[1:-1]:
            self.regions.append(region(list(),eps[i],boundaries[i-1],boundaries[i]))
        self.regions.append(region(list(),eps[-1],boundaries[-1],inf))

    #Used for adding free charges prior to image charges being generated
    def addPointCharge(self,qinit,xinit,yinit,zinit):
        for region in self.regions:
            if(region.zleft < zinit <= region.zright):
                destLeft = region.zleft != -inf
                destRight = region.zright != inf
                region.pointCharges.append(pointCharge(qinit,xinit,yinit,zinit,destLeft,destRight,True))

    def iterate(self):
        #At each iteration new charges are not added to the underlying class until the end of the iteration
        #This prevents multiple reflections from happening at the same iteration in an asymmetric manner
        #Charges to be added at the end are stored here
        chargesAdded = list()
        for i in range(len(self.regions)):
            chargesAdded.append(list())
            
        for i in range(len(self.regions)):
            for p in self.regions[i].pointCharges:
                #Check if charge is to be reflected over left boundary
                if(p.destLeft):
                    #Add charge for current region
                    qNew = -p.q*(self.regions[i-1].epsilon - self.regions[i].epsilon)/(self.regions[i-1].epsilon + self.regions[i].epsilon)
                    zNew = self.regions[i].zleft - (p.z - self.regions[i].zleft)
                    newDestRight = self.regions[i].zright != inf
                    chargesAdded[i].append(pointCharge(qNew,p.x,p.y,zNew,False,newDestRight,False))
                    #Update total surface charges for both regions at boundary
                    self.regions[i].leftSurfaceCharge = self.regions[i].leftSurfaceCharge + (qNew/self.regions[i].epsilon)
                    self.regions[i-1].rightSurfaceCharge = self.regions[i-1].rightSurfaceCharge + (qNew/self.regions[i].epsilon)
                    #Add charge for region to the left (dest depends on left region)
                    qNew = p.q*(2*self.regions[i-1].epsilon)/(self.regions[i-1].epsilon + self.regions[i].epsilon)
                    newLeftDest = self.regions[i-1].zleft != -inf
                    chargesAdded[i-1].append(pointCharge(qNew,p.x,p.y,p.z,newLeftDest,False,False))
                    
                #Check if charge is to be reflected over right boundary
                if(p.destRight):
                    #Add charge for current region
                    qNew = -p.q*(self.regions[i+1].epsilon - self.regions[i].epsilon)/(self.regions[i+1].epsilon + self.regions[i].epsilon)
                    zNew = self.regions[i].zright + (self.regions[i].zright - p.z)
                    newDestLeft = self.regions[i].zleft != -inf
                    chargesAdded[i].append(pointCharge(qNew,p.x,p.y,zNew,newDestLeft,False,False))
                    #Update total surface charges for both regions at boundary
                    self.regions[i].rightSurfaceCharge = self.regions[i].rightSurfaceCharge + (qNew/self.regions[i].epsilon)
                    self.regions[i+1].leftSurfaceCharge = self.regions[i+1].leftSurfaceCharge + (qNew/self.regions[i].epsilon)
                    #Add charge for region to the right (dest dependes on right region)
                    qNew = p.q*(2*self.regions[i+1].epsilon)/(self.regions[i+1].epsilon + self.regions[i].epsilon)
                    newDestRight = self.regions[i+1].zright != inf
                    chargesAdded[i+1].append(pointCharge(qNew,p.x,p.y,p.z,False,newDestRight,False))

                #At this point all necessary reflections have been made, so the charge is updated to avoid duplicated reflections on the next iteration
                p.destLeft = False
                p.destRight = False

        for i in range(len(self.regions)):
            self.regions[i].pointCharges = self.regions[i].pointCharges + chargesAdded[i]

    def checkSolutionAccuracy(self,dz):
        numInteriorCharges = list()
        for i in range(len(self.regions)):
            numCharges = 0
            for j in range(len(self.regions[i].pointCharges)):
                if(self.regions[i].zleft <= self.regions[i].pointCharges[j].z <= self.regions[i].zright):
                    numCharges = numCharges + 1
            numInteriorCharges.append(numCharges)
        print('Number of Charges in Regions:  ' + str(numInteriorCharges))
                    
        for i in range(len(self.regions))[:-1]:
            vLeft = self.regions[i].potential(0,0,self.regions[i].zright)
            vRight = self.regions[i+1].potential(0,0,self.regions[i+1].zleft)
            dvLeft = (self.regions[i].potential(0,0,self.regions[i].zright) - self.regions[i].potential(0,0,self.regions[i].zright-dz))/dz
            dvRight = (self.regions[i+1].potential(0,0,self.regions[i+1].zleft + dz) - self.regions[i+1].potential(0,0,self.regions[i+1].zleft))/dz

            print('Continuity at boundary ' + str(i) + '/' + str(i+1))
            print(str(vLeft))
            print(str(vRight))
            print('Derivative at boundary ' + str(i) + '/' + str(i+1))
            print(self.regions[i].epsilon*dvLeft)
            print(self.regions[i+1].epsilon*dvRight)

           
    def printNumCharges(self):
        print('\n')
        s = str(len(self.regions[0].pointCharges))
        for i in range(len(self.regions))[1:]:
            s = s + " " + str(len(self.regions[i].pointCharges))
            #print(len(self.regions[i].pointCharges))
        print(s)


    def getPotentialFull(self,pointsPerRegion):
        z = list()
        v = list()
        for i in range(len(self.regions)):
            if(i == 0):
                znew = np.linspace(self.regions[i].zright-(self.regions[i+1].zright - self.regions[i+1].zleft),self.regions[i].zright,pointsPerRegion)
            elif(i == len(self.regions)-1):
                znew = np.linspace(self.regions[i].zleft,self.regions[i].zleft+(self.regions[i-1].zright - self.regions[i-1].zleft),pointsPerRegion)
            else:
                znew = np.linspace(self.regions[i].zleft,self.regions[i].zright,pointsPerRegion)

            vnew = [self.regions[i].potential(0,0,j) for j in znew]

            z = z + list(znew)
            v = v + list(vnew)
        return(z,v)
                

    #As intended, this function does NOT include free charges added at the beginning, only image charges 
    def returnImageCharges(self):
        qAll = list()
        xAll = list()
        yAll = list()
        zAll = list()

        for i in range(len(self.regions)):
            qNew = list()
            xNew = list()
            yNew = list()
            zNew = list()
            for p in self.regions[i].pointCharges:
                #if(not p.freeCharge):
                if(True):
                    qNew.append(p.q)
                    xNew.append(p.x)
                    yNew.append(p.y)
                    zNew.append(p.z)

            qAll.append(qNew)
            xAll.append(xNew)
            yAll.append(yNew)
            zAll.append(zNew)

        return(qAll,xAll,yAll,zAll)
