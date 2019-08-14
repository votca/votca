#!/usr/bin/env python
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import numpy as np
import matplotlib.pyplot as plt
import csv
import re
import sys
import os
import argparse as ap
import numpy.linalg as lg


parser=ap.ArgumentParser(description="Making an .ogv video from a set of trajectory.csv files")
parser.add_argument("-t","--trajectories",nargs='+',required=True, help="trajectory.csv files to parse data from")
parser.add_argument('-f',"--frames",type=int,default=400, help="number of frames in movie")
parser.add_argument('-o',"--output",default="carriers.ogv", help="Name of the output video")
args=parser.parse_args()

class carrier(object):
    numberofobjects= 0
    def __init__(self):
        carrier.numberofobjects+=1
        self.id=carrier.numberofobjects
        print("Charge {} has been added".format(self.id))
        self.t=[]
        self.x=[]
        self.y=[]
        self.z=[]

    def tmax(self):
        return self.t[-1]

    def append(self,x,y,z,t):
        self.x.append(x) 
        self.y.append(y) 
        self.z.append(z) 
        self.t.append(t) 

    def returnpos(self,timestep):
        return self.pos[np.argmin(np.absolute(self.pos[:,0]-timestep))]

    def returntraj(self,timestep):
        return self.pos[0:np.argmin(np.absolute(self.pos[:,0]-timestep))+1]
    
    def returntrajdim(self):
        return np.array([np.amin(self.pos[:,1]),np.amax(self.pos[:,1]),np.amin(self.pos[:,2]),np.amax(self.pos[:,2]),np.amin(self.pos[:,3]),np.amax(self.pos[:,3])])     
    
    def array(self):
        self.pos=np.array([self.t,self.x,self.y,self.z]).T   

    def colorcharge(self,color):
        self.color=color   
     
    def info(self):
        print("Carrier No. ",self.id)                
        print(self.array())


class carriercollection(object):
    
    def __init__(self,steps=False):
        self.listofcarriers=[]
        self.foldername="kmc_video"
        self.steps=steps
        self.colors=['#610B0B','#B40404','#DF7401','#AEB404','#3ADF00','#01A9DB','#0101DF','#A901DB','#BE3E3E','#3326BF','#13D2BC','#1AD213','#AFD213','#1AD213','#66EA41','#D641EA']
    def appendcharge(self):
        charge=carrier()
        self.listofcarriers.append(charge)
        

    def addtimestep(self,row):
        for i,carrier in enumerate(self.listofcarriers):
            carrier.append(float(row[3*i+1]),float(row[3*i+2]),float(row[3*i+3]),float(row[0]))
                            
    
    def tmin(self):
        tmin=1E25
        for carrier in self.listofcarriers:
            ttest=carrier.tmax()
            if ttest<tmin:
                tmin=ttest
        return tmin


    def boxdimension(self):
            box=[]
            for carrier in self.listofcarriers:
                box.append(carrier.returntrajdim())
            boxsizes=np.array(box)
            self.xinterval=1.2*np.array([np.min(boxsizes[:,0]),np.max(boxsizes[:,1])])
            self.yinterval=1.2*np.array([np.min(boxsizes[:,2]),np.max(boxsizes[:,3])])
            self.zinterval=1.2*np.array([np.min(boxsizes[:,4]),np.max(boxsizes[:,5])])
            print(self.xinterval,self.yinterval,self.zinterval)
    def mergeonto(self,other):
        self.listofcarriers=self.listofcarriers+other.listofcarriers

    def assigncolors(self):
        i=0
        while len(self.listofcarriers)>len(self.colors):
            self.colors=self.colors+self.colors
        for carrier in self.listofcarriers:
            carrier.colorcharge(self.colors[i])
            i+=1
            
        
    def readinfile(self, infile):   
        with open (infile, "rb") as csvfile:
            print("Opening file {}".format(infile))
            reader = csv.reader(csvfile, dialect="excel-tab")
            for row in reader:
                if  reader.line_num>self.steps and self.steps!=False:
                    break
                if "time" in row[0]:
                    if "charges" in ''.join(row):
                        print(''.join(row))
                        m=re.search('\( ([0-9].?) charges \)',''.join(row))
                        noofcharges=int(m.group(1))
                    else:
                        noofcharges=len(row)/3

                    for i in range(noofcharges):
                        self.appendcharge()  

                else:
                        self.addtimestep(row)
            for carrier in self.listofcarriers:
                carrier.array()

            print("Found {} charges in file.".format(noofcharges))
    

    def plottimestep(self,time,step):
        #markers=r'$\ddot\smile$'
        markers=r'$\ddot\frown$'
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        timetext="t={:3e} s".format(time)
        print(timetext)
        ax.set_title(timetext,loc='center')
        ax.set_xlabel('x in [nm]')
        ax.set_ylabel('y in [nm]')
        ax.set_zlabel('z in [nm]')
        ax.set_xlim(self.xinterval[0],self.xinterval[1])
        ax.set_ylim(self.yinterval[0],self.yinterval[1])
        ax.set_zlim(self.zinterval[0],self.zinterval[1])
        for carrier in self.listofcarriers:
            posarray=carrier.returntraj(time)
            #print(posarray)
            ax.plot(posarray[:,1], posarray[:,2], posarray[:,3],c=carrier.color) 
            ax.scatter(posarray[-1,1], posarray[-1,2], posarray[-1,3],s=250,marker="o",c='#FFFF00')     
            ax.scatter(posarray[-1,1], posarray[-1,2], posarray[-1,3],s=150,marker=markers,c='#FFFF00')  
            ax.scatter(posarray[0,1], posarray[0,2], posarray[0,3],s=400,marker="x",c=carrier.color)
        plt.draw()
        filename = 'frame_{0:04d}.png'.format(step)
        plt.savefig(filename,dpi=100,bbox_inches="tight",pad_inches=0.0)
        plt.close() 
        return 0

    def createfolder(self,i):

        if os.path.isdir(self.foldername+str(i)) and (os.path.isdir(self.foldername+str(i-1)) or (os.path.isdir(self.foldername) and i==1)):
            os.rename(self.foldername+str(i),self.foldername+str(i+1))
            print("backing up folder {} to {}".format(self.foldername+str(i),self.foldername+str(i+1)))
        if os.path.isdir(self.foldername) and i==0:
            os.rename(self.foldername,self.foldername+str(i+1))
            print("backing up folder {} to {}".format(self.foldername,self.foldername+str(i+1)))
        if i==0:
            os.mkdir(self.foldername)
            print("creating folder {}".format(self.foldername))
        elif i>0:
            self.createfolder(i-1)
        else:
            print("Something about the foldernames is not right. Exiting....")
            sys.exit()

    def makemovie(self,filename,length):
        framerate=25
        self.assigncolors()
        self.boxdimension()
        os.chdir(self.foldername)
        tmin=self.tmin()
        print(tmin)
        steps=range(int(length*framerate)+1)
        for step in steps:
            print("evaluating frame {} of {}.".format(step,length*framerate))
            time=step*tmin/(length*framerate)
            self.plottimestep(time,step)
        command='ffmpeg -i frame_%04d.png {}'.format(filename)
        os.system(command)
        


step=args.frames
trajectories=carriercollection(step)

for csvfile in args.trajectories:
    print(csvfile)
    temptraj=carriercollection(step)
    temptraj.readinfile(csvfile)
    trajectories.mergeonto(temptraj)
    
trajectories.createfolder(20)
length=args.frames/25
trajectories.makemovie(args.output,length)
    

