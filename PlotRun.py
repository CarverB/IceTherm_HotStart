#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 18:31:54 2016

@author: carverb
"""

from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np

def bulkdensity(r,rho):
    
    dr= np.gradient(r)
    mass=(r[0]+dr[0]/2.0)**3.0 * rho[0]
    
    for i in range(1,len(r)-1):
        mass+=((r[i]+dr[i]/2.0)**3.0 - (r[i]-dr[i]/2.0)**3.0)* rho[i]
    mass+=((r[-1])**3.0 - (r[-1]-dr[-1]/2.0)**3.0)* rho[-1]
    
    return (mass/r[-1]**3.0)

# Create a model run class
class ThermalOutput():
    def __init__(self, filename):
        self.params={}
        self.time=[]
        self.R=[]
        self.T=[]
        self.fmelt=[]
        self.density=[]
        self.flux=[]
        # load data
        self.ReadData(filename)
        
    # Read data file    
    def ReadData(self, filename):
        with open(filename, 'r') as infile:
            for line in infile:
                if ('&INP' in line or '/' in line):
                    continue
                elif ('=' in line): # Read Header values
                    line=line.replace(',','')
                    line=line.split('=')
                    self.params[line[0].strip()]=line[1] # place in dict
                else: # Read Data!
  #Log10 time (yrs),    R (m),    T (K),    Melt Fraction,    Density (kg/m^3),    Flux (W/m^2)
                    line=line.split()
                    self.time.append(10.0**float(line[0]))
                    self.R.append(float(line[1]))
                    self.T.append(float(line[2]))
                    self.fmelt.append(float(line[3]))
                    self.density.append(float(line[4]))
                    self.flux.append(float(line[5]))
                    
        # Reshape arrays
        #print self.params.keys()
        self.Nlayers=int(self.params['IMAX'])+int(self.params['IMAXC'])              
                    
        self.timeGyr=[x/1.0E9 for x in self.time[0:-1:self.Nlayers] ]    
        self.Ntime=int(len(self.time)/self.Nlayers)
        #self.Nlayers=len(self.R)/Ntime #get from header

        
        self.time2DGyr=np.reshape(self.time,(self.Nlayers,self.Ntime),'F')/1.0E9
        self.R2D=np.reshape(self.R,(self.Nlayers,self.Ntime),'F')/1000.0
        self.T2D=np.reshape(self.T,(self.Nlayers,self.Ntime),'F')
        self.fmelt2D=np.reshape(self.fmelt,(self.Nlayers,self.Ntime),'F')
        self.density2D=np.reshape(self.density,(self.Nlayers,self.Ntime),'F')
        self.flux2D=np.reshape(self.flux,(self.Nlayers,self.Ntime),'F')        

        self.coreI=int(self.params['IMAXC'])-1     
        
        self.OceanThick=np.sum(1000*np.gradient(self.R2D,axis=0)*self.fmelt2D, axis=0)      
        
        # calculate density as a function of time
        self.bulkdensity=np.zeros(self.Ntime)
        for i in range(self.Ntime):
            #print bulkdensity(self.R2D[:,i]*1000.0,
            #                     self.DR2D[:,i],self.density2D[:,i])
            self.bulkdensity[i]= bulkdensity(self.R2D[:,i]*1000.0,
                                 self.density2D[:,i])

        
    def plottemp(self):
        # Plot a 2D heatmap of the temperature
        plt.figure()
        plt.pcolor(self.time2DGyr, self.R2D, self.T2D)
        
        plt.colorbar()
        
        dmin=np.min(self.T2D)
        dmax=np.max(self.T2D)
        
        levels=np.arange(np.round(dmin/100)*100,np.round(dmax/100)*100,100)
        CL=plt.contour(self.time2DGyr, self.R2D, self.T2D, levels)
        plt.clabel(CL, fmt='%2.0f', colors='k', fontsize=14)

        
        
        # Plot a solid line a the core-ice boundary
        plt.plot(self.timeGyr, self.R2D[self.coreI,:], '-',
           color=[0.5020,0.3647,0.1294], linewidth=3)
           
       
         
        # Use try catch in case that level is all 0's
        try:   
            plt.contour(self.time2DGyr, self.R2D,self.fmelt2D,
                    levels=[.99], colors=['k'], linewidths=3)
        except:
            pass
        
        try:
            plt.contour(self.time2DGyr,self.R2D,self.porosity2D, 
                        levels=[1.0E-2], colors=['m'], linewidths=3)
        except:
            pass
        
        plt.title('Temperature')
        plt.ylabel('R (km)')
        plt.xlabel('Time (Gyr)')
        
        plt.axis('tight')
        
    def plotocean(self):
        # plot a heatmap of the melt fraction
        plt.figure()
        plt.pcolor(self.time2DGyr, self.R2D, self.fmelt2D)
        # Plot a solid line a the core-ice boundary
        plt.plot(self.timeGyr, self.R2D[self.coreI,:], '-',
           color=[0.5020,0.3647,0.1294], linewidth=3)
        plt.title('Melt Fraction')
        plt.ylabel('R (km)')
        plt.xlabel('Time (Gyr)')
        plt.colorbar()
        plt.axis('tight')

    def plotoceanthickness(self):
        # plot a line of the ocean thickness vs time
        plt.figure()
        plt.plot(self.timeGyr, self.OceanThick/1000.0, linewidth=3)
        #plt.title('')
        plt.ylabel('Ocean Thickness (km)')
        plt.xlabel('Time (Gyr)')
        #plt.colorbar()
        plt.axis('tight')
        
    def plotdensity(self):
        # Plot the density of each layer vs time
        plt.figure()
        plt.pcolor(self.time2DGyr, self.R2D, self.density2D)
        # Plot a solid line a the core-ice boundary
        plt.plot(self.timeGyr, self.R2D[self.coreI,:], '-',
           color=[0.5020,0.3647,0.1294], linewidth=3)
        plt.title('Density (kg/m^3)')
        plt.ylabel('R (km)')
        plt.xlabel('Time (Gyr)')
        plt.colorbar()
        plt.axis('tight')
        
    def plotbulkdensity(self):
        # Plot the bulk density of the body vs time
        plt.figure()
        plt.plot(self.timeGyr, self.bulkdensity, linewidth=3)
        #plt.title('')
        plt.ylabel('Bulk density (kg/m^3)')
        plt.xlabel('Time (Gyr)')
        #plt.colorbar()
        plt.axis('tight')
        
    def plotstrain(self):
        # Plot the percent strain vs time
        plt.figure()
        plt.plot(self.timeGyr, 100*(self.R2D[-1,:]**2/self.R2D[-1,0]**2-1)/2.0 , linewidth=3)
        
        plt.ylabel('Strain (%)')
        plt.xlabel('Time (Gyr)')
        plt.axis('tight')
        
    def plotradius(self):
        # plot body radius vs time
        plt.figure()
        plt.plot(self.timeGyr, self.R2D[-1,:], linewidth=3)
        #plt.title('')
        plt.ylabel('Radius (km)')
        plt.xlabel('Time (Gyr)')
        #plt.colorbar()
        plt.axis('tight')

    def plotflux(self):
        # Plot the flux in each layer as a function of time
        plt.figure()
        plt.pcolor(self.time2DGyr, self.R2D, self.flux2D)
        plt.plot(self.timeGyr, self.R2D[self.coreI,:], '-',
           color=[0.5020,0.3647,0.1294], linewidth=3)
        plt.title('Flux (W/m^2)')
        plt.ylabel('R (km)')
        plt.xlabel('Time (Gyr)')
        plt.colorbar()
        plt.axis('tight')

if __name__=='__main__':
    
    infile = 'IceTherm.out'

    Run=ThermalOutput(infile)

    #print('Ocean Thickness={:0.2f}'.format(Run.OceanThick[-1]/1000.0))    
    
    # Chose which plots to make
    Run.plottemp()
    #Run.plotocean()
    #Run.plotflux()
    Run.plotoceanthickness()
    Run.plotbulkdensity()
    Run.plotstrain()
    Run.plotradius()
    plt.show()


