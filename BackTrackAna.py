#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to compute eigenvectors, powers, etc. for TMs


@author: wichmann
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import scipy.sparse.linalg as sp_linalg
#import colormaps as cm
from matplotlib import cm as cm2
import datetime
import matplotlib.animation as animation
from AnaObjects import oceanvector
#from numba import jit
#import matplotlib.colors as colors
from scipy import interpolate


datadir = '/Users/wichmann/sources/inertialparticles_David/Simulations/Run1Output/'
plotdir = '/Users/wichmann/sources/inertialparticles_David/Simulations/Run1Output/Plots/'


files02deg = {'P': 'P_Run1KernelP_gammaNone_wriseNone_gridnum',
              'S': 'S_Run1KernelS_gammaNone_wriseNone_gridnum',
              'I1': 'I1_Run1KernelI_gamma5.0_wriseNone_gridnum',
              'I2': 'I2_Run1KernelI_gamma15.0_wriseNone_gridnum',
              }

files5deg = {'P': 'TestRunPKernelP_gammaNone_wriseNone_gridnum',
             'S': 'TestRunSKernelS_gammaNone_wriseNone_gridnum',
             'I1': 'TestRunI1KernelI_gamma5.0_wriseNone_gridnum',
             'I2': 'TestRunI2KernelI_gamma15.0_wriseNone_gridnum',
             'M1': 'TestRunM1KernelM_gammaNone_wrise0.01_gridnum',
             'M2': 'TestRunM2KernelM_gammaNone_wrise0.001_gridnum',
             'M3': 'TestRunM3KernelM_gammaNone_wrise0.03_gridnum'}


files10deg = {'I': 'TestKernelI_gamma5.0_gridnum', 'P': 'TestKernelP_gammaNone_gridnum', 'M': 'TestKernelM_gammaNone_wrise0.0001_gridnum','S': 'TestKernelS_gammaNone_gridnum'} #, 'PN': 'TrajRun1KernelPN_tau1.0_D10.0_gridnum'}

files =files02deg

#For binned plots
ddeg =2.0
minlon=0.
maxlon=360.
minlat=-75.
maxlat=75.

Lons=np.arange(minlon,maxlon,ddeg)
Lats=np.arange(minlat,maxlat,ddeg)
lon_bins=np.linspace(minlon,maxlon,int((maxlon-minlon)/ddeg)+1)        
lat_bins=np.linspace(minlat,maxlat,int((maxlat-minlat)/ddeg)+1)        

#Number of different grid arrays
N=24

class ParticleData(object):
    def __init__(self, pdir, fname, tmin=-2, tmax=-1,indices=None,load_initial=True):
        print 'Setting up ParticleData object for pdir: ', pdir
        
        #Load data
        i = 0
        pfile = pdir + fname + str(i)+'.nc'     
        data = Dataset(pfile,'r')
        print 'Load grid no: ', i
        time=data.variables['time'][:,tmin:tmax]
        lon=data.variables['lon'][:,tmin:tmax]
        lat=data.variables['lat'][:,tmin:tmax]
        
        for i in range(1,N):
            pfile = pdir + fname + str(i)+'.nc'     
            data = Dataset(pfile,'r')
            print 'Load grid no: ', i
            time=data.variables['time'][:,tmin:tmax]
            lon=np.vstack((lon, data.variables['lon'][:,tmin:tmax]))
            lat=np.vstack((lat, data.variables['lat'][:,tmin:tmax]))
        time/=86400. #Convert to days
        time-=time[0,0]

        Time=[]        

        for t in range(len(time[0])):
            for i in range(len(time)):
                if not np.ma.is_masked(time[i,t]):
                    Time.append(time[i,t])
                    break
        
        self.Time = Time
        
        if indices is not None: #Only several particles (e.g. from coasts) are considered
            lon=lon[indices]
            lat=lat[indices]
            time=time[indices]
            
        self.lat = lat
        self.lon = lon
        self.time = time
        self.time_origin = datetime.datetime(2000, 1, 1) #When particles were released

        if load_initial:
            print 'Setting distribution at time zero'        
            i = 0
            pfile = pdir + fname + str(i)+'.nc'
            data = Dataset(pfile,'r')
            lon0=data.variables['lon'][:,0]
            lat0=data.variables['lat'][:,0]        
            for i in range(1,N):
                pfile = pdir + fname + str(i)+'.nc'
                data = Dataset(pfile,'r')
                lon0=np.append(lon0, data.variables['lon'][:,0])
                lat0=np.append(lat0, data.variables['lat'][:,0])
            
            self.lons0=lon0
            self.lats0=lat0
            distr , _, _ = np.histogram2d(lat0, lon0, [lat_bins, lon_bins])
            D=oceanvector(distr,minlon=minlon,maxlon=maxlon,minlat=minlat,maxlat=maxlat,ddeg=ddeg)
            self.distr0=D
            m=np.max(D.V1d)
            oc = np.array([1. if D.V1d[i]==m else 0. for i in range(len(D.V1d))])
            self.Ocean=oceanvector(oc,minlon=minlon,maxlon=maxlon,minlat=minlat,maxlat=maxlat,ddeg=ddeg)
            co = np.array([1. if (D.V1d[i]!=0 and self.Ocean.V1d[i]==0.) else 0. for i in range(len(D.V1d))])
            self.Coast=oceanvector(co,minlon=minlon,maxlon=maxlon,minlat=minlat,maxlat=maxlat,ddeg=ddeg)

    def avg_distribution(self):
        avg_distr = np.zeros((len(Lats),len(Lons)))
        for k in range(len(self.lon[0])):
            distrk, _, _ = np.histogram2d(self.lat[:,k], self.lon[:,k], [lat_bins, lon_bins])
            avg_distr+=distrk
        avg_distr/=len(self.lon[0])
        avg_distr=oceanvector(avg_distr,minlon=minlon,maxlon=maxlon,minlat=minlat,maxlat=maxlat,ddeg=ddeg)
        return avg_distr

    def get_beached(self):
        lonbeached = []
        latbeached=[]
        for i in range(len(self.lon)):
             if np.ma.is_masked(self.lon[i]):
                 ind=np.where(np.ma.getmask(self.lon[i]))[0][0]-1
                 lonbeached.append(self.lon[i,ind])
                 latbeached.append(self.lat[i,ind])
        
        return (lonbeached, latbeached)
    
    @classmethod    
    def from_coastal_particles(cls,FlowType, tmin=-2, tmax=-1):        
        f = TRAJfiles[FlowType]
        datadir =TRAJpath + FlowType + '/'
        i = 0
        pfile = datadir + f + str(i) + '.nc'
        data = Dataset(pfile, 'r')
        lon0=data.variables['lon'][:,0]
        lat0=data.variables['lat'][:,0]
        
        for i in range(1,N):
            pfile = datadir + f + str(i) + '.nc'
            data = Dataset(pfile, 'r')
            lon0=np.append(lon0, data.variables['lon'][:,0])
            lat0=np.append(lat0, data.variables['lat'][:,0])
        
        distr , _, _ = np.histogram2d(lat0, lon0, [lat_bins, lon_bins])
        D=oceanvector(distr)
        m=np.max(D.V1d)
        Ocean=np.array([1. if D.V1d[i]==m else 0. for i in range(len(D.V1d))])
        Coast=np.array([1. if (D.V1d[i]!=0 and Ocean[i]==0.) else 0. for i in range(len(D.V1d))])
        bindex = np.array([int(((la-np.min(Lats))//ddeg)*len(Lons)+(lo-np.min(Lons))//ddeg) for lo,la in zip(lon0,lat0)])
        oncoast = np.array([i if Coast[bindex[i]] ==1. else 0 for i in range(len(bindex))])
        oncoast = oncoast[oncoast!=0]
        return cls(FlowType,tmin=tmin,tmax=tmax,indices=oncoast)
    
    @classmethod
    def from_westcoast_particles(cls,FlowType, tmin=-2, tmax=-1):
        f = TRAJfiles[FlowType]
        datadir =TRAJpath + FlowType + '/'
        i = 0
        pfile = datadir + f + str(i) + '.nc'
        data = Dataset(pfile, 'r')
        lon0=data.variables['lon'][:,0]
        lat0=data.variables['lat'][:,0]        
        for i in range(1,N):
            pfile = datadir + f + str(i) + '.nc'
            data = Dataset(pfile, 'r')
            lon0=np.append(lon0, data.variables['lon'][:,0])
            lat0=np.append(lat0, data.variables['lat'][:,0])
        
        distr , _, _ = np.histogram2d(lat0, lon0, [lat_bins, lon_bins])
        D=oceanvector(distr)
        m=np.max(D.V1d)
        Ocean=np.array([1. if D.V1d[i]==m else 0. for i in range(len(D.V1d))])
        Coast=np.array([1. if (D.V1d[i]!=0 and Ocean[i]==0.) else 0. for i in range(len(D.V1d))])        
        bindex = np.array([int(((la-np.min(Lats))//ddeg)*len(Lons)+(lo-np.min(Lons))//ddeg) for lo,la in zip(lon0,lat0)])
        oncoast = np.array([i if Coast[bindex[i]] ==1. and lon0[i]<-40 else 0 for i in range(len(bindex))])
        oncoast = oncoast[oncoast!=0]
        return cls(FlowType,tmin=tmin,tmax=tmax,indices=oncoast)    
    
    def movie_me(self,FPS=15,outname='nonamemovie',land=True):
        fig = plt.figure(figsize=(15, 10))
        Time=self.Time
        
        def settitle(time):
            titstr = (self.time_origin + datetime.timedelta(days=time)).strftime("%Y-%m-%d")
            return titstr
        
        ax = fig.add_subplot(1, 1, 1)
        ttl = ax.set_title(settitle(int(Time[0])), size=22)
        m = Basemap(projection='mill', llcrnrlat=-75., urcrnrlat=75., llcrnrlon=0., urcrnrlon=360., resolution='l')
        m.drawcoastlines()
        if land:
            m.fillcontinents(color='dimgrey')
        m.drawmapboundary(fill_color='lightgray')
        m.drawparallels(np.array([-60,-40,-20,0,20,40,60]), labels=[True, False, False, True])
        m.drawmeridians(np.array([10,50,100,150,200,250,300,350]), labels=[False, False, False, True])
        xs, ys = m(self.lon[:,0], self.lat[:,0])
        scat = m.scatter(xs, ys, marker='.', s=1)
        
        def animate(t):
            if t%100 ==0:
                print 'time: ', t
            scat.set_offsets(np.matrix(m(self.lon[:,t], self.lat[:,t])).transpose())
            ttl.set_text(settitle(int(Time[t])))
            return scat
    
        anim = animation.FuncAnimation(fig, animate, frames=range(0,len(Time)-2),blit=False)
        anim.save(outname + '.mp4', fps=FPS, extra_args=['-vcodec', 'libx264'])      

    def Loss_regions(self):
        #See where the particles that got deleted start for trajectory files (start and end used)
        keep = np.invert(np.isnan(self.lon[:,-1]))
        distr , _, _ = np.histogram2d(self.lat[keep,0], self.lon[keep,0], [lat_bins, lon_bins])
        distr_i , _, _ = np.histogram2d(self.lat[:,0], self.lon[:,0], [lat_bins, lon_bins])    
        Loss = oceanvector(distr/distr_i*100)
        return Loss

    def get_connection_times(self,Lonmin, Lonmax, Latmin, Latmax):
        PosTime = []
        for i in range(len(self.lon)):
            if i%1000==0:
                print i
            for t in range(len(self.lon[0,:])):
                if (self.lon[i,t]>Lonmin and self.lon[i,t]<Lonmax and self.lat[i,t]>Latmin and self.lat[i,t]<Latmax):
                    lo=self.lon[i,0]
                    la=self.lat[i,0]
                    bindex = int(((la-np.min(Lats))//ddeg)*len(Lons)+(lo-np.min(Lons))//ddeg)
                    PosTime.append([bindex,self.time[i,t]])
                    break
        PosTime = np.array(PosTime)
        print 'pos: ', PosTime.shape
        BinTime = np.zeros(len(Lons)*len(Lats))        
        for k in range(len(BinTime)):
            C=PosTime[np.where(PosTime[:,0]==k)]
            print 'Cshape: ', C.shape
            if not C.shape[0]==0:
                m=np.nanmin(C,axis=0)
                BinTime[k]=m[1]
            else:
                BinTime[k]=np.nan
            
        B=oceanvector(BinTime)
        return B

    def Particle_Locations(self):
        distr, _, _ = np.histogram2d(self.lat[:,-1], self.lon[:,-1], [lat_bins, lon_bins])
        D = oceanvector(distr)
        N_coast = np.dot(D.V1d,self.Coast.V1d)
        N_ocean = np.sum(distr) - N_coast
        N_lost=np.sum(self.distr0.V1d)-N_coast-N_ocean        
        return (N_coast, N_ocean, N_lost)



def movie():
    for n, d in files.iteritems():
        print 'moviefile: ', n
        pdir = datadir #+ n + '/'
        tmin=0
        tmax=-1
        P = ParticleData(pdir,d,tmin=tmin,tmax=tmax,load_initial=False)
        P.movie_me(outname=plotdir + 'Movies/FullSim_' + n,land=False)

movie()
    
def SimpleStatsParticles():
    tmin=0
    tmax=-1
        
    plt.figure()
    for n, d in files.iteritems():
        P = ParticleData(datadir, d ,tmin=tmin,tmax=tmax)        
        lonmean = [np.nanmean(P.lon[:,t]) for t in range(len(P.lon[0]))]
        plt.plot(range(len(P.lon[0])),lonmean,label=n)
    plt.legend()
    plt.show()
        
#SimpleStatsParticles()
        
def plot_final_avg_distr():
    """
    Average distribution over the defined time
    """
    tmin=1700 #For 1 year
    tmax=1820
 
    for n, d in files.iteritems():
        pdir = datadir
        P = ParticleData(pdir,d,tmin=tmin,tmax=tmax)
        D=P.avg_distribution()
        P.distr0.plot_me(outname = plotdir + 'Distr0_'+n,save=True)
        Ocean = P.Ocean
        Docean=oceanvector(np.multiply(D.V1d,Ocean.V1d),minlon=minlon,maxlon=maxlon,minlat=minlat,maxlat=maxlat,ddeg=ddeg)
        Dnorm=Docean.normalize_area()
        Dnorm.plot_me(aaoutname = plotdir + 'Distr_'+n,cbartitle=r'#/$km^2$',title = n,land=False,save=True,pl='contour')

#plot_final_avg_distr()


def deleted_particles():
    tmin=0
    tmax=-1

    nsurv=[]
#    for n, d in files.iteritems():
    n='P'
    d = files[n]

    P = ParticleData(datadir, d ,tmin=tmin,tmax=tmax)
    for t in range(len(P.time[0])):
        lo=P.lon[:,t]
        lons = lo[~lo.mask]
        nsurv.append(len(lons))
    plt.plot(nsurv)        
    
#deleted_particles()
    
def beached_particles_map():
    tmin=0
    tmax=-1

#    for n, d in files.iteritems():
    n='S'
    d = files[n]
    P = ParticleData(datadir, d ,tmin=tmin,tmax=tmax)
    (lons,lats)=P.get_beached()
    
    plt.figure(figsize=(12,8))
    m = Basemap(projection='mill', llcrnrlat=-80., urcrnrlat=80., llcrnrlon=0., urcrnrlon=360., resolution='l')
    m.drawcoastlines()
    m.drawmapboundary(fill_color='lightgray')
    m.drawparallels(np.array([-60,-40,-20,0,20,40,60]), labels=[True, False, False, True])
    m.drawmeridians(np.array([10,50,100,150,200,250,300,350]), labels=[False, False, False, True])
    xs, ys = m(lons, lats)
    m.scatter(xs, ys, marker='.') #, s=1)
    
#beached_particles_map()

def Particles_over_time():
    """
    Bar-plot for particles in ocean, coasts and lost for Traj
    """
    for n, d in files.iteritems():
        Trange=np.linspace(0,1815,16)
        N_Coast =[]
        N_Ocean =[]
        N_Lost =[]
        
        for t in Trange[0:-1]:
            print 'time: ', t
            tmin=int(t)
            tmax=tmin+1
            P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
            (c,o,l)=P.Particle_Locations()
            N_Coast.append(c)
            N_Ocean.append(o)
            N_Lost.append(l)
    
        N0=N_Coast[0] + N_Ocean[0] + N_Lost[0]
        N_Ocean=np.array(N_Ocean)/N0*100.
        N_Coast=np.array(N_Coast)/N0*100.
        N_Lost=100-N_Ocean-N_Coast
       
        Time = range(15)    
        plt.figure(figsize=(12,9))
        p1=plt.bar(Time, N_Ocean)
        p2=plt.bar(Time, N_Coast,bottom=N_Ocean)
        p3=plt.bar(Time, N_Lost,bottom=N_Coast+N_Ocean)    
        plt.legend((p1[0], p2[0],p3[0]), ('Ocean', 'Coast','Lost'),prop={'size': 18},loc=3)
        plt.title('Particle Locations',size=18)
        plt.xlabel('Years',size=12)
        plt.ylabel('% of particles',size=12)
        plt.savefig(plotdir + 'ParticlesOverTime/' + FlowType,bbox_inches='tight' )



    

def SimpleStats():
    conf=1.
    tmin=0
    tmax=-1        
    for n, d in files.iteritems():
        P = ParticleData(datadir + d ,tmin=tmin,tmax=tmax)
        
        D=P.avg_distribution()
        Docean=oceanvector(np.multiply(D.V1d,P.Ocean.V1d)).normalize_area()

        DV2=Docean.V2d
        LonDistr=np.sum(DV2,axis=0)/np.sum(DV2)
        LatDistr=np.sum(DV2,axis=1)/np.sum(DV2)
        LonMean = np.dot(Lons,LonDistr)
        LatMean = np.dot(Lats,LatDistr)    
        LonStd=np.sqrt(np.dot(LonDistr,(Lons-LonMean) ** 2))
        LatStd=np.sqrt(np.dot(LatDistr,(Lats-LatMean) ** 2))    
    
        fig = plt.figure(figsize=(15,6))
        ax = fig.add_subplot(121)
        ax.grid(True)
        ax.plot(Lons,LonDistr,'k')
        ax.yaxis.set_ticks(np.arange(0.0, 0.08, 0.02))
        ax.xaxis.set_ticks(np.arange(-100, 30, 20.))
        ax.set_xlim([-100,30])
        ax.set_ylim([0,0.08])
        ax.set_xlabel('Longitude [deg]',size=14)
        ax.set_title('Zonal distribution',size=15)
        ax.axvline(x=LonMean,label='Mean')
        plt.axvline(x=LonMean+conf*LonStd,c='g',ls='--',label=r'Mean $\pm$ Std')    
        plt.axvline(x=LonMean-conf*LonStd,c='g',ls='--')
        ax.legend()
        
        ax = fig.add_subplot(122)
        ax.grid(True)
        ax.plot(Lats,LatDistr,'k')
        ax.yaxis.set_ticks(np.arange(0.0, 0.18, 0.04))
        ax.xaxis.set_ticks(np.arange(0, 75, 10.))
        ax.set_xlim([0,75])
        ax.set_ylim([0,0.18])
        ax.set_xlabel('Latitude [deg]',size=14)
        ax.set_title('Meridional distribution',size=15)
        ax.axvline(x=LatMean,label='Mean')
        plt.axvline(x=LatMean+conf*LatStd,c='g',ls='--',label=r'Mean $\pm$ Std')    
        plt.axvline(x=LatMean-conf*LatStd,c='g',ls='--')
        ax.legend()
        
        #Save Mean and Std
        np.savez(plotdir + 'GP_connection/' + FlowType + '_stats', LonMean=LonMean, LonStd=LonStd, LatMean=LatMean, LatStd=LatStd)
        fig.savefig(plotdir + 'SimpleStats/'+ FlowType + '.eps',bbox_inches='tight',dpi=100)

    
    


def GP_connection():
    """
    Function that plots the connection time for each grid cell
    """
    #    for FlowType in TRAJdirs:
    print 'Running GP connection'
    
    FlowType = 'Traj_PN'
    statfile = np.load(plotdir + 'GP_connection/Traj_PN_stats' + '.npz')
    LonMean = statfile['LonMean']
    LatMean = statfile['LatMean']
    LonStd = statfile['LonStd']   
    LatStd = statfile['LatStd']
    tmin=0
    tmax=-1
    Lonmin=LonMean-0.5*LonStd
    Lonmax=LonMean+0.5*LonStd
    Latmin=LatMean-0.5*LatStd    
    Latmax=LatMean+0.5*LatStd
    P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
    B=P.get_connection_times(Lonmin,Lonmax,Latmin,Latmax)
    levels = np.arange(0, 1000, 50)
    B.save_me(plotdir + 'GP_connection/'+FlowType + 'ConnTime')
    B.plot_me(outname=plotdir+'GP_connection/V4'+FlowType,pl='contour',cbartitle='days',title='Time needed to reach garbatch patch',land=True,save=True,levels=levels)

#GP_connection()


#def share_in_gp():
#    tmin=-10
#    tmax=-1
#    FlowType='Traj_PN'
#    P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
    


def fit_data():
    import scipy.io as sio
    from scipy import interpolate
    data = sio.loadmat('/Users/wichmann/paper_v3/Observations/Observations.mat')
    npcell = np.array(data['data_npkm2'])
    lon = np.array(data['lon'])
    lat = np.array(data['lat'])
    lon = lon[:,0]
    lat = lat[:,0]
    npcell = npcell[:,0]
    
    #Cut out North Atlantic. Note that we use thate data between 0 and 30 degrees is only in the mediterranian (which we are not interested in)    
    northatlantic=[False if lo<260. or la < 0. or la >75. or (lo<280 and la<11) or (lo<270 and la<20) else True for lo,la in zip(lon,lat)]
    npcell=npcell[northatlantic]
    lon=lon[northatlantic]
    lat=lat[northatlantic]    
    lon=np.array([lo-360. if lo>180. else lo for lo in lon])

    assert (len(lon)==len(lat))
    assert (len(lon)==len(npcell))

    tmin=400
    tmax=500
    
    for FlowType in [TRAJdirs[0]]:
        #Compute average over different average distributions for the 14 years        
        D = ParticleData(FlowType,tmin=tmin,tmax=tmax)
        Davg=D.avg_distribution()
#        Docean = oceanvector(np.multiply(D.Ocean.V1d,Davg.V1d))
        Distr = Davg.normalize_area()
        Lo = Davg.Lons
        La = Davg.Lats
        print 'plankv1d: '
        Plank = D.avg_plankton()
        Plank=oceanvector(np.log(Plank.V1d))
        Plank.plot_me()
        
        distr_interp = interpolate.interp2d(Lo, La, Distr.V2d, kind='linear')
        print 'test: ', distr_interp(-40,30)
#        
#
#        plank_interp = interpolate.interp2d(Lo, La, Plank.V1d, kind='linear')
#        
#        model_distr = np.array([distr_interp(lon[i],lat[i])[0] for i in range(len(lon))])
#        model_plankton = np.array([plank_interp(lon[i],lat[i])[0] for i in range(len(lon))])
#        print 'model plank: '
#        
#        for i in range(len(model_plankton)):
#            print model_plankton[i]
#        
#        from sklearn import linear_model
#        from sklearn.preprocessing import Imputer
#       
#        inds = [True if (model_plankton[i] is not np.nan and model_plankton[i]!=0 and model_distr[i] is not np.nan and model_distr[i]!=0.)  else False for i in range(len(model_plankton))]
#       
#        y=np.log(npcell[inds])
#        distr = np.log(model_distr[inds])
#        plank = np.log(model_plankton[inds])
#        print distr[distr==float('nan')]
#        print plank[plank==float('nan')]
#        
#        print np.max(y)
#        print np.max(distr)
#        print plank
#        
#        X=[[distr[i],plank[i]] for i in range(len(plank))]
#    
#        reg = linear_model.LinearRegression()
#        
#        reg.fit (X, y)
#        print reg.coef_
#        
#        
#        
#        diff = model_data - npcell          
#        
#        assert(np.sum(diff)==0)
#        
#        fig = plt.figure(figsize=(10,15))
#        ax1 = fig.add_subplot(211)
#        ax1.grid(True)
#        ax1.set_yscale('log')
#        ax1.set_xscale('log')
#        ax1.scatter(npcell,model_data,c=diff)
#        ax1.set_xlabel('Data',size=14)
#        ax1.set_ylabel(TRAJnames[FlowType],size=14)
#        ax1.set_title('Comparison',size=15)
#
#        ax2 = fig.add_subplot(212)
#        m = Basemap(projection='merc', llcrnrlat=0, urcrnrlat=75., llcrnrlon=-100, urcrnrlon=30, resolution='l')
#        m.drawcoastlines()
#        m.drawparallels(np.arange(np.min(Lats),np.max(Lats),10.),labels=[True,False,True,False])
#        m.drawmeridians(np.arange(np.min(Lons),np.max(Lons),20.),labels=[False,False,False,True])
#        m.drawmapboundary(fill_color='lightgray')
#        xs,ys = m(lon,lat)
#        sc=m.scatter(xs,ys,c=diff)
#        fig.colorbar(sc, ax=ax2, format='%.0e', shrink=0.8)
#        ax2.set_title('Model - Data',size=15)
#
#        plt.savefig(plotdir + 'ComparisonData/FllTimeAvg' + FlowType,bbox_inches='tight',dpi=100)
#        plt.close()    
    
    
#fit_data()   








def emtpy_NA():
    plt.figure(figsize=(19.2,10.80))
    m = Basemap(projection='merc', llcrnrlat=np.min(Lats), urcrnrlat=np.max(Lats), llcrnrlon=np.min(Lons), urcrnrlon=np.max(Lons), resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(np.min(Lats),np.max(Lats),10.),labels=[True,False,True,False])
    m.drawmeridians(np.arange(np.min(Lons),np.max(Lons),20.),labels=[False,False,False,True])
    m.fillcontinents(color='dimgrey')
    m.drawmapboundary(fill_color='lightgray')
    plt.savefig(plotdir + 'NorthAtlantic',bbox_inches='tight',dpi=100)

#emtpy_NA()

def plot_data():
    """
    Simple plot of data
    """
    import scipy.io as sio
    data = sio.loadmat('/Users/wichmann/paper_v3/Observations/Observations.mat')
    
    npcell = np.array(data['data_npkm2'])
    lon = np.array(data['lon'])
    lat = np.array(data['lat'])
    lon = lon[:,0]
    lat = lat[:,0]
    npcell = npcell[:,0]
    
    #Cut out North Atlantic. Note that we use thate data between 0 and 30 degrees is only in the mediterranian (which we are not interested in)    
    northatlantic=[False if lo<260. or la < 0. or la >75. or (lo<280 and la<11) or (lo<270 and la<20) else True for lo,la in zip(lon,lat)]
    npcell=npcell[northatlantic]
    lon=lon[northatlantic]
    lat=lat[northatlantic]    
    lon=np.array([lo-360. if lo>180. else lo for lo in lon])

    assert (len(lon)==len(lat))
    assert (len(lon)==len(npcell))
    
    plt.figure(figsize=(19.2,10.80))
    m = Basemap(projection='merc', llcrnrlat=0, urcrnrlat=75., llcrnrlon=-100, urcrnrlon=30, resolution='l')
    m.drawcoastlines()
    m.drawparallels(np.arange(0,75,10.),labels=[True,False,True,False])
    m.drawmeridians(np.arange(260,390,20.),labels=[False,False,False,True])
    m.drawmapboundary(fill_color='lightgray')
    m.fillcontinents(color='dimgrey')
    xs,ys = m(lon,lat)
    plt.title('Plastic data',size=18)   
    m.scatter(xs, ys, c=npcell, cmap=cm2.coolwarm)
    cbar = plt.colorbar(orientation='vertical',shrink=0.8,format='%.1e')
    cbar.set_label(r'#/$km^2$',size=15)
    plt.savefig(plotdir + 'Datakm2',bbox_inches='tight',dpi=100)
    
#plot_data()

def LonLat_scatter():
    minlat = 25
    maxlat= 35
    minlon = -50
    maxlon=-30
    tmin=-121
    tmax=-1
    
    #Draw square
    Region = oceanvector(np.array([1 if (Lats[i]<=maxlat and Lats[i]>=minlat) or (Lons[j]<=maxlon and Lons[j]>=minlon) else np.nan for i in range(len(Lats)) for j in range(len(Lons))]))
    Region.plot_me(land=True,colbar=False)
    
    import scipy.io as sio
    data = sio.loadmat('/Users/wichmann/paper_v3/Observations/Observations.mat')

    N = np.array(data['data_npkm2'])
    lon = np.array(data['lon'])
    lat = np.array(data['lat'])
    lon = lon[:,0]
    lat = lat[:,0]
    N = N[:,0]
    
    northatlantic=[False if lo<260. or la < 0. or la >75. or (lo<280 and la<11) or (lo<270 and la<20) else True for lo,la in zip(lon,lat)]
    N=N[northatlantic]
    lon=lon[northatlantic]
    lat=lat[northatlantic]    
    lon=np.array([lo-360. if lo>180. else lo for lo in lon])
    
    #Take only the ones in the correct range
    #Longitudal plot
    LonSelect = [True if (lat[i]>=minlat and lat[i]<=maxlat) else False for i in range(len(lat))]
    Lonlon=lon[LonSelect]
    Lonlat=lat[LonSelect]
    LonN=N[LonSelect]    

    #Latitudal plot
    LatSelect = [True if (lon[i]>=minlon and lon[i]<=maxlon) else False for i in range(len(lat))]
    Latlon=lon[LatSelect]
    Latlat=lat[LatSelect]
    LatN=N[LatSelect]  
    
    #Get the particles for Lon plot
    LonParticles=[]
    for i in range(len(Lonlon)):         
        bindex = int(((Lonlat[i]-np.min(Lats))//ddeg)*len(Lons)+(Lonlon[i]-np.min(Lons))//ddeg)
        LonParticles.append([bindex,LonN[i]])        
    LonParticles = np.array(LonParticles)

    #Get the particles for Lat plot
    LatParticles=[]
    for i in range(len(Latlon)):         
        bindex = int(((Latlat[i]-np.min(Lats))//ddeg)*len(Lons)+(Latlon[i]-np.min(Lons))//ddeg)
        LatParticles.append([bindex,LatN[i]])        
    LatParticles = np.array(LatParticles)

    PsLon = np.zeros(len(Lons)*len(Lats))    
    for k in range(len(PsLon)):
        C=LonParticles[np.where(LonParticles[:,0]==k)]
        mean = np.nanmean(C)
        PsLon[k]=mean
    
    PsLat = np.zeros(len(Lons)*len(Lats))    
    for k in range(len(PsLat)):
        C=LatParticles[np.where(LatParticles[:,0]==k)]
        mean = np.nanmean(C)
        PsLat[k]=mean
    
    DLon=oceanvector(PsLon)
    DLat=oceanvector(PsLat)
    DLonV2=DLon.V2d
    DLatV2=DLat.V2d

    LonDistrData=np.nanmean(DLonV2,axis=0)
    LatDistrData=np.nanmean(DLatV2,axis=1)
    
    LonRegion = np.array([1 if (Lats[i]<=maxlat and Lats[i]>=minlat) else 0 for i in range(len(Lats)) for _ in range(len(Lons))])
    LatRegion = np.array([1 if (Lons[i]<=maxlon and Lats[i]>=minlon) else 0 for i in range(len(Lats)) for _ in range(len(Lons))])
    RLon=oceanvector(LonRegion)    
    RLat=oceanvector(LatRegion)    
    
    for FlowType in TRAJdirs:
        P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
        D=P.avg_distribution()
        Docean=oceanvector(np.multiply(D.V1d,P.Ocean.V1d)).normalize_area()
        
        DLon = oceanvector(np.multiply(RLon.V1d,Docean.V1d))
        DLat = oceanvector(np.multiply(RLat.V1d,Docean.V1d))
        
        LonDistrSimulation=np.nanmean(DLon.V2d,axis=0)
        LatDistrSimulation=np.nanmean(DLat.V2d,axis=1)
        
        Lons_interpolate = interpolate.interp1d(Lons, LonDistrSimulation, kind='linear')
        Lats_interpolate = interpolate.interp1d(Lats, LatDistrSimulation, kind='linear')

        LonSim_points = np.array([Lons_interpolate(Lons[i]) if Lons[i] is not np.nan else np.nan for i in range(len(Lons))])
        LatSim_points = np.array([Lats_interpolate(Lats[i]) if Lats[i] is not np.nan else np.nan for i in range(len(Lats))])
        
        LonSim_points=LonSim_points/np.nansum(LonSim_points)*np.nansum(LonDistrData)
        LatSim_points=LatSim_points/np.nansum(LatSim_points)*np.nansum(LatDistrData)
        
        assert(abs(np.nansum(LonSim_points)-np.nansum(LonDistrData))<0.01)
        assert(abs(np.nansum(LatSim_points)-np.nansum(LatDistrData))<0.01)

        
        fig = plt.figure(figsize=(15,6))
        ax = fig.add_subplot(121)
        ax.grid(True)
        ax.plot(Lons,LonSim_points,c='k',label='Simulation')
        ax.scatter(Lons,LonDistrData, c='r',label='Data')
        ax.set_xlabel('Longitude [deg]',size=14)
        ax.set_ylabel(r'Particles per $km^2$',size=14)
        ax.set_title(r'Zonal distribution, Lat = $[%s^\circ,%s^\circ]$'%(minlat,maxlat),size=15)
        ax.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
        ax.legend()
        
        ax = fig.add_subplot(122)
        ax.grid(True)
        ax.plot(Lats,LatSim_points,c='k',label='Simulation')
        ax.scatter(Lats,LatDistrData, c='r',label='Data')
        ax.set_xlabel('Latitude [deg]',size=14)
        ax.set_ylabel(r'Particles per $km^2$',size=14)
        ax.set_title(r'Meridional distribution, Lon = $[%s^\circ,%s^\circ]$'%(minlon,maxlon),size=15)
        ax.ticklabel_format(style='sci',axis='y', scilimits=(0,0))
        ax.legend()
        plt.show()
    
        fig.savefig(plotdir + 'LonLatDistributions/'+ FlowType,bbox_inches='tight',dpi=100)
    
#LonLat_scatter()
    



def test_me():
    FlowType='Traj_PN'
    tmin=0
    tmax=1
    P = ParticleData.from_westcoast_particles(FlowType,tmin=tmin,tmax=tmax)  
    D=P.avg_distribution()
    D.plot_me()
    
    
 




def loss_regions():
    """
    Share of particles that is lost from each grid cell
    """    
    for FlowType in TRAJdirs:
        tmin=0 #For 1 year
        tmax=121
        P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
        D=P.Loss_regions()
        D.plot_me(outname=plotdir+'LossRegions/'+FlowType,cbartitle='%',title='Share of lost particles after 1 year',land=True,save=True)


def plot_final_avg_distr():
    """
    Average distribution over the defined time
    """
    for FlowType in TRAJdirs:
        n=TRAJnames[FlowType]
        tmin=-121 #For 1 year
        tmax=-1
        P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
        D=P.avg_distribution()
        Ocean = P.Ocean
        Docean=oceanvector(np.multiply(D.V1d,Ocean.V1d))
        Dnorm=Docean.normalize_area()
        Dnorm.plot_me(outname = plotdir + 'AvgDistributions/'+FlowType,cbartitle=r'#/$km^2$',title = n,land=True,save=True,pl='contour')

def SimpleStats():
    conf=1.
    tmin=-121
    tmax=-1
        
    for FlowType in TRAJdirs:
        P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
        D=P.avg_distribution()
        Docean=oceanvector(np.multiply(D.V1d,P.Ocean.V1d)).normalize_area()

        DV2=Docean.V2d
        LonDistr=np.sum(DV2,axis=0)/np.sum(DV2)
        LatDistr=np.sum(DV2,axis=1)/np.sum(DV2)
        LonMean = np.dot(Lons,LonDistr)
        LatMean = np.dot(Lats,LatDistr)    
        LonStd=np.sqrt(np.dot(LonDistr,(Lons-LonMean) ** 2))
        LatStd=np.sqrt(np.dot(LatDistr,(Lats-LatMean) ** 2))    
    
        fig = plt.figure(figsize=(15,6))
        ax = fig.add_subplot(121)
        ax.grid(True)
        ax.plot(Lons,LonDistr,'k')
        ax.yaxis.set_ticks(np.arange(0.0, 0.08, 0.02))
        ax.xaxis.set_ticks(np.arange(-100, 30, 20.))
        ax.set_xlim([-100,30])
        ax.set_ylim([0,0.08])
        ax.set_xlabel('Longitude [deg]',size=14)
        ax.set_title('Zonal distribution',size=15)
        ax.axvline(x=LonMean,label='Mean')
        plt.axvline(x=LonMean+conf*LonStd,c='g',ls='--',label=r'Mean $\pm$ Std')    
        plt.axvline(x=LonMean-conf*LonStd,c='g',ls='--')
        ax.legend()
        
        ax = fig.add_subplot(122)
        ax.grid(True)
        ax.plot(Lats,LatDistr,'k')
        ax.yaxis.set_ticks(np.arange(0.0, 0.18, 0.04))
        ax.xaxis.set_ticks(np.arange(0, 75, 10.))
        ax.set_xlim([0,75])
        ax.set_ylim([0,0.18])
        ax.set_xlabel('Latitude [deg]',size=14)
        ax.set_title('Meridional distribution',size=15)
        ax.axvline(x=LatMean,label='Mean')
        plt.axvline(x=LatMean+conf*LatStd,c='g',ls='--',label=r'Mean $\pm$ Std')    
        plt.axvline(x=LatMean-conf*LatStd,c='g',ls='--')
        ax.legend()
        
        #Save Mean and Std
        np.savez(plotdir + 'GP_connection/' + FlowType + '_stats', LonMean=LonMean, LonStd=LonStd, LatMean=LatMean, LatStd=LatStd)
        fig.savefig(plotdir + 'SimpleStats/'+ FlowType + '.eps',bbox_inches='tight',dpi=100)

        
def plankton_sinking():
    """
    Function to plot total number of particles for different Plankton parameters s
    """
    f=[]
    srange=np.logspace(0.,1.5,30)
    FlowType='Traj_PN'
    tmin=-121
    tmax=-1
    P = ParticleData(FlowType,tmin=tmin,tmax=tmax)

    for s in srange:
        print 's: ', s        
        D_npp=P.avg_npp_int(s)
        DOcean=np.multiply(D_npp.V1d,P.Ocean.V1d)
        f.append(np.sum(DOcean))
    plt.plot(srange,f)
    

def difference_plots():
    """
    Difference plots for different distributions
    """
    tmin=-121
    tmax=-1
    P_passive = ParticleData(TRAJdirs[0],tmin=tmin,tmax=tmax)
    P_inertial = ParticleData(TRAJdirs[1],tmin=tmin,tmax=tmax)    
    D_passive=P_passive.avg_distribution()
    D_inertial=P_inertial.avg_distribution()
    Diff = oceanvector(D_passive.V1d-D_inertial.V1d)
    Diff.plot_me()


def mean_over_time():
    """
    Variation of the mean over time
    """
    FlowType=TRAJdirs[0]
    LonMeans =[]
    LatMeans =[]    
    Trange=range(1820)
    
    for t in Trange:
        tmin=int(t)
        tmax=tmin+1        
        P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
        D=P.avg_distribution()
        Docean=oceanvector(np.multiply(D.V1d,P.Ocean.V1d)).normalize_area()
        DV2=Docean.V2d
        LonDistr=np.sum(DV2,axis=0)/np.sum(DV2)
        LatDistr=np.sum(DV2,axis=1)/np.sum(DV2)
        LonMean = np.dot(Lons,LonDistr)
        LatMean = np.dot(Lats,LatDistr)    
        LonMeans.append(LonMean)
        LatMeans.append(LatMean)
    np.savez(plotdir + 'MeanOverTime/' + FlowType, LonMeans = LonMeans, LatMeans = LatMeans)

def Plot_mean_over_time():
    """
    Plot data from previous function
    """
    FlowType = 'Traj_PN'
    statfile = np.load(plotdir + 'MeanOverTime/' + FlowType + '.npz')
    LonMeans = statfile['LonMeans']
    LatMeans = statfile['LatMeans']
    
    plt.plot(LonMeans)
    plt.show()
    plt.plot(LatMeans)

def Compare_distributions():
    """
    Scatter plot for different distributions
    """
    tmin=-121
    tmax=-1
    P_passive = ParticleData('Traj_PN',tmin=tmin,tmax=tmax)    
    D_passive=P_passive.avg_distribution()
    Docean_passive=oceanvector(np.multiply(D_passive.V1d,P_passive.Ocean.V1d))
    Docean_passive=Docean_passive.normalize()
    
    print 'Sum (should be 1.0): ', np.sum(Docean_passive.V1d)

    for FlowType in TRAJdirs:
        print 'Compare to: ', FlowType
        P = ParticleData(FlowType,tmin=tmin,tmax=tmax)
        D=P.avg_distribution()
        Docean=oceanvector(np.multiply(D.V1d,P_passive.Ocean.V1d))
        Docean=Docean.normalize()      
        print 'Sum (should be 1.0) comp: ', np.sum(Docean.V1d)
        plt.figure(figsize=(10.8,10.80))
        plt.scatter(Docean_passive.V1d,Docean.V1d,c='r',s=3,label='Simulations')
        plt.plot(Docean_passive.V1d,Docean_passive.V1d,c='k',label='y=x')
        plt.xlabel(TRAJnames['Traj_PN'],size=14)
        plt.ylabel(TRAJnames[FlowType],size=14)
        plt.xlim(xmin=0)
        plt.ylim(ymin=0)        
        plt.title('Year-averaged final distributions in the ocean',size=16)
        plt.legend(loc=4,prop={'size': 14 })
        plt.savefig(plotdir + 'DistributionComparison/' + FlowType,bbox_inches='tight',dpi=100)
        plt.close()


def Compare_to_data():
    import scipy.io as sio
    from scipy import interpolate
    data = sio.loadmat('/Users/wichmann/paper_v3/Observations/Observations.mat')
    npcell = np.array(data['data_npcell'])
    lon = np.array(data['lon'])
    lat = np.array(data['lat'])
    lon = lon[:,0]
    lat = lat[:,0]
    npcell = npcell[:,0]
    
    #Cut out North Atlantic. Note that we use thate data between 0 and 30 degrees is only in the mediterranian (which we are not interested in)    
    northatlantic=[False if lo<260. or la < 0. or la >75. or (lo<280 and la<11) or (lo<270 and la<20) else True for lo,la in zip(lon,lat)]
    npcell=npcell[northatlantic]
    lon=lon[northatlantic]
    lat=lat[northatlantic]    
    lon=np.array([lo-360. if lo>180. else lo for lo in lon])

    assert (len(lon)==len(lat))
    assert (len(lon)==len(npcell))
    
    for FlowType in TRAJdirs:
        
        #Compute average over different average distributions for the 14 years        
        P=np.empty(15,dtype=oceanvector) 
        for i in range(15):
            print 'time: ', i
            tmin=i*121
            tmax=tmin+121
            D = ParticleData(FlowType,tmin=tmin,tmax=tmax)
            Davg=D.avg_distribution()
            Docean = np.multiply(D.Ocean.V1d,Davg.V1d)
            P[i]=oceanvector(Docean)
        Davg = np.zeros(len(P[0].V1d))        
        for i in range(15):
            Davg+=P[i].V1d
            
        Davg/=len(P)
        Docean = oceanvector(Davg)
        Distr = Docean.normalize()
       
        Lo = Docean.Lons
        La = Docean.Lats
        distr_interp = interpolate.interp2d(Lo, La, Distr.V1d, kind='linear')
        model_data = np.array([distr_interp(lon[i],lat[i])[0] for i in range(len(lon))])        
        model_data = model_data/np.sum(model_data) * np.sum(npcell)
        diff = model_data - npcell          
        
        assert(np.sum(diff)==0)
        
        fig = plt.figure(figsize=(10,15))
        ax1 = fig.add_subplot(211)
        ax1.grid(True)
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.scatter(npcell,model_data,c=diff)
        ax1.set_xlabel('Data',size=14)
        ax1.set_ylabel(TRAJnames[FlowType],size=14)
        ax1.set_title('Comparison',size=15)

        ax2 = fig.add_subplot(212)
        m = Basemap(projection='merc', llcrnrlat=0, urcrnrlat=75., llcrnrlon=-100, urcrnrlon=30, resolution='l')
        m.drawcoastlines()
        m.drawparallels(np.arange(np.min(Lats),np.max(Lats),10.),labels=[True,False,True,False])
        m.drawmeridians(np.arange(np.min(Lons),np.max(Lons),20.),labels=[False,False,False,True])
        m.drawmapboundary(fill_color='lightgray')
        xs,ys = m(lon,lat)
        sc=m.scatter(xs,ys,c=diff)
        fig.colorbar(sc, ax=ax2, format='%.0e', shrink=0.8)
        ax2.set_title('Model - Data',size=15)

        plt.savefig(plotdir + 'ComparisonData/FllTimeAvg' + FlowType,bbox_inches='tight',dpi=100)
        plt.close()        
  
#Compare_to_data()




#def TRAJplot_trajectories():
#    #Function to plot some trajectories as an example for Lagrangian Spaghetti
#    
#    d = 'Traj_SN'
#    f = 'FullRun2_2Traj3daily_kernelSN_tau0.1_gridnum'
#    datadir =TRAJpath+d+'/'
#    
#    k=3
#    r=1110
#    pfile = datadir + f + str(k) + '.nc'
#    data = Dataset(pfile, 'r')
#    lon=data.variables['lon'][100::r,:]
#    lat=data.variables['lat'][100::r,:]
#    
#    plt.figure(figsize=(12,10))
#    m = Basemap(projection='merc', llcrnrlat=0, urcrnrlat=60, llcrnrlon=280, urcrnrlon=360, resolution='l')
#    m.drawcoastlines()
#    m.fillcontinents(color='dimgrey')
#    m.drawmapboundary(fill_color='lightgray')
#    m.drawcoastlines()
#    m.drawparallels(np.arange(0,60,10.),labels=[True,False,True,False])
#    m.drawmeridians(np.arange(280,360,20.),labels=[False,False,False,True])
#    plt.title('Lagrangian Spaghetti',size=15)
#    
#    for i in range(len(lon)):
#        lo=lon[i]
#        la=lat[i]
#        xs,ys = m(lo,la)
#        m.plot(xs,ys)
#    
#    plt.savefig(plotdir + 'Spaghetti.jpg',bbox_inches='tight')
