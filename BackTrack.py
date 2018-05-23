#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Backwards-simulation of particles ending up at Hausgarten observatory site with nemo 1/12 degree data.
"""

import numpy as np
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, Field, ErrorCode, AdvectionRK4
from argparse import ArgumentParser
from datetime import timedelta
from netCDF4 import Dataset
from datetime import datetime
from glob import glob

datadir = '/data2/imau/oceanparcels/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/'
outputdir = "/scratch/wichm003/"

timestep=-5 #Time step in minutes
Nyears = 4 #Simulation time
Nparticles = 120 #Number of particles (released at different times)

#Coordinates [lon, lat, depth] in [deg, deg, m]
coords = [4.25,79,2500]

def rise(particle, fieldset, time, dt): #Kernel to let the particle rise with rise velocity wrise
    d = particle.depth + fieldset.wrise * dt
    if d<0.6: #This is roughly the smallest depth for nemo data
        particle.depth = 0.6
    else:
        particle.depth =d
    
def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()

def periodicBC(particle, fieldset, time, dt):
    if particle.lon > 180:
        particle.lon -= 360

def p_advect(ptype='JITParticle',outname='noname',wrise=0.001):

    print '-------------------------'
    print 'Start run... Rise velocity: ', wrise
    print '-------------------------'
    
    [lons,lats,depths] = coords 
    outfile = outputdir + outname + "wrise_" + str(wrise) + '_lons0_' + str(lons) + '_lats0_' + str(lats) + '_depth_' + str(depths)

    time = [datetime(2005,12,1)-timedelta(days=5)*i for i in range(120)]
    lons = [lons]*len(time)
    lats = [lats]*len(time)
    depths = [depths]*len(time)

    ufiles = sorted(glob(datadir+'means/ORCA0083-N06_200?????d05U.nc'))
    vfiles = sorted(glob(datadir+'means/ORCA0083-N06_200?????d05V.nc'))

    filenames = {'U': ufiles,
                 'V': vfiles,
                 'mesh_mask': '/data2/imau/oceanparcels/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/domain/coordinates.nc'}

    variables = {'U': 'uo', 'V': 'vo'}

    dimensions = {'U': {'lon': 'nav_lon', 'lat': 'nav_lat', 'time': 'time_centered', 'depth': 'depthu'},
                  'V': {'lon': 'nav_lon', 'lat': 'nav_lat', 'time': 'time_centered', 'depth': 'depthv'}}

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)    
    fieldset.wrise=wrise
    pset = ParticleSet(fieldset=fieldset, pclass=ptype, lon=lons, lat=lats, depth=depths, time=time) 

    kernels = AdvectionRK4
    kernels+= pset.Kernel(rise)    
    kernels+= pset.Kernel(periodicBC)    
    
    pset.execute(kernels, runtime=timedelta(days=Nyears*365), dt=timedelta(minutes=timestep), output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(days=5)),recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    
    
if __name__=="__main__":
    ptype = {'scipy': ScipyParticle, 'jit': JITParticle}
    p = ArgumentParser(description="""Global advection of different particles""")
    p.add_argument('-ptype', '--ptype',choices=('scipy', 'jit'), nargs='?', default='jit',help='Execution mode for performing computation')
    p.add_argument('-name', '--name', default='noname',help='Name of output folder')
    p.add_argument('-wrise', '--wrise', type=float,default=0.01,help='Rise velocity')
    args = p.parse_args()
    p_advect(ptype=ptype[args.ptype],outname=args.name,wrise=args.wrise)
