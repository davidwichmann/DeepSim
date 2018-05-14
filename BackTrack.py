#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Global simultion of different flow models.

@author: wichmann
"""
import numpy as np
from parcels import FieldSet, ParticleSet, ScipyParticle, JITParticle, Field, ErrorCode, AdvectionRK4_3D
from argparse import ArgumentParser
from datetime import timedelta
from netCDF4 import Dataset
from datetime import datetime

datadir = '/Users/wichmann/sources/PhysicalEffects/data/'
outputdir = "./Output/"
griddir = '/Users/wichmann/sources/PhysicalEffects/ParticleGrid/Global10_grid/'
landfile=datadir+'vel-3day/uvel/uvel02432.nc'
dt=-0.1 #Time step in hours

Nyears = 4 #Simulation time

coords = [4.25,79,2500]


def GetOFESLandArray(filename, fieldname):
    """
    Function to return a Field with 1's at land points and 0's at ocean points
    :param f: a field of the .nc file, but not a parcels field object! For OFES, land points are masked. This is used here!
    """    
    pfile = Dataset(filename, 'r')
    Lon = pfile.variables['lon'][:]
    Lat = pfile.variables['lat'][:]
    f = pfile.variables[fieldname][:]
    f = f[0,0,:,:]
    Landmask= np.ma.getmask(f)
    Land=Field('Land',Landmask,transpose=False,lon=Lon,lat=Lat)
    return Land

def rise(particle, fieldset, time, dt):
    d = particle.depth + fieldset.wrise * dt
    if d>0:
        particle.depth = 0
    else:
        particle.depth =d
    
def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()

def periodicBC(particle, fieldset, time, dt):
    if particle.lon < 0.:
        particle.lon += 360.
    elif particle.lon > 360.:
        particle.lon -= 360.

def p_advect(ptype='JITParticle',outname='noname',vrise=0.01):
    """
    Main function for execution
    :outname: name of the output file. Note that all important parameters are also in the file name.
    :kernel: Choices are PN, PB, MN, MB, IN, IB, SN, SB. This is passive, mixing, inertial, stokes, each in normal and brownian mode
    :D: Diffusion constant
    :stokestau: specify tau for inertial particles
    :posidx: Execution is manually parallelized over different initial position grids. These are indexed.
    """
    
    print '-------------------------'
    print 'Start run... Parameters: '
    print '-------------------------'
    print 'ptype: ', ptype
    print 'Rise velocity: ', vrise
    print 'griddir: ', griddir
    print 'landfilename: ', landfile
    print '-------------------------'
    
    [lons,lats,depths] = coords 
    time = [datetime(2005,1,1)-timedelta(days=3)*i for i in range(Nyears*120)]
    lons = [lons]*len(time)
    lats = [lats]*len(time)
    depths = [depths]*len(time)

    outfile = outputdir + outname + "vrise_" + str(vrise) + '_lons0_' + lons + '_lats0_' + lats

    filenames = {'U': datadir + 'U_purely_zonal-ORCA025_grid_U.nc4',
                 'V': datadir + 'V_purely_zonal-ORCA025_grid_V.nc4',
                 'W': datadir + 'V_purely_zonal-ORCA025_grid_V.nc4',
                 'mesh_mask': datadir + 'mesh_mask.nc4'}

    variables = {'U': 'U', 'V': 'V','W': 'W'}

    dimensions = {'U': {'lon': 'nav_lon_u', 'lat': 'nav_lat_u', 'depth': 'depthu'},
                  'V': {'lon': 'nav_lon_v', 'lat': 'nav_lat_v', 'depth': 'depthv'},
                  'W': {'lon': 'nav_lon_v', 'lat': 'nav_lat_v', 'depth': 'depthw'}}

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)
    
    kernels = AdvectionRK4_3D
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)
    
    fieldset.add_field(GetOFESLandArray(landfile, 'uvel'))
    fieldset.U.set_scaling_factor(0.01)
    fieldset.V.set_scaling_factor(0.01)
    
    fieldset.add_periodic_halo(zonal=True)

    #Add constant
    fieldset.vrise=vrise
    
    pset = ParticleSet(fieldset=fieldset, pclass=ptype, lon=lons, lat=lats, depth=depths, time=time)
    kernels+= pset.Kernel(periodicBC)
    kernels+= pset.Kernel(rise)
    
    pset.execute(kernels, runtime=timedelta(days=Nyears*365), dt=timedelta(minutes=60), output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(days=3)),recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    
    
if __name__=="__main__":
    ptype = {'scipy': ScipyParticle, 'jit': JITParticle}
    p = ArgumentParser(description="""Global advection of different particles""")
    p.add_argument('-ptype', '--ptype',choices=('scipy', 'jit'), nargs='?', default='jit',help='Execution mode for performing computation')
    p.add_argument('-name', '--name', default='noname',help='Name of output folder')
    p.add_argument('-stokestau', '--stokestau', type=float,default=1.0,help='Characteristic time')
    p.add_argument('-kernel', '--kernel', default='PN',help='Way of advecting particles')
    p.add_argument('-Diff', '--Diff', type=float,default=10.,help='Diffusion Constant')
    p.add_argument('-posidx', '--posidx', type=int,default=0,help='Label of Lon/Lat initial array')
    p.add_argument('-loc', '--loc', default=None,type=int,help='Initial location. 1...41') 
    args = p.parse_args()
    p_advect(ptype=ptype[args.ptype],outname=args.name,kernel=args.kernel,stokestau=args.stokestau, D=args.Diff, posidx=args.posidx, loc=args.loc)
