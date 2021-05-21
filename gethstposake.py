import numpy as np
from astropy import coordinates,units
def gethstposake(header,t):

    meananom=header['MEANANOM']+2*np.pi*(header['FDMEANAN']*(t-header['EPCHTIME'])+.5*header['SDMEANAN']*(t-header['EPCHTIME'])**2)
    nu=meananom+np.sin(meananom)*(header['ECCENTX2']+header['ECBDX3']*np.cos(meananom)**2-header['ECBDX4D3']*np.sin(meananom)**2+header['ESQDX5D2']*np.cos(meananom))
    r=header['SEMILREC']/(1+header['ECCENTRY']*np.cos(nu))
    bigom=2*np.pi*(header['RASCASCN']+header['RCASCNRV']*(t-header['EPCHTIME']))
    littleom=2*np.pi*(header['ARGPERIG']+header['RCARGPER']*(t-header['EPCHTIME']))

    x=r*(np.cos(bigom)*np.cos(littleom+nu)-header['COSINCLI']*np.sin(bigom)*np.sin(littleom+nu))
    y=r*(np.sin(bigom)*np.cos(littleom+nu)+header['COSINCLI']*np.cos(bigom)*np.sin(littleom+nu))
    z=r*header['SINEINCL']*np.sin(littleom+nu)
    c=coordinates.EarthLocation(x*units.m,y*units.m,z*units.m)
    return c
