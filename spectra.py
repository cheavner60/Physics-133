from numpy import *
from math import *
from scipy import *
from matplotlib.pyplot import *
import scipy.stats as stat
from astropy.io import ascii

def analyze(file,c=1.4463489,color=False):
    data=ascii.read(file)
    data=array(data)
    deg=data['Theta']
    min=data['min']
    side=data['side']
    order=data['Order']
    lh=len(min)
    rmin=min*0.000290888
    rdeg=deg*(pi/180)
    theta=rmin+rdeg
    theta=abs(theta-c)
    lefto=[]
    righto=[]
    lor=[]
    ror=[]
    for i in range(lh):
        if side[i]==0:
            x=theta[i]
            z=order[i]
            lor.append(z)
            lefto.append(x)
        if side[i]==1:
            x=theta[i]
            z=order[i]
            ror.append(z)
            righto.append(x)
    if file=='helium.txt':
        wl=data['Lambda']
        lwl=[]
        rwl=[]
        for i in range(lh):
            if side[i]==0:
                lwl.append(wl[i])
            if side[i]==1:
                rwl.append(wl[i])
        left=[lefto,lwl,lor]
        right=[righto,rwl,ror]
        left=array(left)
        right=array(right)
        return right, left
    else:
        if color==False:
            color=array([1,2,3,1,2,3,1,2,3])
        else:
            color=array(color)
        left=[lefto,lor,color]
        right=[righto,ror,color]
        left=array(left)
        right=array(right)
        return right, left

def spacing(file):
    right, left =analyze(file)
    lh=len(right[0])
    sig=pi/5400.
    sig=sig/2
 #   print 'Error in theta is '+str(sig)
    theta=0.5*(right[0]+left[0])
    delta=0.5*(right[0]-left[0])
    error=[]
    d=[]
    for i in range(lh):
        dv=(right[2][i]*right[1][i])/sin(theta[i])
        dw=(((right[2][i]*right[1][i])*(cos(theta[i])))/((sin(theta[i]))**2))**2*sig**2
        error.append(dw)
        d.append(dv)
    error=array(error)
    d=array(d)
    d=d/error
    d=sum(d)
    error=1/error
    error=sum(error)
    d=d/error
    error=1/error
    d=d*1e-10
    error=sqrt(error)
    error=error*1e-10
#    print 'Error in diffraction constant is '+str(error)
    return d, error

def spectra(file='hydrogen.txt', file1='helium.txt',c=1.4463489,color=False):
    right, left=analyze(file)
    d, error=spacing(file1)
    sig=pi/5400
    print 'Error in theta is '+str(sig)
    print 'Diffration constant is '+str(d)
    print 'Error in diffraction constant is '+str(error)
    lh=len(right[0])
    theta=0.5*(right[0]+left[0])
    lam1=[]
    lam2=[]
    lam3=[]
    error1=[]
    error2=[]
    error3=[]
    for i in range(lh):
        if right[2][i]==1:
            lam=d*sin(theta[i])/right[1][i]
            dw=(((d*cos(theta[i])/right[1][i])))**2*sig**2+(((sin(theta[i])/right[1][i])))**2*error**2
            error1.append(dw)
            lam1.append(lam)
        if right[2][i]==2:
            lam=d*sin(theta[i])/right[1][i]
            dw=(((d*cos(theta[i])/right[1][i])))**2*sig**2+(((sin(theta[i])/right[1][i])))**2*error**2
            error2.append(dw)
            lam2.append(lam)
        if right[2][i]==3:
            lam=d*sin(theta[i])/right[1][i]
            dw=(((d*cos(theta[i])/right[1][i])))**2*sig**2+(((sin(theta[i])/right[1][i])))**2*error**2
            error3.append(dw)
            lam3.append(lam)
    error1=array(error1)
    error1=mean(error1)
    error1=sqrt(error1)
    lam1=array(lam1)
    lam2=array(lam2)
    lam3=array(lam3)
    lam1=mean(lam1)
    lam2=mean(lam2)
    lam3=mean(lam3)
    lam=array([lam1,lam2,lam3])
    print lam
    n2=2.0
    n1=array([5.0,4.0,3.0])
    hl=len(n1)
    rh=[]
    dx=[]
    for i in range(hl):
        v=((1/n2)**2-(1/n1[i])**2)
        #rhv=1/(v*lam[i])
        rhv=(right[1][i]/(v*d*sin(theta[i])))
        mom=(right[1][i]/(v*d**2*sin(theta[i])))**2*error**2+(((right[1][i]*cos(theta[i]))/(v*d*(sin(theta[i]))**2)))**2*sig**2
        dx.append(mom)
        rh.append(rhv)
    dx=array(dx)
    rh=array(rh)
    rh=rh/dx
    rh=sum(rh)
    dx=1/dx
    dx=sum(dx)
    rh=rh/dx
    dx=1/dx
    dx=sqrt(dx)
    print 'The Rydberg constant is '+str(rh)
    print 'Error in rydberg is '+str(dx)

def laser(ltheta=47.0*pi/108., rtheta=6001.0*pi/10800.0,c=1.5563489):
    d, error=spacing()
    sig=pi/5400.
    ltheta=abs(ltheta-c)
    rtheta=abs(rtheta-c)
    theta=0.5*(rtheta+ltheta)
    lam=d*sin(theta)
    print 'Make sure valuse are floats'
    print 'Wavelength of laser beam is '+str(lam)
