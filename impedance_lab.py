from numpy import *
from math import *
from scipy import *
from matplotlib.pyplot import *
import scipy.stats as stat
from astropy.io import ascii

def dc(file='dc_data.txt'):
    data=ascii.read(file)
    data=array(data)
    r=data['R']
    v=data['V']
    i=v/r
    figure()
    plot(i[6:],v[6:])
    pofit,si=polyfit(i[6:],v[6:],1,cov=True)
    xinter=-pofit[1]/pofit[0]
    xint=pofit
    ylabel('Voltage')
    xlabel('Current')
    title('IV for Battery')
    print pofit
    print si
    print 'the x intercepts is '+str(xint)

def ac(file,type,xlimit=False,ylimit=False):
    r=1000.
    data=ascii.read(file)
    data=array(data)
    freq=data['Freq']
    vb=data['VB']
    vr=data['VR']
    p=data['Phase']
    datasize=len(freq)
    I=[]
    for i in range(datasize):
        n=vr[i]/r
        m=vb[i]/n
        I.append(m)
    z=array(I)
    figure()
    if type is 'c':
# Inverts Impedance
        plot(freq,1/z,label='Observed')
        if xlim!=False:
            xlim(0,xlimit)
        if ylim!=False:
            ylim(0,ylimit)
        pofit,si=polyfit(freq,1/z,1,cov=True)
        print 'the capacitor value is '+str(pofit[0])
        print max(1/z)
        w=[]
        for m in range(datasize):
            x=freq[m]
            rand=pofit[0]*x + pofit[1]
            w.append(rand)
        w=array(w)
        plot(freq,w,label='Theoretical')
        title('$\omega^2$ vs 1/Z')
        xlabel('$\omega^2$')
        ylabel('1/Z')
        legend(loc='upper right')
        print pofit
        print si
    if type is 'inv':
# Inverts Frequency
        cs=(1/2*freq)
        cs=cs**2
        plot(cs[1:-1],z[1:-1], label='Observed')
        if xlim!=False:
            xlim(0,xlimit)
        if ylim!=False:
            ylim(0,ylimit)
        xlabel('$\omega^2$')
        ylabel('Z')
        title('Capacitence In Series With A Resistor')
        pofit,si=polyfit(cs[1:-1],z[1:-1],1,cov=True)
        cap=1/pofit[0]
        lw=cs[1:-1]
        lz=z[1:-1]
        pl=sum(lw*lz)
        lx2=sum(lw*lw)
        lx1=sum(lw)
        ly1=sum(lz)
        la=(pl-lx1*ly1)/(lx2-lx1**2)
        print 'error in Capacitence'
        print 1/(lx2-lx1**2)
        print lx2/(lx2-lx1**2)
        print 'the capacitor value is '+str(cap)
        y=[]
        for m in range(datasize):
            x=1/(2*freq[m])
            rand=2.608e6*x + 1.146e2
            y.append(rand)
        y=array(y)
        plot(cs[1:-1],y[1:-1],label='Theoretical')
        xlabel('1/$\omega^2$')
        ylabel('Z')
        legend(loc='upper right')
        b=y[1:-1]
        d=z[1:-1]
        lh=len(d)
        stand=std(d)
        chi=[]
        for i in range(lh):
            n=(d[i]-b[i])**2/stand
            chi.append(n)
        chi=array(chi)
        chisqr=sum(chi)
        print 'Chi Squared is '
        print chisqr
        print pofit
        print si
    if type is 'norm':
        plot(freq,z)
        if xlim!=False:
            xlim(0,xlimit)
        if ylim!=False:
            ylim(0,ylimit)
        pofit,si=polyfit(freq,z,1,cov=True)
        print pofit
        print si
    if type is 'e':
# Specifically for e
        plot(freq,1/z)
        xlabel('$\omega^2$')
        ylabel('1/Z')
        title('Inductor and Capacitor in Series')
        invz=1/z
        print max(1/z)
        xlim(0,100000)
        ylim(0,0.0075)
        lh=len(z)
        cz=[]
        cw=[]
        lz=[]
        lw=[]
        for i in range(lh):
            if freq[i]<15000 and freq[i]>1500:
                cz.append(z[i])
                cw.append(freq[i])
            if freq[i]>70000:
                lz.append(z[i])
                lw.append(freq[i])
        cz=array(cz)
        lz=array(lz)
        cw=array(cw)
        lw=array(lw)
        cpofit,csi=polyfit(cw,1/cz,1,cov=True)
        lpofit,lsi=polyfit(lw,lz,1,cov=True)
        llen=len(lw)
        clen=len(cw)
        pl=sum(lw*lz)
        lx2=sum(lw*lw)
        lx1=sum(lw)
        ly1=sum(lz)
        la=(pl-lx1*ly1)/(lx2-lx1**2)
        print 'error in inductance'
        print 1/(lx2-lx1**2)
        cx=[]
        lx=[]
        for m in range(clen):
            crand=cpofit[0]*cw[m] + cpofit[1]
            cx.append(crand)
        for n in range(llen):
            lrand=lpofit[0]*lw[n]+lpofit[1]
            lx.append(lrand)
        lx=array(lx)
        cx=array(cx)
        figure()
        plot(cw,1/cz,label='Observed')
        plot(cw,cx,label='Theoretical')
        xlabel('$\omega^2$')
        ylabel('1/Z')
        ylim(0,0.00016)
        title('Capacitence $\omega^2$ vs 1/Z')
        legend()
        figure()
        plot(lw,lz,label='Observed')
        plot(lw,lx,label='Theoretical')
        ylim(1500,4300)
        xlabel('$\omega^2$')
        ylabel('Z')
        title('Inductance $\omega^2$ vs Z')
        legend()
        print 'CAPACITENCE'
        print cpofit
        print csi
        print 'INDUCTANCE'
        print lpofit
        print lsi

def diode(file='diode.txt'):
    data=ascii.read(file)
    data=array(data)
    r=data['R']
    v=data['V']
    lh=len(r)
    i0=[]
    for n in range(lh):
        i=19.2/r[n]
        i0.append(i)
    figure()
    i0=array(i0)
    plot(v[3:-3],i0[3:-3])
    xlabel('V')
    ylabel('ln(I)')
    title('Diode Voltage vs ln(I)')
