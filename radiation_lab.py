from numpy import *
from scipy import *
from matplotlib.pyplot import *
import scipy.stats as stat
from math import *

def prelim():
    x=([10,7,11,8,6,6,3,5,8,6,7,6,7,10,9,3,6,13,6,11,10,5,4,9,8,7,9,12,5,4,11,2,6,8,1,2,5,5,6,4,10,8,8,11,8,12,13,7,9,5,5,4,6,9,6,9,9,6,9,11,9,7,11,10,6,9,10,5,7,3,6,13,11,7,10,7,6,6,6,13,4,8,11,10,12,8,7,11,6,5,4,13,5,10,11,12,6,7,9,7])
    x.sort()
    l=mean(x)
    print 'Mean of x is '+str(l)
    hist, edges1=histogram(x,bins=13)
    lh=len(hist)
    stand=std(x)
    poison=[]
    for i in range(lh):
        b=1+i
        v=exp(-l)*l**b/factorial(b)
        poison.append(v)
    #for i in range(lh):
     #   b=hist[i]
      #  v=(1/(sqrt(2*pi)*stand))*exp(-(b-l)**2/(2*stand**2))
       # poison.append(v)
    chi=[]
    base=linspace(1,13,13)
    poison=array(poison)
    npoison=poison*100
    for i in range(lh):
        n=(hist[i]-npoison[i])**2/npoison[i]
        chi.append(n)
    chi=array(chi)
    chisqr=sum(chi)
#chisqr=sum((hist-npoison)**2)/stand
    print 'Chi Squared is '
    print chisqr
    alpha=chisqr/10
    print 'Alpha is '+str(alpha)
    figure(1)
    bar(edges1[:-1], hist, width=edges1[1:]-edges1[:-1])
    plot(base,npoison,'r')
    xlim([1,13])
    ylabel("Counts")
    xlabel('Number Counted')
    title('Histogram of Counts')

    

def deadtime():
    t=300.
    taum=520e-6
    n1=([19138.,19332.,19363.,19239.,19269.])
    n12=([122410.,121496.,122531.,121892.,122075.])
    n2=([109306.,108943.,109655.,108458.,111237.])
    n1a=mean(n1)
    n12a=mean(n12)
    n2a=mean(n2)
    taus=[]
    for s in range(5):
        n1[s]=n1[s]/t
        n2[s]=n2[s]/t
        n12[s]=n12[s]/t
        tau=(1/n12[s])*(1-sqrt(1-((n12[s]*(n1[s]+n2[s]-n12[s]))/(n1[s]*n2[s]))))
        taus.append(tau)
    taus=array(taus)
    standtaus=std(taus)
    na=n1a/t
    nb=n2a/t
    nab=n12a/t
    tau=(1/nab)*(1-sqrt(1-((nab*(na+nb-nab))/(na*nb))))
    print 'Tau is '+str(tau)
    print standtaus
    error_tau=abs(tau-taum)/taum
    print 'the error is tau is'+str(error_tau)


def lead():
    x=([0,1.22e-1,5.17e-1,12.00e-1,15.54e-1,22.26e-1,25.63e-1,30.34e-1],[0,97.21e-3,422.75e-3,1199.59e-7,1429.21e-7,2206.86e-3,2432.71e-3,2418.26e-3],[1010.,1011.,999.,1000.,999.,1000.,1001.,1000.],[123.,115.,185.,349.,463.,814.,1011.,1304.])
    background=43/100
    gamma=0.662e6
    y=len(x[2])
    io=[]
    for i in range(y):
        m=x[2][i]/x[3][i]
        m=m-background
        io.append(m)
    io=array(io)
    z=[]
    for i in range(y):
        lnx=log(io[i])
        z.append(lnx)
    z=array(z)
    figure(2)
    plot(x[0],z)
    xlabel('Thickness (cm)')
    ylabel('ln(I)')
    title('Graph of All Measurements')
    figure(3)
    plot(x[0][1:],z[1:],label='Observed')
    pofit,si=polyfit(x[0],z,1,cov=True)
    pofit=array(pofit)
    print si
    print pofit
    r=pofit[0]
    w=linspace(0,3,8)
    fitline=[]
    for i in w:
        q=pofit[0]*i+pofit[1]
        fitline.append(q)
    fake=[]
    for i in w:
        q=-1.16*i+2.1
        fake.append(q)
    fitline=array(fitline)
    stand=std(z)
    dmax=((z[0]+stand)-(z[-1]-stand))/(x[0][0]-x[0][-1])
    dmin=((z[0]-stand)-(z[-1]+stand))/(x[0][0]-x[0][-1])
    dslope=(dmax-dmin)/2
    print dslope
    print "Standars deviation is "+str(stand)
    print "slope is "+str(pofit)
    plot(w,fitline, label='Theoretical')
    xlabel('Thickness (cm)')
    ylabel('ln(I)')
    title('Thickness > 0')
    legend(loc='upper right')
    #plot(w,fake,'r')
    chisqr=sum((fitline-z)**2)/stand**2
    print "Chi Squared is "+str(chisqr)
    #print z
    print io


def zshield():
    x=(['shiny up with lead','shiny down with lead','shiny up no lead','shiny down no lead',],[1008.,1000.,1004.,1002.,1006.,999.],[169.,131.,151.,128.,171.,172.])
    z=[]
    for i in range(6):
        w=x[1][i]/x[2][i]
        z.append(w)
    print z

def whole():
    print 'PRELIM DATA'
    print ''
    prelim()
    print ''
    print 'DEADTIME DATA'
    print ''
    deadtime()
    print ''
    print 'LEAD DATA'
    print ''
    lead()
