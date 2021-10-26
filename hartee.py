#!/usr/bin/python3
import numpy as np
from numpy import tile
from scipy import linalg as LA
import math,cmath
pi = math.pi
cos = math.cos
sin = math.sin
exp=np.exp   
sqrt = math.sqrt
log=np.log

def hartee(ef, ec, ufc, v, nc, nf, uff, m, Nm, flag=True):
    from specfn import specfn
    from commonfun import commonfn
    t=1.0
    a1=t**2/4.0
    a2=t**4/8.0
    alpha=a1
    abeta=a2-a1**2
    dw=0.0001
    mconv=1.0
    miter=1
    pmag=m
    pnf=nf
    pnc=nc
    tol=1**(-5)
    tol2=1**(-3)
    #--------------------------------------
    dcup=specfn(2*Nm).spec.real
    dcdn=specfn(2*Nm).spec.real
    dfup=specfn(2*Nm).spec.real
    dfdw=specfn(2*Nm).spec.real
    dc=specfn(2*Nm).spec.real
    df=specfn(2*Nm).spec.real
    #---------------------------------------
    gcup=specfn(2*Nm).spec
    gcdn=specfn(2*Nm).spec
    gfup=specfn(2*Nm).spec
    gfdw=specfn(2*Nm).spec
    gc=specfn(2*Nm).spec
    gf=specfn(2*Nm).spec
    
    #---------------------------------------
    hyb=specfn(2*Nm).spec
    sbar=specfn(2*Nm).spec
    
    while ((mconv > tol) and (miter < 50)):
        r=0.0
        r1=0.0
        r2=0.0
        r3=0.0
        r4=0.0
        r5=0.0
        for i in range(0, 2*Nm):
            omega=(i-Nm)*dw
            p=commonfn(omega,tol)
            z=p.mycomplex()
            zz=p.sem()
            hyb[i]=((t**2)/4.0)*zz
            gup=z-(ec+ufc*nf)-v**2/(z-(ef+uff*nf/2.0+ufc*nc-uff*m/2.0))
            gdn=z-(ec+ufc*nf)-v**2/(z-(ef+uff*nf/2.0+ufc*nc+uff*m/2.0))
            g0=z-(ec+ufc*nf)
            conv=1.0
            iter=1
            while ((conv > tol) and (iter < 50)):
                phyb=hyb[i]
                gamma=hyb[i]+(2.0*(g0-hyb[i])*(gup*gdn-hyb[i]*(gup+gdn)+hyb[i]**2))/((1.0-dop)*(g0-hyb[i])(gup+gdn-2.0*hyb[i])+2.0*dop*(gup*gdn-hyb[i]*(gup+gdn)+hyb[i]**2))
                if(abs(gamma) < tol2 ):
                    p=commanfn(gamma,0.0)
                    zz2=p.sem()
                    nhyp=gamma-1.0/zz2
                else:
                    nhyb=alpha/gamma-abeta/gamma**3
                hyb[i]=0.5*(phyb+nhyb)
                conv=abs(hyb[i]-phyb)
                iter=iter+1
            sbar[i]=gamma
            gc[i]=1.0/(sbar[i]-hyb[i])
            gfup(i)=1.0/(z-ef-uff*nf/2.0-ufc*nc+uff*m/2.0-v**2/(z-(ec+ufc*nf)-hyb(i)))
            gfdn(i)=1.0/(z-ef-uff*nf/2.0-ufc*nc-uff*m/2.0-v**2/(z-(ec+ufc*nf)-hyb(i)))
            gf[i]=0.5*(gfup[i]+gfdn[i])
            dc[i]=-gc[i].imag/pi
            dfup[i]=-gfup[i].imag/pi
            dfdw[i]=-gfdw[i].imag/pi
            df[i]=-gf[i].imag/pi
            p=commanfn(-gamma,0.0)
            r=r+2.0*dc(i)*dw*p.theta()
            r1=r1+dfup[i]*dw*theta(-w)
            r2=r2+dfdn[i]*dw*theta(-w)
            r3=r3+dc[i]*dw
            r4=r4+df[i]*dw
            
        f=r
        f1=r1
        f2=r2
        nc=f
        nf=f1+f2
        m=f1-f2
        dctot=r3
        dftot=r4
        mconv=c2
        miter=miter+1
        pmag=m
        pnf=nf
        pnc=nc
        
    if(flag==True)
        print('mconv=',mconv)
        p
       if(flag.eq.1) write(6,*)'spectral weight'
       if(flag.eq.1)write(6,*)'nc=',nc,'nf=',nf
       if(flag.eq.1) write(6,*)dctot,dftot
      
                    
            