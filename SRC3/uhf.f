c	THE UHF-2 PART of the ONE-LOOP LEVEL
c	**************************************************************

	real*8 function solveq(flag,flag1)

	include 'Glob_decl'
	real*8 xa,mu,nc0,fv,metins,wi,nup,ndn,
     .		galpha,falpha,ep_f
	integer iter,sp,spj,flag,flag1
	complex*16 gfhf(-Nm:Nm),gfsp,pgc
	character ch

	external metins

	include 'Glob_cons'
	
c	If flag=0 symmetric case
c	else	  asymmetric case
c	If flag1=0 short stdout
c	else	   full output

c	THE UHF C-GREEN'S FUNCTION
	call uhfgc

c	NOW THE F-GREEN'S FUNCTION
	
	r1=0.d0
	r2=0.d0
	r3=0.d0
	do i=-N,N
	   z1=w(i)+ii*1.D-10
	   gfhf(i)=zero
	   gfsp=zero
	   do sp=-1,1,2
	     spj=aint(sp/2.d0+1.5d0)
	     z2=1.0d0/(z1-ei+dfloat(sp)*x
     .	     -V**2/(z1-ec-hyb(i)))
	     Gfscript(spj,i)=z2
	     gfhf(i)=gfhf(i)+0.5d0*z2
	     gfsp=gfsp+sp*z2
	   end do
	   r1=r1+dw(i)*ferm(i)*(-dimag(gfsp)/pi)
	   r2=r2+2.d0*dw(i)*(-dimag(gfhf(i))/pi)*ferm(i)
	   r3=r3+2.d0*dw(i)*(-dimag(gc0(i))/pi)*ferm(i)
	end do
	mu=r1
	mu0=r1
	n0=r2
	nc0=r3

	
	r1=0.d0
	r2=0.d0
	do i=-N,0
	  r1=r1+(-dimag(Gfscript(1,i))/pi)*dw(i)*ferm(i)
	  r2=r2+(-dimag(Gfscript(2,i))/pi)*dw(i)*ferm(i)
	end do
	ndn=r1
	nup=r2

	galpha=nup-ndn
	falpha=nup+ndn
	U0=2.d0*x/galpha
	ch='I'
	wi=0.d0
	fv=metins(0.d0)
	if(fv.ge.0.d0) ch='M'

	ep_f=ei-U0*n0/2.d0
c	r1=1.d0+2.d0*ep_f/U0   ! (r1=asym_p at UHF level)

	write(6,*) 'ei=',ei,' ep_f=',ep_f,' eta=',r1
	

	if(flag1.eq.1) then
	open(unit=11,file='Rho_uhf.dat',status='unknown')
	write(6,*) 'UHF f,c-spec. fn in fort.11'
	do i=-N,N
	   write(11,*) w(i),-dimag(gfhf(i))/pi,-dimag(gc0(i))/pi
	end do
	close(11)


	write(6,400) x,ec,ei,ch
	write(6,500) U0,mu,nc0,n0
	
	open(unit=20,file='Uhf_par.dat',status='unknown')
	write(20,400) x,ec,ei,ch
	write(20,500) U0,mu,nc0,n0
	close(20)

	else

	write(6,500) U0,mu,nc0,n0


	end if


 400	format('x=',f10.6,'  ec=',f10.6,'  ei=',f10.6,'  ch=',a3)
 500	format('U0=',f10.6,'  mu=',f10.6,'  nc0=',f10.6,'  n0=',f10.6)
 600	format(f7.3,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,a3)
 601	format(f10.6,f10.6,f10.6)


	if(flag.eq.0) then
	    solveq=galpha	! Symmetric case
	else
	    solveq=ei*galpha-x*(falpha-1.d0)  !Asymmetric case
        end if


	return
	end

c	************************************************************

	subroutine uhfgc

	include 'Glob_decl'
	real*8 mu
	integer iter,sp
	real*8 wi,pfv,fv
	real*8 conv,a1,a2,alpha,abeta
	complex*16 gcinit(-Nm:Nm),pgc,sbar(-Nm:Nm),gup,gdn,
     .	gamma,phyb,hilbtr,gauss,nhyb,g0

	external hilbtr,gauss

	include 'Glob_cons'

	dfac=0.5d0

	do i=-N,N
	  z1=w(i)+ii*1.d-10
          gcinit(i)=gauss(z1)
          hyb(i)=t**2*gcinit(i)/4.d0
          hyb_G(i)=t**2*gcinit(i)/4.d0
        end do
        open(unit=53,file='hyn_G.dat',status='unknown')
        do i=-N,N
          write(53,*)w(i),-dimag(hyb_G(i))
        enddo
        close(53)
         

c       Calculate second and fourth moments of the Gaussian DOS

        if(idos.eq.1) then
          a1=t**2/2.d0
          a2=t**4*3.d0/4.d0
        else if(idos.eq.2) then
          a1=t**2/4.d0
          a2=t**4/8.d0
        end if

        alpha=a1
        abeta=a2-a1**2



c	UHF
	
c	first find the UHF Gc
	open(unit=10,file='gc0.dat',status='unknown')
	do i=-N,N
	   iter=1
	   conv=1.d0
	   z1=w(i)+ii*1.D-10
	   phyb=hyb(i)
           gup=z1-ec-V**2/(z1-ei+x)
           gdn=z1-ec-V**2/(z1-ei-x)
!           g0=z1-ec
	   do while(conv.ge.1.D-10.and.iter.le.1000000)
            gamma=(2.d0*gup*gdn-hyb(i)*(gup+gdn))/
     .            (gup+gdn-2.d0*hyb(i))
            
!           gamma=hyb(i)+(2.d0*(g0-hyb(i))*(gup*gdn 
!     .              -hyb(i)*(gup+gdn)+hyb(i)**2))/  
!     .            ((1.d0-dop)*(g0-hyb(i))* 
!     .             (gup+gdn-2.d0*hyb(i))+2.d0*dop* 
!     .           (gup*gdn-hyb(i)*(gup+gdn)+hyb(i)**2))
          if(zabs(gamma).le.1.D3) then
               nhyb=(gamma-1.d0/hilbtr(gamma,idos))
            else
               nhyb=alpha/gamma+abeta/gamma**3  ! From moment expansion
            end if
            hyb(i)=dfac*phyb+(1.d0-dfac)*nhyb
            conv=zabs(hyb(i)-phyb)
            phyb=hyb(i)
c           write(6,*) i,iter,conv
            iter=iter+1
          end do
          sbar(i)=gamma
            gc0(i)=1.d0/(sbar(i)-hyb(i))
!          gc0(i)=0.5d0*(1.d0/(gup-hyb(i))+1.d0/(gdn-hyb(i)))
          write(10,*) w(i),-dimag(gc0(i))/pi
        end do
	close(10)



 	return
	end
	

c	************************************************************
	real*8 function bisection(lx,rx,lfv,rfv,fun)
	real*8 lx,rx,lfv,rfv,mx,mfv,fun,h
	integer iter
	include 'Glob_decl'

	external fun
	include 'Glob_cons'
	  
	iter=1
	i=1
	do while(dabs(lfv).ge.1.D-12.and.dabs(rfv).ge.1.D-12
     .		.and.i.le.100000)
	  if(iter.gt.2) iter=1
          mx=0.5d0*(lx+rx)
	  if(iter.eq.2) mx=lx-lfv*(rx-lx)/(rfv-lfv)
	  mfv=fun(mx)
	  if((mfv*lfv).le.0.d0) then
	     rx=mx
	     rfv=mfv
	  else
	     lx=mx
	     lfv=mfv
	  end if
	  iter=iter+1
	  i=i+1
	  if(i.eq.100000) then
		write(6,*) 'bisection max reached'
	  end if
 	end do

        bisection=lx
        if(dabs(rfv).lt.dabs(lfv)) bisection=rx

        return
        end


c	************************************************************

	real*8 function metins(wi)
	include 'Glob_decl'
	real*8 wi,a00,a1,a2,q,r,fv

	r1=wi-ec-V**2/(wi-ei+x)
	r2=wi-ec-V**2/(wi-ei-x)

	a2=-4.d0*(r1+r2)/t**2
	a1=16.d0*(r1*r2+t**2/4.d0)/t**4
	a00=2.d0*a2/t**2
	q=a1/3.d0-a2**2/9.d0
	r=(a1*a2-3.d0*a00)/6.d0-a2**3/27.d0
	
	if(wi.ne.0.d0) then
	metins=3.d0*dlog(dabs(q))-2.d0*dlog(dabs(r))
	else
	metins=q**3+r**2
	end if


	return
	end
	
c	************************************************************

