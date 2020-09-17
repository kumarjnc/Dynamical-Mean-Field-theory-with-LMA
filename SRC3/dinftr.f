	Subroutine  Dinftransp(flag)

c	flag=0  => only DC 

	include 'Glob_decl'
	integer locate,lj,ulim,N1,jll,jul,iter,flag
	character*72 linemc
	real*8 cond1(2*Nm),cond2(2*Nm),wl(2*Nm),fermic,rho,wij
	real*8 conv,a2,a4,zx,zy
	complex*16 Gc_hc(-Nm:Nm),Gc_bl(-Nm:Nm),hybr,phyb,nhyb
	complex*16 gamma(-Nm:Nm),gammal,Gcl_hc,Gcl_bl,gup,gdn,
     .		  sem,gauss,g0
	complex*16 dgammaj,wz
	real*8 sgn,m2_h,m4_h

	external sem,gauss,fermic,sgn,locate


	include'Glob_cons'	
	dfac=0.5d0


c	HCL	
	a2=t**2/2.d0
	a4=t**4/2.d0
	m2_h=a2
	m4_h=3.d0*t**2/4.d0

	x0=0.5d0*U*mu0
	ep_f=-0.5d0*U*(1.d0-asym_p)
	ei0=ep_f+U*n0/2.d0

        do i=-N,N
          iter=1
          conv=1.d0
          z1=w(i)+ii*1.D-10
	  hybr=hyb(i)
          phyb=hybr
          z3=sigma(1,i)
          gdn=z1-ec-V**2/(z1-ei0-x0-z3)
          z3=sigma(2,i)
          gup=z1-ec-V**2/(z1-ei0+x0-z3)
          g0=z1-ec
          do while(conv.gt.1.D-10.and.iter.le.100000)
!            gammal=(2.d0*gup*gdn-hybr*(gup+gdn))/
!     .            (gup+gdn-2.d0*hybr)
 
           gammal=hybr+(2.d0*(g0-hybr)*(gup*gdn 
     .              -hybr*(gup+gdn)+hybr**2))/  
     .            ((1.d0-dop)*(g0-hybr)* 
     .             (gup+gdn-2.d0*hybr)+2.d0*dop* 
     .           (gup*gdn-hybr*(gup+gdn)+hybr**2))
            if(zabs(gammal).le.1.D3) then
              nhyb=(gammal-1.d0/gauss(gammal))
            else
               nhyb=a2/gammal+a4/gammal**3   ! From moment expansion
            end if
            hybr=dfac*phyb+(1.d0-dfac)*nhyb
            conv=zabs(hybr-phyb)
            phyb=hybr
            if(iter.eq.100000) write(6,*) 'MAX LIM',i,conv
            iter=iter+1
          end do
          Gc_hc(i)=1.d0/(gammal-hybr)
	  gamma(i)=gammal
	end do
        open(unit=27,file='sigcpa_test.dat',status='unknown')
        do i=-N,N
        if(w(i).gt.(0.001d0*V**2*Zfac).and.w(i).le.0.1*V**2*Zfac)then
        write(27,*)w(i)/(V**2*Zfac),dimag(gamma(i)),dreal(gamma(i))
        endif
        enddo
        close(27)

	write(6,*) Gc_hc(0),gamma(0)

	  write(6,*) 'HC done'
c	BL	
	a2=t**2/4.d0
	a4=t**4/16.d0

        do i=-N,N
          iter=1
          conv=1.d0
          z1=w(i)+ii*1.D-10
	  hybr=hyb(i)
          phyb=hybr
          z3=sigma(1,i)
          gdn=z1-ec-V**2/(z1-ei0-x0-z3)
          z3=sigma(2,i)
          gup=z1-ec-V**2/(z1-ei0+x0-z3)
          g0=z1-ec
          do while(conv.gt.1.D-10.and.iter.le.100000)
!            gammal=(2.d0*gup*gdn-hybr*(gup+gdn))/
!     .            (gup+gdn-2.d0*hybr)

           gammal=hybr+(2.d0*(g0-hybr)*(gup*gdn 
     .              -hybr*(gup+gdn)+hybr**2))/  
     .            ((1.d0-dop)*(g0-hybr)* 
     .             (gup+gdn-2.d0*hybr)+2.d0*dop* 
     .           (gup*gdn-hybr*(gup+gdn)+hybr**2))
            if(zabs(gammal).le.1.D3) then
              nhyb=(gammal-1.d0/sem(gammal))
            else
               nhyb=a2/gammal+a4/gammal**3   ! From moment expansion
            end if
            hybr=dfac*phyb+(1.d0-dfac)*nhyb
            conv=zabs(hybr-phyb)
            phyb=hybr
            if(iter.eq.100000) write(6,*) 'MAX LIM',i,conv
            iter=iter+1
          end do
          Gc_bl(i)=1.d0/(gammal-hybr)
!	  gamma(i)=gammal
	end do

	if(temp.ge.0.d0) then	!  DC conductivity

        if(temp.eq.0.d0) then
           r4=(-dimag(Gc_hc(0)))/dimag(gamma(0))
           r4=r4/(2.d0*pi**2)
           r6=0.d0
	   !write(6,*) 'Loop execution'
	else
 	r1=0.d0
	r4=0.d0
        r6=0.d0
        r10=0.0d0
	do i=-N,N
          z1=w(i)+ii*1.D-10 
	  rho=-dimag(gauss(z1))/pi
	  r2=fermic(w(i),beta)
	  r3=beta*r2*(1.d0-r2)
	  r1=r1+rho**2*r3*dw(i)*t**2
	  r5=zabs(gamma(i))
          dgammaj=dconjg(gamma(i))
	  gammal=gamma(i)
	  rho=-dimag(Gc_hc(i))/pi
!	  if(r5.le.1.D3) then
!	  r2=dreal(pi*rho/dimag(gamma(i)) + 2.d0*(1.d0-gamma(i)*
!     .		  Gc_hc(i)))
!          else
!	  z2=(1.d0+m2_h*(1.d0/gammal**2+1.d0/(dgammaj*gammal)+
!     .		1.d0/dgammaj**2))/(gammal*dgammaj)-
!     .	   	2.d0*(m2_h+m4_h/gammal**2)/gammal**2
!          r2=dreal(z2)
!          
!     	  end if
          
	  r2=(pi*rho/dimag(gamma(i)))
          
          zx=dreal(gamma(i))
          zy=dimag(gamma(i))
          wz=ii*Gc_hc(i)/dsqrt(pi)
          r8=3.D0*pi/(4.D0*zy)*(r2*0.5D0*pi)+zy/2.D0+
     .           dsqrt(pi)/2.D0*dreal(((gamma(i))**2-0.5D0)*wz)
          r8=3.0D0*zx*r8-dsqrt(pi)*zy**2*
     .            dreal(gamma(i)*wz)-
     .           dsqrt(pi)*zy**2*dreal(((gamma(i))**2-0.5D0)*
     .           (ii/dsqrt(pi)-gamma(i)*wz))
           r8=r8/(2.D0*pi)/(1.5D0+zy**2)
          
          r4=r4+r2*r3*dw(i)*t**2/(2.d0*pi**2)
         
          r6=r6+w(i)*r2*r3*dw(i)*t**2/(2.d0*pi**2)
          r10=r10+r3*r8*0.5D0*dw(i)
          r11=r11+(w(i))**2*r3*r2*dw(i)
	end do
	end if
	open(unit=11,file='dccond.dat',status='unknown')
	write(11,*) temp,r4,r6/(temp*r4)
	close(11)

	end if


	if(flag.eq.0) return

	write(6,*) 'Now optical conductivity'
	do i=-N,N
	  j=i+N+1
	  wl(j)=w(i)
	end do
	ulim=locate(2*N+1,wl,1,2*N+1,4.0d0)-N-1
	write(6,*) 'ulim=',ulim


	open(unit=10,file='optcond.dat',status='unknown')
	do 200 i=1,ulim
	 z1=zero
	 z3=zero
	 if(temp.eq.0.d0) then
	   jll=-i
	   jul=0
	 else
	   jll=-N
	   jul=locate(2*N+1,wl,1,2*N+1,w(N)-w(i))-N-1
	 end if
	 llim=-N
	 do 201 j=jll,jul
	  wij=w(i)+w(j)
	  do while(wij.gt.w(llim+1))
	  llim=llim+1
	  end do
c	  lj=locate(2*N+1,wl,1,2*N+1,wij)-N-1
	  lj=llim
	  Gcl_hc=Gc_hc(lj)+(Gc_hc(lj+1)-Gc_hc(lj))*((wij-w(lj))/
     .		   (w(lj+1)-w(lj)))
	  Gcl_bl=Gc_bl(lj)+(Gc_bl(lj+1)-Gc_bl(lj))*((wij-w(lj))/
     .		   (w(lj+1)-w(lj)))
	  gammal=gamma(lj)+(gamma(lj+1)-gamma(lj))*   !HC
     .		   ((wij-w(lj))/(w(lj+1)-w(lj)))		! HC
          dgammaj=dconjg(gamma(j))
	  r1=zabs(gammal)
	  r2=zabs(gamma(j))
	  r3=zabs(gammal-gamma(j))
	!  if(r1.le.1.D3.or.r2.le.1.D3) then
	  z2=(dconjg(Gc_hc(j))-Gcl_hc)/(gammal-dgammaj)	! HC
     .		 - (Gc_hc(j)-Gcl_hc)/(gammal-gamma(j))	! HC
        !  else
	 ! z2=(1.d0+m2_h*(1.d0/gammal**2+1.d0/(dgammaj*gammal)+
!     .		1.d0/dgammaj**2))/(gammal*dgammaj)-
!     .	     (1.d0+m2_h*(1.d0/gammal**2+1.d0/(gamma(j)*gammal)+
!     .		1.d0/gamma(j)**2))/(gammal*gamma(j))
!     	  end if
	  z4=(-dimag(Gc_bl(j)))*(-dimag(Gcl_bl))/pi**2	! Bethe
	  r1=fermic(w(j),beta)-fermic(wij,beta)
	  z1=z1+z2*dw(j)*r1
	  z3=z3+z4*dw(j)*r1
 201	 continue
	 cond1(i)=dreal(z1)/(2.d0*pi*pi*w(i))		! HC
	 cond2(i)=dreal(z3)/(w(i))		! Bethe
	 if(w(i).ge.1.D-10) then
	   write(10,*) w(i),cond1(i),cond2(i)
	 end if
 200	continue
	close(10)


	return
	end

c	************************************************************
