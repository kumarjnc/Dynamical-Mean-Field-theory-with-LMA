	Program ASPAM_UHFPi0

c	ASYMMETRIC case
c	USES UHF-2 ALGORITHM in uhf.f
c	Auxiliary routines are in funct.f
c	FIXED ei,ec,x, ep_f=-U/2.0(1.d0-asym_p), UHF Pi0
c	T >= 0

	include 'Glob_decl'
	integer iter,sp,unit1,nloop,nloopmax,init,imax,flag
	integer sflag
	real*8 ep_f,wi,conv,pw_m,lval,findgcf
	real*8 solveq,tolf
	complex*16  pGf(-Nm:Nm)
	character*72 linemc

	external solveq,findgcf

	include'Glob_cons'	
!        dop=0.0d0
	t=1.0d0
	!V=t*dsqrt(0.6d0)
	dfac=0.5d0
        
c	READ FILES


        open(unit=12,file='par.dat',status='unknown')
        read(12,*) ep0,ep2,ep3,ep4,ep5,ep_wm
        read(12,*) b1,b2,db2,b3_n,db3_n,b3_p,db3_p
	read(12,*) b4,db4,wm_l,wm_r,asym_p
	read(12,*) x,dop,ec,ei,init,idos,flag
	read(12,*) temp,frac
	read(12,*) dU,Ufact
        close(12)
        write(6,*) 'Enter temp,frac'
	read(5,*) temp,frac
         V=t*dsqrt(0.6d0+dop*0.0d0)
        ! if(dop.le.0.5d0)then
         !ec=0.5d0
         !else
!         ec=0.3d0+dop*0.0d0
         !endif
!	read(5,*) ei
        !init=2
        !flag=1
	sflag=1
	if((ec.eq.0.d0.and.ei.eq.0.d0).and.asym_p.eq.0.d0) sflag=0

	if(temp.eq.0.d0) then
	  beta=0.d0
	  tolf=1.D-7
	else
	  beta=1.d0/temp
	  tolf=1.D-4
	end if

        call makegrid(1)

	unit1=60
	pw_m=0.d0
	w_m=0.d0
	conv=1.d0	
	nloop=1
	nloopmax=10
	Npi0=N
	do i=-N,N  
	   Gc(i)=zero
	   Gf(i)=zero
	   ssig(i)=zero
	   w_pi(i)=w(i)
	   Pi0(i)=zero
	   Pi0_G(i)=zero
	   do sp=1,2
	     Gfscript(sp,i)=zero
	     sigma(sp,i)=zero
	   end do
	end do



c	LOOPING STARTS HERE  OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
	do while(conv.ge.tolf.and.nloop.le.nloopmax)
c           do while(nloop.le.nloopmax)
	if(nloop.eq.1.and.init.eq.0) then

	  call initialise(flag)

	else if(nloop.eq.1.and.(init.eq.1.or.init.eq.2)) then
	
	  call readinput(flag)

	  if((flag.eq.0.and.temp.eq.0.d0).or.
     .	  (flag.eq.1.and.U.gt.1.d0/dreal(pi0(0)))) then
	     U=1.d0/dreal(pi0(0))-dU
	  end if

	  call grid_tpi
	  write(6,*) 'w_m=',w_m
	  pw_m=w_m

	  ep_f=-U*(1.d0-asym_p)/2.d0
	  lval=findgcf(ep_f)
	  !call Grid_Sig
	  !call Grid_Gf(0)
	  !call interpol(1)
	  lval=findgcf(ep_f)


	  if(init.eq.2) then
	     call dinftransp(0)
	     stop
	  end if

	  call findgfscript(conv)
	  call writeGfscript

	  call takeoutpole(1)
	  call writeGfsclp

	else

 	  call takeoutpole(1)
	  call writeGfsclp

	end if

 	write(6,*) 'FINDING SELF ENERGY'

	call selfenergy(nloop,init,conv)

	do i=-N,N
	  pGf(i)=Gf(i)
	end do

	if(temp.eq.0.d0) then
 	  write(6,*) 'FINDING EP_F'
 	  write(6,*) '				'
 	  write(6,*) '				'

	  ep_f=-U*(1.d0-asym_p)/2.d0
          if(sflag.ne.0) call findepf(ep_f)
	end if

	ep_f=-U*(1.d0-asym_p)/2.d0
	lval=findgcf(ep_f)	! Gives Gf,Gc
	if(nloop.eq.1.and.init.eq.0.and.sflag.ne.0) then
	  !call Grid_Sig
	  !call Grid_Gf(0)
	  !call interpol(1)
	  ep_f=-U*(1.d0-asym_p)/2.d0
	  lval=findgcf(ep_f)
	end if

	if(temp.eq.0.d0.and.sflag.ne.0) then
 	  write(6,*) 'FINDING EP_F AGAIN'
 	  write(6,*) '				'
 	  write(6,*) '				'
	  ep_f=-U*(1.d0-asym_p)/2.d0
          call findepf(ep_f)
	end if


	ep_f=-U*(1.d0-asym_p)/2.d0
	lval=findgcf(ep_f)	! Gives Gf,Gc


	write(6,*) 'FINDING GFSCRIPT'
        call findgfscript(conv)
	call writeGfscript

 1000	format(a72)

	    
c       WRITE INTO FILES

        open(unit=21,file='par.dat',status='unknown')
        write(21,'(e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,
     .  e10.3)') ep0,ep2,ep3,ep4,ep5,ep_wm
	write(21,'(f10.3,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5,
     .  1x,f10.5,1x,f10.5)') b1,b2,db2,b3_n,db3_n,b3_p,db3_p
        write(21,'(f10.5,1x,f10.5,1x,e10.4,1x,e10.4,1x,f10.5)')
     .  b4,db4,wm_l,wm_r,asym_p
	write(21,'(f10.5,1x,f10.5,1x,f10.5,1x,f20.10,1x,i5,1x,i5,1x,i5)')
     .	x,dop,ec,ei,init,idos,flag
	write(21,'(e15.6,1x,f10.6)')  temp,frac
	write(21,'(f10.6,1x,f10.6)')  dU,Ufact
	write(21,*) 'ep0,ep2,ep3,ep4,ep5,ep_wm'
	write(21,*) 'b1,b2,db2,b3_n,db3_n,b3_p,db3_p'
	write(21,*) 'b4,db4,wm_l,wm_r,asym_p'
	write(21,*) 'x,dop,ec,ei,init,idos,flag'
	write(21,*) 'temp,frac'
	write(21,*) 'dU,Ufact'
        close(21)

	if(temp.eq.0.d0) then
	  conv=dabs(w_m-pw_m)/w_m
	  pw_m=w_m
	else
	  r1=0.d0
	  r2=0.d0
	  do i=-N,N
	    r1 = r1 + dabs(-dimag(Gf(i))+dimag(pGf(i)))*dw(i)
	    r2 = r2 + (-dimag(Gf(i)))*dw(i)
	  end do
	  if(r2.le.0.d0) then
	    write(6,*) 'negative spectral density'
	    stop
	  end if
	  conv=r1/r2
	end if


	write(6,210) nloop,conv,w_m,asym_p
	open(unit=10,file='nloop.dat',status='unknown')
	write(10,210) nloop,conv,w_m,asym_p
	close(10)

	call writeoutput

	nloop=nloop+1

	end do
c	LOOPING ENDS HERE  OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
	write(6,*) 'NOW TRANSPORT PROPERTIES'
c        call interpol2(1)
	call dinftransp(1)

        call intener ()   
 202    format(e16.9,1x,e16.9,1x,e16.9,1x,e16.9)
 210	format('nloop=',i4,'  conv=',f10.6,'  w_m=',e16.9,'  
     .	eta=',f10.6)
 207	format('# nloop=',i4)

	stop
	end

c	************************************************************

        subroutine Grid_Gf(flag)
        integer i1,i2,i3,flag
        include 'Glob_decl'
        real*8 nb_l,nb_r,pb_l,pb_r,yy(-Nm:Nm)
        real*8 l1,l2

        include 'Glob_cons'

c       W>0; Gf pole

        i=N
        r3=-dimag(Gf(N))
        do while(w(i).gt.b1)
           r4=-dimag(Gf(i))
           if(r4.ge.r3) then
              r3=r4
              i1=i
           end if
           yy(i)=-dimag(Gf(i))/pi
           i=i-1
        end do

        i=i1-5
        do while(yy(i).gt.0.1d0.and.yy(i).gt.yy(i-1))
           i=i-1
        end do
        i2=i
        pb_l=w(i2)

        i=i1+5
        do while(yy(i).gt.0.1d0.and.yy(i).gt.yy(i+1))
           i=i+1
        end do
        i3=i
        pb_r=w(i3)

c       W>0; Gf pole

        i=-N
        r3=-dimag(Gf(-N))
        do while(w(i).lt.-b1)
           r4=-dimag(Gf(i))
           if(r4.ge.r3) then
              r3=r4
              i1=i
           end if
           yy(i)=-dimag(Gf(i))/pi
           i=i+1
        end do

        i=i1-5
        do while(yy(i).gt.0.01d0.and.yy(i).gt.yy(i-1))
           i=i-1
        end do
        i2=i
        nb_r=-w(i2)

        i=i1+5
        do while(yy(i).gt.0.01d0.and.yy(i).gt.yy(i+1))
           i=i+1
        end do
        i3=i
        nb_l=-w(i3)

        write(6,*) 'Poles'
        write(6,*) 'Up',nb_l,nb_r,' Dn',pb_l,pb_r


        if(flag.eq.1) then

          l1=min(b3_n-db3_n,nb_l)
          l2=max(b3_n+db3_n,nb_r)

          r1=min(b3_p-db3_p,pb_l)
          r2=max(b3_p+db3_p,pb_r)

        else

          l1=nb_l
          l2=nb_r

          r1=pb_l
          r2=pb_r

        end if


        b3_n=0.5d0*(l1+l2)
        db3_n=0.5d0*(l2-l1)

        b3_p=0.5d0*(r1+r2)
        db3_p=0.5d0*(r2-r1)

        write(6,*) 'n-b3,db3  ',b3_n,db3_n
        write(6,*) 'p-b3,db3  ',b3_p,db3_p

        return
        end

c       ************************************************************

        subroutine Grid_Gfscript
        integer i1,i2,i3
        include 'Glob_decl'
        real*8 nb_l,nb_r,pb_l,pb_r,yy(-Nm:Nm)
        real*8 l1,l2

        include 'Glob_cons'

c       W>0; Gfscript pole

        i=N
        r3=-dimag(Gfscript(1,N))
        do while(w(i).gt.b1)
           r4=-dimag(Gfscript(1,i))
           if(r4.ge.r3) then
              r3=r4
              i1=i
           end if
           yy(i)=-dimag(Gfscript(1,i))/pi
           i=i-1
        end do

        i=i1-5
        do while(yy(i).gt.0.01d0.and.yy(i).gt.yy(i-1))
           i=i-1
        end do
        i2=i
        pb_l=w(i2)

        i=i1+5
        do while(yy(i).gt.0.01d0.and.yy(i).gt.yy(i+1))
           i=i+1
        end do
        i3=i
        pb_r=w(i3)

        r1=pb_l
        r2=pb_r

        b3_p=0.5d0*(r1+r2)
        db3_p=0.5d0*(r2-r1)


c       W<0; Gfscript pole

        i=-N
        r3=-dimag(Gfscript(2,-N))
        do while(w(i).lt.-b1)
           r4=-dimag(Gfscript(2,i))
           if(r4.ge.r3) then
              r3=r4
              i1=i
           end if
           yy(i)=-dimag(Gfscript(2,i))/pi
           i=i+1
        end do

        i=i1-5
        do while(yy(i).gt.0.01d0.and.yy(i).gt.yy(i-1))
           i=i-1
        end do
        i2=i
        nb_r=-w(i2)

        i=i1+5
        do while(yy(i).gt.0.01d0.and.yy(i).gt.yy(i+1))
           i=i+1
        end do
        i3=i
        nb_l=-w(i3)

        write(6,*) 'Poles'
        write(6,*) 'Up',nb_l,nb_r,' Dn',pb_l,pb_r


        l1=nb_l
        l2=nb_r

        b3_n=0.5d0*(l1+l2)
        db3_n=0.5d0*(l2-l1)


	if(Ap(1).gt.0.d0) then
	   b3_p=omp(1)
	   db3_p=0.2d0
	end if
	if(Ap(2).gt.0.d0) then
	   b3_n=-omp(2)
	   db3_n=0.2d0
	end if

        write(6,*) 'n-b3,db3  ',b3_n,db3_n
        write(6,*) 'p-b3,db3  ',b3_p,db3_p

        return
        end


c       ************************************************************

        subroutine Grid_Sig
        integer i1,i2,i3
        include 'Glob_decl'

        include 'Glob_cons'

	r1=-dimag(ssig(-N))
	do i=-N+1,0
	  r2=-dimag(ssig(i))
	  if(r2.gt.r1) then
	     r1=r2
	     i1=i
	  end if
	end do
	b2=dabs(w(i1))
	db2=0.05d0

        return
        end

c       ************************************************************



	subroutine selfenergy(nloop,init,conv)

	include 'Glob_decl'
     	real*8 UU,evalSg,fact,sgn,conv
	integer imax,nloop,init,findU_crude,flag

	external evalsg,sgn
	external locate,findU_crude

	include'Glob_cons'	

	write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

	if(temp.eq.0.d0) then	!  symmetry restoration for T=0 only

	write(6,*) 'FINDING U FROM SYMMETRY RESTORATION'

	if(nloop.eq.1.and.init.eq.0) then
	  U=1.d0/dreal(pi0(0))-dU
	else if(U.gt.1.d0/dreal(pi0(0))) then
	  U=1.d0/dreal(pi0(0))-dU
	end if
	write(6,*) 'U=',U

	call grid_tpi
	write(6,*) 'w_m=',w_m

	fact=Ufact
	if(conv.gt.1.0d0) then
	flag=findU_crude(fact)
	do while(flag.eq.1)
	   write(6,*) '___________________________________'
	   write(6,*) 'Initial Guess for U not good'
	   U=U-fact
	   flag=findU_crude(fact)
	end do
	end if

	if(conv.gt.1.0d0) then
	  call grid_tpi
	  flag=findU_crude(fact)
	end if

	open(unit=33,file='tpi.dat',status='unknown')
	do i=-N,N
          z1=dcmplx(dreal(pi0(i)),sgn(w(i))*dimag(pi0(i)))
	  tpi(i)=z1/(1.d0-U*z1)  
	  write(33,*) w(i),dimag(tpi(i)),dreal(tpi(i))
	end do
	close(33)

	call grid_tpi
	call findU(fact)	! Find U from Sym. Rest.

	end if		!  T=0

	open(unit=33,file='tpi.dat',status='unknown')
	do i=-N,N
          z1=dcmplx(dreal(pi0(i)),sgn(w(i))*dimag(pi0(i)))
	  tpi(i)=z1/(1.d0-U*z1)  
	  write(33,*) w(i),dimag(tpi(i)),dreal(tpi(i))
	end do
	close(33)

	r3=dabs(dimag(tpi(-N)))
	do i=-N+1,N
	   if(dabs(dimag(tpi(i))).gt.r3) then
	     r3=dabs(dimag(tpi(i)))
	     imax=i
	   end if
	end do

	w_m=w(imax)
	
	write(6,*) 'EVALUATING SIGMA WITH THE ABOVE U'
	UU=U
	r1= evalSg(UU,1,0)     ! With this U, find the full Self Energy

	
	write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

	return
	end
c	************************************************************

	subroutine grid_tpi
	include 'Glob_decl'
	integer ilf,irt,imax,M1
	real*8 ep,sgn

        external sgn


	include'Glob_cons'	

	do i=-N,N
          z1=dcmplx(dreal(pi0(i)),sgn(w(i))*dimag(pi0(i)))
	  tpi(i)=z1/(1.d0-U*z1)  
	end do

	r3=dabs(dimag(tpi(-N)))
	do i=-N+1,N
	   if(dabs(dimag(tpi(i))).gt.r3) then
	     r3=dabs(dimag(tpi(i)))
	     imax=i
	   end if
	end do

	w_m=w(imax)
	write(6,*) 'w_m=',w_m

        i=imax
        do while(dabs(dimag(tpi(i))).gt.10.d0.or.w(i).gt.w_m/10.d0)
          i=i-1
        end do
        ilf=i+1
	wm_l=min(w(ilf),b1-0.0002d0)

        i=imax
        do while(dabs(dimag(tpi(i))).gt.15.d0.or.w(i).lt.w_m*3.d0)
          i=i+1
        end do
        irt=i-1
	wm_r=min(w(irt),b1-0.0002d0)

	
        ep=min(dw(ilf),ep_wm*(wm_r-wm_l))
        M1=nint((wm_r-wm_l)/ep)
        do while(M1.gt.13000)
          ilf=ilf+1
          irt=irt-1
          wm_l=w(ilf)
          wm_r=w(irt)
          ep=min(dw(ilf),ep_wm*(wm_r-wm_l))
          M1=nint((wm_r-wm_l)/ep)
        end do

	call interpol(0)

	return
	end
c	************************************************************

	subroutine evalpi0(flag)

	include 'Glob_decl'
	real*8 rho1(-Nm:Nm),rho2(-Nm:Nm),impi0(-Nm:Nm),solveq
	real*8 lb4,rb4
	integer i1,i2,i3,N1,locate,lj,flag
	real*8 wl(2*Nm),wi,sgn
	complex*16 pi0l(2*Nm)

	external solveq,locate,sgn
	include 'Glob_cons'


	write(6,*) 'flag=',flag
	if(flag.ne.1.and.temp.eq.0.d0) then

	  do i=-N,N
	    rho1(i)=-dimag(gfscript(1,i))/pi
	    rho2(i)=-dimag(gfscript(2,i))/pi
	  end do

	  call imagpi0_zt(rho1,rho2,impi0)


	else if(flag.eq.1) then

	  i=1
	  open(unit=10,file='Pi0.dat',status='unknown')
 7	  read(10,*,end=8) r1,r2,r3
 	  wl(i)=r1
	  pi0l(i)=dcmplx(r3,r2)
	  i=i+1
	  goto 7
 8	  close(10)
 	  N1=i-1

	  Npi0=N1/2
	  do i=-Npi0,Npi0
	    w_pi(i)=wl(i+Npi0+1)
	    Pi0_G(i)=pi0l(i+Npi0+1)
	  end do

	  do i=-N,N
	    wi=w(i)
	    j=locate(N1,wl,1,N1,wi)
	    pi0(i)=pi0l(j)+(pi0l(j+1)-pi0l(j))/(wl(j+1)-wl(j))*
     .		   (wi-wl(j))
     	    impi0(i)=dimag(pi0(i))
	  end do

	end if

	if(temp.eq.0.d0.or.flag.eq.1) then

	r3=dabs(impi0(-N))
	do i=-N,N
	  if(dabs(impi0(i)).gt.r3) then
	    r3=dabs(impi0(i))
	    i1=i
	  end if
	end do

	i=i1
	do while(dabs(impi0(i)).gt.0.1d0.and.
     .        dabs(impi0(i)).ge.dabs(impi0(i-1)))
	  i=i-1
	end do
	i2=i

	i=i1
	do while(dabs(impi0(i)).gt.0.1d0.and.
     .        dabs(impi0(i)).ge.dabs(impi0(i+1)))
	  i=i+1
	end do
	i3=i

	lb4=w(i2)
	rb4=w(i3)
	b4=0.5d0*(lb4+rb4)
	db4=b4-lb4

	write(6,*) 'b4,db4  ',b4,db4

	call interpol(1)

	end if
	
 20	if(flag.ne.1) then

	  do i=-N,N
	    rho1(i)=-dimag(gfscript(1,i))/pi
	    rho2(i)=-dimag(gfscript(2,i))/pi
	  end do


c	  Convolution of down and up spec. functions
	  call convolpi(rho1,rho2)

	  do i=-N,N
	     w_pi(i)=w(i)
	     Pi0_G(i)=pi0(i)
	  end do
	  Npi0=N

	end if

	return
	end


c	************************************************************
	
	real*8 function evalsg(UU,flag,flag1)

	include 'Glob_decl'
	integer flag,ilf,irt,spin,sp,spbar,imax,flag1
	real*8 rho1(-Nm:Nm),rho2(-Nm:Nm),tmax,sgn,UU,spj
     	complex*16 tmpsig(-Nm:Nm),sig0,sig1,sig2,
     .	y1(-Nm:Nm),y2(-Nm:Nm)
     	real*8 wl(2*Nm),kktsing,npc,pc
	integer llim,ulim,locate,lj

	external sgn,locate,kktsing
	include 'Glob_cons'

c	This function evaluates self energy at
c	w=0   for flag.eq.0
c	All w for flag.ne.0
c	flag1=1	  Change Grid at w_m
c	     =0	  Grid as it is
c	Returns the symmetry restoration condition value.

	U=UU	

	if(flag1.eq.1) call grid_tpi

	do i=-N,N
          z1=dcmplx(dreal(pi0(i)),sgn(w(i))*dimag(pi0(i)))
	  tpi(i)=z1/(1.d0-U*z1)  
	end do

	r3=dabs(dimag(tpi(-N)))
	do i=-N+1,N
	   if(dabs(dimag(tpi(i))).gt.r3) then
	     r3=dabs(dimag(tpi(i)))
	     imax=i
	   end if
	end do

	w_m=w(imax)

	do i=-N,N
	   wl(i+N+1)=w(i)
	end do

        if(flag.eq.2) then
        do i=-N,N
           z3=dcmplx(dreal(Gfsc_lp(1,i)),sgn(w(i))*
     .     dimag(Gfsc_lp(1,i)))
           y1(i)=tpi(i)*z3
           z3=dcmplx(dreal(Gfsc_lp(2,i)),sgn(w(i))*
     .     dimag(Gfsc_lp(2,i)))
           y2(i)=z3*tpi(-i)
        end do

        z1=zero
        z2=zero
        do i=-N,N
           z1=z1+dw(i)*y1(i)
           z2=z2+dw(i)*y2(i)
        end do

        sig1=U**2*(z1)/(2.d0*pi*ii)     ! up
        sig2=U**2*(z2)/(2.d0*pi*ii)     ! down

c          Large Pole contribution
           do i=-N,N
              rho1(i)=dimag(tpi(i))
           end do

           if(Ap(2).gt.0.d0) then
           sp=1
           spbar=2
	   r2=dabs(omp(spbar))
           r1=kktsing(rho1,r2)
           sig2=sig2-U**2*Ap(2)*r1/pi
           end if

           if(Ap(1).gt.0.d0) then
           sp=2
           spbar=1
	   r2=dabs(omp(spbar))
           r1=kktsing(rho1,r2)
           sig1=sig1+U**2*Ap(1)*r1/pi
           end if


	sig0=sig1-sig2

	evalsg=dreal(sig0)-U*mu0

	return
	   
	end if
	  

	do i=-N,N
	   rho1(i)=dimag(tpi(i))	! Convolsg_zt
	   rho2(i)=-dimag(Gfsc_lp(1,i))/pi
	end do
	
	spin=1
	if(temp.eq.0.d0) then
	  call convolsg_zt(rho1,rho2,tmpsig,flag,spin)
	else
	  call convolsg_ft(rho1,rho2,tmpsig,flag,spin)
	end if

	if(flag.eq.0) then	! Zero frequency value
	  ulim=0
	  llim=0
	else			! flag=1 Full self energy
	  ulim=N
	  llim=-N
	end if

	do i=llim,ulim
	   sigma(2,i)=U**2*tmpsig(i)
	end do

	spin=-1
	do i=-N,N
	   rho1(i)=dimag(tpi(-i))
	   rho2(i)=-dimag(Gfsc_lp(2,i))/pi
	end do

	if(temp.eq.0.d0) then
	  call convolsg_zt(rho1,rho2,tmpsig,flag,spin)
	else
	  call convolsg_ft(rho1,rho2,tmpsig,flag,spin)
	end if

	do i=llim,ulim
	   sigma(1,i)=U**2*tmpsig(i)
	end do

	if(flag.eq.1) then
	write(6,*) 'diffsig',sigma(2,0)-sigma(1,0)-U*mu0
	end if

	evalsg=dreal(sigma(2,0)-sigma(1,0))-U*mu0

 299	return
	end

c	************************************************************

	subroutine convolpi(rho1,rho2)

	include 'Glob_decl'
	integer ljvec(-Nm:Nm),llim,ulim,lj,flag,locate
	real*8 wji,sgn,tmax
	real*8 rho1(-Nm:Nm),rho2(-Nm:Nm),d2(-Nm:Nm),rho2ji,
     .	       impi0(-Nm:Nm),repi0(-Nm:Nm)
	real time1,time2,time3

	external sgn,locate
	include 'Glob_cons'

	if(temp.eq.0.d0) then
	   call imagpi0_ft(rho1,rho2,impi0)
	else
	   call imagpi0_ft(rho1,rho2,impi0)
	end if

	do i=-N,N
	   d2(i)=-impi0(i)/pi
	end do

	time1=secnds(0.0)
	
	flag=1
	call kktransf(d2,repi0,flag)

	time2=secnds(time1)

	write(6,*) 'kktransf',time2

	do i=-N,N
	   pi0(i)=repi0(i)+ii*impi0(i) 
	end do

	return
	end


c	************************************************************


	subroutine imagpi0_zt(rho1,rho2,impi0)

	include 'Glob_decl'
	integer ljvec(-Nm:Nm),llim,ulim,lj,flag,locate
	real*8 wji,sgn
	real*8 rho1(-Nm:Nm),rho2(-Nm:Nm),d2(-Nm:Nm),rho2ji,
     .	       impi0(-Nm:Nm),wl(2*Nm)
	real time1,time2,time3

	external sgn,locate
	include 'Glob_cons'

	time1=secnds(0.0)

	do i=-N,N-1
	    wl(i+N+1)=w(i)
    	    d2(i)=(rho2(i+1)-rho2(i))/(w(i+1)-w(i))
	end do
	wl(2*N+1)=w(N)
	d2(N)=d2(N-1)

	do 98 i=-N,-1
	  llim=0
	  do j=i+1,0
	    wji=w(j)-w(i)
	    ulim=min(llim+1,N)
	    do while(wji.ge.w(ulim).and.ulim.ne.N)
	      llim=llim+1
	      ulim=min(llim+1,N)
	    end do
	    ljvec(j)=llim
	  end do
	  r1=0.d0
	  do j=i+1,0
	    wji=w(j)-w(i)
	    lj=ljvec(j)
	    llim=lj
c	    This implies that  ==> w(lj)<= w_j-w_i <=w(lj+1)
	    rho2ji=rho2(lj)+(rho2(lj+1)-rho2(lj))*
     .	    ((wji-w(lj))/(w(lj+1)-w(lj)))
	    r1=r1+rho1(j)*rho2ji*dw(j)
	  end do
	  impi0(i)=-pi*r1 
 98	continue

	do 99 i=1,N
	  r1=0.d0
	  llim=-i
	  do j=1,i
	    wji=w(j)-w(i)
	    do while(wji.ge.w(llim+1))
	      llim=llim+1
	    end do
	    ljvec(j)=llim
	  end do
	  do j=1,i
	    wji=w(j)-w(i)
	    lj=ljvec(j)
	    rho2ji=rho2(lj)+(rho2(lj+1)-rho2(lj))*
     .	    ((wji-w(lj))/(w(lj+1)-w(lj)))
	    r1=r1+rho1(j)*rho2ji*dw(j)
	  end do
	  impi0(i)=pi*r1 
 99	continue

 	do i=-N,N
	  impi0(i)=impi0(i)+pi*0.5d0*(dw(0)*rho1(0)*
     .		   rho2(-i) - dw(i)*rho1(i)*rho2(0))
	end do

 	impi0(0)=0.d0


	time2=secnds(time1)
	write(6,*) 'convol',time2

	return
	end
c	**************************************************************

	subroutine convolsg_zt(rho1,rho2,totsig,flag,spin)

	include 'Glob_decl'
	integer ljvec(-Nm:Nm),llim,ulim,lj,flag,spin,locate,
     .	spin1,li,ui,sp,spbar
	real*8 wi,wji,sgn,sign,piji,spi
	real*8 rho1(-Nm:Nm),rho2(-Nm:Nm),d2(-Nm:Nm),rho2ji,
     .	       imsig(-Nm:Nm),wl(2*Nm),resig(-Nm:Nm),kktsing,
     .		rho11(-Nm:Nm)
     	complex*16 totsig(-Nm:Nm)

	external sgn,locate,kktsing
	include 'Glob_cons'

	do i=-N,N-1
	    wl(i+N+1)=w(i)
    	    d2(i)=(rho2(i+1)-rho2(i))/(w(i+1)-w(i))
	end do
	wl(N+N+1)=w(N)
	d2(N)=d2(N-1)

	do 198 i=-N,-1
	  llim=i
	  do j=1,-i
	    wji=w(j)+w(i)
	    do while(wji.ge.w(llim+1))
	      llim=llim+1
	    end do
	    ljvec(j)=llim
	  end do
	  r1=0.d0
	  do j=1,-i
	    wji=w(j)+w(i)
	    lj=ljvec(j)
c	    This implies that  ==> w(lj)<= w_j+w_i <=w(lj+1)
	    rho2ji=rho2(lj)+(rho2(lj+1)-rho2(lj))*((wji-w(lj))/
     .                                          (w(lj+1)-w(lj)))
	    r1=r1+rho1(j)*rho2ji*dw(j)
	  end do
	  imsig(i)=-r1+0.5d0*dw(i)*rho1(-i)*rho2(0)
 198	continue

	do 199 i=1,N
	  llim=0
	  do j=-i+1,0
	    wji=w(j)+w(i)
	    ulim=min(llim+1,N)
	    do while(wji.ge.w(ulim).and.ulim.ne.N)
	      llim=llim+1
	      ulim=min(llim+1,N)
	    end do
	    ljvec(j)=llim
	  end do
	  r1=0.d0
	  do j=-i+1,0
	    wji=w(j)+w(i)
	    lj=ljvec(j)
c	    This implies that  ==> w(lj)<= w_j+w_i <=w(lj+1)
	    rho2ji=rho2(lj)+(rho2(lj+1)-rho2(lj))*((wji-w(lj))/
     .                                          (w(lj+1)-w(lj)))
	    r1=r1+rho1(j)*rho2ji*dw(j)
	  end do
	  imsig(i)=-r1-0.5d0*dw(i)*rho1(-i)*rho2(0)
 199	continue

 	imsig(0)=0.d0

	do i=-N,N
	   resig(i)=0.d0
	   d2(i)=-imsig(i)/pi
	end do

	call kktransf(d2,resig,flag)

	if(flag.eq.0) then
	  llim=0
	  ulim=0
	else
	  llim=-N
	  ulim=N
	end if

c       Large pole contribution

        if(spin.eq.-1) then
	  if(Ap(2).gt.0.d0) then
          li=-N
          ui=1+locate(2*N+1,wl,1,N,omp(2))-N-1
	  end if
          sp=1
          spbar=2
        end if
	if(spin.eq.1) then
	  if(Ap(1).gt.0.d0) then
          li=locate(2*N+1,wl,N,2*N+1,omp(1))-N-1
          ui=N
	  end if
          sp=2
          spbar=1
        end if

c	Imaginary part

        if(Ap(spbar).gt.0.d0) then
        spi=dfloat(spin)
        do i=li,ui
          wi=(omp(spbar)-w(i))
          j=locate(2*N+1,wl,1,2*N+1,wi)-N-1
          r1=rho1(j)+(rho1(j+1)-rho1(j))*((wi-w(j))/
     .	            (w(j+1)-w(j)))
          if(wi.ne.0.d0)imsig(i)=imsig(i)-Ap(spbar)*r1
        end do

c	Real part

	if(spin.eq.-1) then
	  do i=-N,N
	    rho11(-i)=rho1(i)
	  end do
	  do i=-N,N
	    rho1(i)=rho11(i)
	  end do
	end if

	do i=llim,ulim
	  r2=dabs(omp(spbar))-spi*w(i)
	  r1=kktsing(rho1,r2)
	  resig(i)=resig(i)+spi*Ap(spbar)*r1/pi
	end do

        end if

	do i=-N,N
	   d2(i)=-imsig(i)/pi
	end do

	if(spin.eq.1) then
	  do i=-N,N
	  write(40,*) w(i),d2(i),resig(i)
	  end do
	  close(40)
	else
	  do i=-N,N
	  write(41,*) w(i),d2(i),resig(i)
	  end do
	  close(41)
	end if


	do i=llim,ulim
	   totsig(i)=resig(i)+ii*imsig(i)   ! Retarded Sigma
	end do

	return
	end

c	**************************************************************

	real*8 function findgcf(ep_f)
	include 'Glob_decl'
	integer spj,iter
	real*8 sp,ep_f,sgn,conv,nf,nd,xx,yy,a1,a2,alpha,abeta,
     .	sig0,ef_s,slope,Lutt_val,lval
	complex*16 num,den,hilbtr,phyb,nhyb,gup,gdn,gamma,g0,
     .  sigcpa(-N:N),ssigc(-N:N),ssigf(-N:N),delta(-N:N)
	complex*16 gdd(-N:N),delta0(-N:N),gcinit(-N:N),gauss
        external hilbtr,Lutt_val,sgn,gauss
	
	include 'Glob_cons'


	dfac=0.5d0

	x0=0.5d0*U*mu0
	ei0=ep_f+U*n0/2.d0
        asym_p=1.d0+2.d0*ep_f/U
        write(6,*)'asym_p=',asym_p

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
           do i=-N,N
           z1=w(i)+ii*1.D-10
           delta0(i)=V**2/(z1-ec-hyb(i))
           write(51,*)w(i),-dimag(delta0(i))/pi
           enddo
           close(51)


	r2=0.d0
	r3=0.d0
	do i=-N,N
           z1=w(i)+ii*1.D-10
           iter=1
           conv=1.d0
	   phyb=hyb(i)
	   z3=sigma(1,i)
	   gdn=z1-ec-V**2/(z1-ei0-x0-z3)
	   z3=sigma(2,i)
	   gup=z1-ec-V**2/(z1-ei0+x0-z3)
           g0=z1-ec
           do while(conv.ge.1.D-10.and.iter.le.100000)
!	    gamma=(2.d0*gup*gdn-hyb(i)*(gup+gdn))/
!     .            (gup+gdn-2.d0*hyb(i))
            
           gamma=hyb(i)+(2.d0*(g0-hyb(i))*(gup*gdn 
     .              -hyb(i)*(gup+gdn)+hyb(i)**2))/  
     .            ((1.d0-dop)*(g0-hyb(i))* 
     .             (gup+gdn-2.d0*hyb(i))+2.d0*dop* 
     .           (gup*gdn-hyb(i)*(gup+gdn)+hyb(i)**2))
             if(zabs(gamma).le.1.D3) then
               nhyb=(gamma-1.d0/hilbtr(gamma,idos))
            else
               nhyb=alpha/gamma+abeta/gamma**3  ! From moment expansion
            end if
            hyb(i)=dfac*phyb+(1.d0-dfac)*nhyb
            conv=zabs(hyb(i)-phyb)
            phyb=hyb(i)
            if(iter.eq.100000) write(6,*) 'MAX LIM',i,conv
            iter=iter+1
          end do
          sigcpa(i)=gamma
          delta(i)=V**2/(z1-ec-hyb(i))
          Gccpa(i)=1.d0/(gamma-hyb(i))
           Gc(i)=0.5d0*(1.d0/(gup-hyb(i))+1.d0/(gdn-hyb(i))) 
          ssigf(i)=z1-ep_f-V**2/(z1-ec-gamma)
          ssig(i)=z1-ep_f-V**2/(z1-ec-hyb(i)-1.d0/Gc(i))
          Gf(i)=1.d0/(z1-ep_f-ssig(i)-V**2/(z1-ec-hyb(i)))
          ssigc(i)=V**2/(z1-ep_f-ssig(i))
          Gfcpa(i)=(1.d0-dop)/(z1-ep_f-ssig(i)
     .                 -V**2/(z1-ec-hyb(i)))
          r2=r2+2.d0*dw(i)*(-dimag(gf(i))/pi)*ferm(i)
          r3=r3+2.d0*dw(i)*(-dimag(gc(i))/pi)*ferm(i)
        end do
        ef_s=ep_f+dreal(ssig(0))
      
        !ef_s=ep_f+dreal(ssig(0))+dop*dreal(V**2/(-ec-hyb(0)))
         nf=r2
        nd=r3
        write(6,*) 'nf=',nf,'  nd=',nd
        

	open(unit=10,file='hyb.dat',status='unknown')
        do i=-N,N
        write(10,*) w(i),-dimag(hyb(i)),dreal(hyb(i))
        end do
        close(10)
        open(unit=11,file='Gc_cpa.dat',status='unknown')
        do i=-N,N
        write(11,*)w(i),-dimag(Gccpa(i)),dreal(Gccpa(i))
        enddo
        close(11)
        open(unit=12,file='Gf_cpa.dat',status='unknown') 
        do i=-N,N
        write(12,*)w(i),-dimag(Gfcpa(i)),dreal(Gfcpa(i))
        enddo
        close(12)
        !open(unit=13,file='sigcpa.dat',status='unknown')
        !do i=-N,N
        !write(13,*)w(i),dimag(sigcpa(i)),dreal(sigcpa(i))
        !enddo
        !close(13)
        open(unit=14,file='sigc.dat',status='unknown')
        do i=-N,N
        write(14,*)w(i),dimag(ssigc(i)),dreal(ssigc(i))
        enddo
        close(14)
        open(unit=15,file='sigcpaf.dat',status='unknown')        
        do i=-N,N
        write(15,*)w(i),dimag(ssigf(i)),dreal(ssigf(i))
        enddo
        close(15)
	r1=0.d0
	r2=0.d0
	r3=0.d0
	r4=0.d0
	do i=-100,100
	  xx=w(i)
	  yy=dreal(ssig(i))
	  r1=r1+xx
	  r2=r2+xx**2
	  r3=r3+yy
	  r4=r4+xx*yy
	end do

	slope=(r1*r3-201.d0*r4)/(r2-201.d0*r2)

	Zfac=1.d0/(1.d0-slope)

	open(unit=25,file='outpar.dat',status='unknown')
	write(25,500) x,ei,ec,U
	write(25,501) x0,ei0,n0,mu0
	write(25,502) nf,nd,asym_p
	write(25,503) ep_f,ef_s,Zfac,w_m
	write(25,504) V**2,temp,beta,frac
	close(25)
 500	format('x=',f6.3,' ei=',f16.10,' ec=',f6.3,' U=',f20.12)
 501	format('x0=',f12.6,' ei0=',f12.6,' n0=',f12.6,' mu0=',f12.6)
 502	format(' nf=',f12.6,' nd=',f12.6,' asym_p=',f12.6)
 503	format('ep_f=',f10.6,' ef_s=',f10.6,' Z=',f10.6,' w_m=',e20.8)
 504	format('V2=',f10.6,' T=',e15.6,'  beta=',e15.6,' frac=',f10.3)
        open(unit=13,file='sigcpa.dat',status='unknown')
        do i=-N,N
        if(w(i).gt.(0.001d0*V**2*Zfac).and.w(i).le.0.1*V**2*Zfac)then  
        write(13,*)w(i)/(V**2*Zfac),dimag(sigcpa(i)),dreal(sigcpa(i))
        endif
        enddo
        close(13)
        do i=-N,N
        z1=w(i)+ii*1.D-10
        gdd(i)=V**2*(1.d0-dop)*Zfac/(z1-ep_s*Zfac-Zfac*delta0(i)*dop)
        end do
        do i=-N,N
        write(42,*)w(i),-dimag(gdd(i))
        enddo
        close(42)

	lval=Lutt_val()

	write(6,600) ep_f,asym_p,lval

 600	format('ep_f=',f12.7,'eta=',f6.3,'lval=',f12.7)

	findgcf=lval
 
	return
	end


c	**************************************************************

	complex*16 function hilbtr(z,idos)
	complex*16 z,gauss,sem
	integer idos

	if(idos.eq.1) then
  	  hilbtr=gauss(z)
	else
	  hilbtr=sem(z)
	end if

	return
	end

c	**************************************************************

	real*8 function Lutt_val()
	include 'Glob_decl'
	
	real*8 dwi
	complex*16 dsig
	include 'Glob_cons'


        r1=0.d0
        i=-N+1
        do while(w(i).lt.0.d0)
          dwi=0.5d0*(w(i+1)-w(i-1))
          dsig=0.5d0*((ssig(i)-ssig(i-1))/(w(i)-w(i-1))+
     .    (ssig(i+1)-ssig(i))/(w(i+1)-w(i)))
          r1=r1+dimag(Gf(i)*dsig)*dwi
          i=i+1
        end do


	Lutt_val=r1
c        write(6,*) 'Lutt thm =',r1


	
	return
	end


c	************************************************************

	subroutine findgfscript(conv)
	include 'Glob_decl'
	integer spj
	real*8 sp,conv
	complex*16 hybr

	include 'Glob_cons'


        do i=-N,N
           z1=w(i)+ii*1.D-10
	   hybr=V**2/(z1-ec-hyb(i))
           do spj=1,2
             sp=2.d0*(dfloat(spj)-1.5d0)
             z2=z1-ei+sp*x-hybr
	     gfscript(spj,i)=1.d0/z2
           end do
        end do

	!call takeoutpole(0)
	!call Grid_Gfscript
	!call interpol(1)

	r1=0.d0
        r2=0.d0
        do i=-N,N
           z1=w(i)+ii*1.D-10
	   hybr=V**2/(z1-ec-hyb(i))
           do spj=1,2
             sp=2.d0*(dfloat(spj)-1.5d0)
             z2=z1-ei+sp*x-hybr
	     gfscript(spj,i)=1.d0/z2
             r1=r1+dw(i)*(-dimag(gfscript(spj,i))/pi)
             r2=r2+dw(i)*(-dimag(gfscript(spj,i))/pi)*ferm(i)
           end do
        end do

	write(6,*) 'Gfscript Normalisation=',r1/2.d0
	write(6,*) 'n0=',r2

	
	return
	end


c	************************************************************

	subroutine findepf(ep_f)
	real*8 ep_f,fact,findgcf,lepf,repf,mepf,sign
	integer iter,flag
	include 'Glob_decl'

	external findgcf

           fact=U*0.002d0

	   write(6,*) 'Guess ep_f'
           write(6,*)'ep_f=',ep_f

	   write(6,599) 
	   sign=1.d0

	   r1=findgcf(ep_f)
	   lepf=ep_f
	   repf=ep_f
	   r2=r1
	   r3=r1
           write(6,*) r1, r2, r3
	   do while(r1*r2.gt.0.d0.and.r1*r3.gt.0.d0)
	    lepf=lepf-sign*fact
	    r2=findgcf(lepf)
	    repf=repf+sign*fact
	    r3=findgcf(repf)
	   end do

	   if(r1*r2.lt.0.d0) then
	     ep_f=lepf
	   else
	     ep_f=repf
	   end if

	   lepf=ep_f
	   r1=findgcf(lepf)
	   write(6,600) lepf,asym_p,r1
	   repf=lepf+fact*sign
	   r2=findgcf(repf)
	   write(6,600) repf,asym_p,r2
           if(dabs(r1).le.1.D-6.or.dabs(r2).le.1.D-6) then
             if(dabs(r1).le.1.D-6) then
               mepf=lepf
               goto 601
             else if (dabs(r2).le.1.D-6) then
               mepf=repf
               goto 601
             end if
           end if

	   if((r1*r2.ge.0.d0).and.(dabs(r1).le.dabs(r2))) then
	      sign=-sign
	      r3=r1
	      r1=r2
	      r2=r3
	      mepf=lepf
	      lepf=repf
	      repf=mepf
	   end if
	   iter=1
	   do while(r1*r2.ge.0.d0)
	       r1=r2
	       repf=repf+sign*fact
	       r2=findgcf(repf)
	       write(6,600) repf,asym_p,r2
	       iter=iter+1
	   end do
 599	   format("       ep_f             r2  ")
 600	   format('ep_f=',f12.7,'eta=',f6.3,'lval=',f12.7)
 	
	   write(6,*) 'now find the actual ep_f'
	   i=1
	   lepf=repf-sign*fact
           if(dabs(r1).le.1.D-6.or.dabs(r2).le.1.D-6) then
             if(dabs(r1).le.1.D-6) then
               mepf=lepf
               goto 601
             else if (dabs(r2).le.1.D-6) then
               mepf=repf
               goto 601
             end if
           end if

	   do while(dabs(r1).gt.1.D-6.and.dabs(r2).gt.1.D-6.
     .		   and.i.le.25)
	      mepf=lepf-r1*(repf-lepf)/(r2-r1)
	      r3=findgcf(mepf)
	      write(6,600) mepf,asym_p,r3
	      if(r3*r1.ge.0.d0) then
	        lepf=mepf
		r1=r3
	      else
	        repf=mepf
		r2=r3
	      end if
	      i=i+1
	   end do

 601       asym_p=1.d0+2.d0*mepf/U
	   write(6,*) 'FINAL asym_p=',asym_p,'  ep_f=',mepf
	   write(6,*) '________________________________________'

	return
	end


c	************************************************************
	
	subroutine findU(dfact)
	real*8 dfact,fact,evalSg,lU,rU,mU,pU,sign,pw_m,U00
	integer iter,flag
	include 'Glob_decl'

	external evalSg

	   write(6,*) 'Guess U'
	   U00=U
	   write(6,599) 
	   sign=-1.d0
	   fact=dfact
	   lU=U
	   r1= evalSg(lU,0,0)
	   write(6,600) U,r1,w_m
	   rU=lU+fact*sign
	   r2= evalSg(rU,0,0)
	   write(6,600) U,r2,w_m
	   if(w_m.le.0.d0) then
	      write(6,*) 'W_m negative'
	      stop
	   end if
	   if((r1*r2.ge.0.d0).and.(dabs(r1).le.dabs(r2))) then
	      sign=-sign
	      r3=r1
	      r1=r2
	      r2=r3
	      mU=lU
	      lU=rU
	      rU=mU
	   end if
	   iter=1
	   do while(r1*r2.ge.0.d0)
	       r1=r2
	       rU=rU+sign*fact
	       r2= evalSg(rU,0,0)
	       write(6,600) U,r2,w_m
	       if(w_m.le.0.d0) then
	          write(6,*) 'W_m negative'
		  stop
	       end if
	       if(dabs(r2).gt.dabs(r1).and.r1*r2.gt.0.d0) then
		  goto 601
	       end if
	       iter=iter+1
	   end do
 599	   format("       U             r1           w_m")
 600	   format(f12.7,1x,f12.7,1x,f12.8)

 	   goto 602 		! Found U

 601	   write(6,*) 'Changing Sign because of local extremum'
	   iter=1
	   sign=-sign
	   rU=U00
	   do while(r1*r2.ge.0.d0)
	       r1=r2
	       rU=rU+sign*fact
	       r2= evalSg(rU,0,0)
	       write(6,600) U,r2,w_m
	       if(w_m.le.0.d0) then
	          write(6,*) 'W_m negative'
		  stop
	       end if
	       iter=iter+1
	   end do
 	
 602	   write(6,*) 'now find the actual U'
	   iter=1
	   i=1
	   lU=rU-sign*fact
	   pU=U
	   flag=0
	   do while(dabs(r1).ge.1.D-6.and.dabs(r2).ge.1.D-6.
     .		   and.i.le.25.and.flag.eq.0)
              if(iter.eq.4) iter=1
	      mU=lU-r1*(rU-lU)/(r2-r1)
	      r3= evalSg(mU,0,0)
	      write(6,600) U,r3,w_m
	      if(w_m.le.0.d0) then
	          write(6,*) 'W_m negative'
	          stop
	      end if
	      if(r3*r1.ge.0.d0) then
	        lU=mU
		r1=r3
	      else
	        rU=mU
		r2=r3
	      end if
	      i=i+1
	      iter=iter+1
	      if(dabs(pU-U).le.1.D-7) flag=1
	      pU=U
	   end do

	   U=mU
	   write(6,*) 'SR --U=',U,' w_m=',w_m

	return
	end


c	************************************************************

	integer function findU_crude(dfact)
	real*8 dfact,fact,evalSg,lU,rU,mU,pU,sign,pw_m
	integer iter,flag
	include 'Glob_decl'

	external evalSg

	   write(6,*) 'Guess U'
	   write(6,599) 
	   sign=-1.d0
	   fact=dfact
	   lU=U
	   r1= evalSg(lU,2,0)
	   write(6,600) U,r1,w_m
	   if(w_m.le.0.d0) then
	     write(6,*) 'w_m negative'
	     stop
	   end if
	   rU=lU+fact*sign
	   r2= evalSg(rU,2,0)
	   write(6,600) U,r2,w_m
	   if(w_m.le.0.d0) then
	     write(6,*) 'w_m negative'
	     stop
	   end if
	   if((r1*r2.ge.0.d0).and.(dabs(r1).le.dabs(r2))) then
	      sign=-sign
	      r3=r1
	      r1=r2
	      r2=r3
	      mU=lU
	      lU=rU
	      rU=mU
	   end if
	   iter=1
	   do while(r1*r2.ge.0.d0)
	       r1=r2
	       rU=rU+sign*fact
	       r2= evalSg(rU,2,0)
	       write(6,600) U,r2,w_m
	       if(w_m.le.0.d0) then
	         write(6,*) 'w_m negative'
	         stop
	       end if
	       iter=iter+1
	       if(iter.eq.100) then
	         findU_crude=1
		 return
	       end if
	   end do
 599	   format("       U             r1           w_m")
 600	   format(f12.7,1x,f12.7,1x,f12.8)

 	   if(dabs(r1).lt.1.D-6) then
	        lU=rU-sign*fact
	        r1= evalSg(lU,2,0)
	        write(6,600) U,r1,w_m
	   end if
 	
	   write(6,*) 'now find the actual U'
	   iter=1
	   i=1
	   lU=rU-sign*fact
	   pU=U
	   flag=0
	   do while(dabs(r1).ge.1.D-6.and.dabs(r2).ge.1.D-6.
     .		   and.i.le.25.and.flag.eq.0)
              if(iter.eq.4) iter=1
	      mU=lU-r1*(rU-lU)/(r2-r1)
	      r3= evalSg(mU,2,0)
	      write(6,600) U,r3,w_m
	      if(w_m.le.0.d0) then
	         write(6,*) 'w_m negative'
	         stop
	      end if
	      if(r3*r1.ge.0.d0) then
	        lU=mU
		r1=r3
	      else
	        rU=mU
		r2=r3
	      end if
	      i=i+1
	      iter=iter+1
	      if(dabs(pU-U).le.1.D-7) flag=1
	      pU=U
	   end do

	   U=mU
c	   r1=evalSg(mu,2,1)
	   write(6,*) 'SR --U=',U,' w_m=',w_m
	   findU_crude=0

	return
	end


c	************************************************************
	
	subroutine interpol2(flag)
	include 'Glob_decl'
	real*8 wl(2*Nm),wi,wl_p(2*Nm)
	complex*16 pi0l(2*Nm),gcl(2*Nm),gfl(2*Nm)
	integer sp,locate,N1,spj,flag,N2,j1
        complex*16 sigmal(2,2*Nm),gfscl_lp(2,2*Nm),sigtl(2*Nm),
     .		   Gfscl(2,2*Nm),hybl(2*Nm)

	external locate
	include 'Glob_cons'
	do i=-Npi0,Npi0
	  wl_p(i+Npi0+1)=w_pi(i)
	  pi0l(i+Npi0+1)=Pi0_G(i)
	end do
	N2=2*Npi0+1

	do i=-N,N
	  wl(i+N+1)=w(i)
	  gcl(i+N+1)=gc(i)
	  gfl(i+N+1)=gf(i)
	  sigtl(i+N+1)=ssig(i)
	  hybl(i+N+1)=hyb(i)
          do sp=1,2
	   sigmal(sp,i+N+1)=sigma(sp,i)
           gfscl_lp(sp,N+i+1)=Gfsc_lp(sp,i)
           Gfscl(sp,N+i+1)=Gfscript(sp,i)
          end do
	end do
	N1=N

	call makegrid(3)

	do i=-N,N
	  wi=w(i)
	  j=min(2*N1+1,locate(2*N1+1,wl,1,2*N1+1,wi))
	  j1=min(N2,locate(N2,wl_p,1,N2,wi))
	  pi0(i)=(wi-wl_p(j1))*(pi0l(j1+1)-pi0l(j1))/
     .		 (wl_p(j1+1)-wl_p(j1))+pi0l(j1)
	  gc(i)=(wi-wl(j))*(gcl(j+1)-gcl(j))/
     .		 (wl(j+1)-wl(j))+gcl(j)
	  gf(i)=(wi-wl(j))*(gfl(j+1)-gfl(j))/
     .		 (wl(j+1)-wl(j))+gfl(j)
	  ssig(i)=(wi-wl(j))*(sigtl(j+1)-sigtl(j))/
     .		 (wl(j+1)-wl(j))+sigtl(j)
	  hyb(i)=(wi-wl(j))*(hybl(j+1)-hybl(j))/
     .		 (wl(j+1)-wl(j))+hybl(j)
          do sp=1,2
            Gfsc_lp(sp,i)=gfscl_lp(sp,j)+(wi-wl(j))*
     .      (gfscl_lp(sp,j+1)-gfscl_lp(sp,j))/(wl(j+1)-
     .      wl(j))
            Gfscript(sp,i)=Gfscl(sp,j)+(wi-wl(j))*
     .      (Gfscl(sp,j+1)-Gfscl(sp,j))/(wl(j+1)-
     .      wl(j))
            sigma(sp,i)=sigmal(sp,j)+(wi-wl(j))*
     .      (sigmal(sp,j+1)-sigmal(sp,j))/(wl(j+1)-
     .      wl(j))
          end do
        end do



	return
	end

c	**************************************************************
c	************************************************************
	
	subroutine interpol(flag)
	include 'Glob_decl'
	real*8 wl(2*Nm),wi,wl_p(2*Nm)
	complex*16 pi0l(2*Nm),gcl(2*Nm),gfl(2*Nm)
	integer sp,locate,N1,spj,flag,N2,j1
        complex*16 sigmal(2,2*Nm),gfscl_lp(2,2*Nm),sigtl(2*Nm),
     .		   Gfscl(2,2*Nm),hybl(2*Nm)

	external locate
	include 'Glob_cons'
	do i=-Npi0,Npi0
	  wl_p(i+Npi0+1)=w_pi(i)
	  pi0l(i+Npi0+1)=Pi0_G(i)
	end do
	N2=2*Npi0+1

	do i=-N,N
	  wl(i+N+1)=w(i)
	  gcl(i+N+1)=gc(i)
	  gfl(i+N+1)=gf(i)
	  sigtl(i+N+1)=ssig(i)
	  hybl(i+N+1)=hyb(i)
          do sp=1,2
	   sigmal(sp,i+N+1)=sigma(sp,i)
           gfscl_lp(sp,N+i+1)=Gfsc_lp(sp,i)
           Gfscl(sp,N+i+1)=Gfscript(sp,i)
          end do
	end do
	N1=N

	call makegrid(flag)

	do i=-N,N
	  wi=w(i)
	  j=min(2*N1+1,locate(2*N1+1,wl,1,2*N1+1,wi))
	  j1=min(N2,locate(N2,wl_p,1,N2,wi))
	  pi0(i)=(wi-wl_p(j1))*(pi0l(j1+1)-pi0l(j1))/
     .		 (wl_p(j1+1)-wl_p(j1))+pi0l(j1)
	  gc(i)=(wi-wl(j))*(gcl(j+1)-gcl(j))/
     .		 (wl(j+1)-wl(j))+gcl(j)
	  gf(i)=(wi-wl(j))*(gfl(j+1)-gfl(j))/
     .		 (wl(j+1)-wl(j))+gfl(j)
	  ssig(i)=(wi-wl(j))*(sigtl(j+1)-sigtl(j))/
     .		 (wl(j+1)-wl(j))+sigtl(j)
	  hyb(i)=(wi-wl(j))*(hybl(j+1)-hybl(j))/
     .		 (wl(j+1)-wl(j))+hybl(j)
          do sp=1,2
            Gfsc_lp(sp,i)=gfscl_lp(sp,j)+(wi-wl(j))*
     .      (gfscl_lp(sp,j+1)-gfscl_lp(sp,j))/(wl(j+1)-
     .      wl(j))
            Gfscript(sp,i)=Gfscl(sp,j)+(wi-wl(j))*
     .      (Gfscl(sp,j+1)-Gfscl(sp,j))/(wl(j+1)-
     .      wl(j))
            sigma(sp,i)=sigmal(sp,j)+(wi-wl(j))*
     .      (sigmal(sp,j+1)-sigmal(sp,j))/(wl(j+1)-
     .      wl(j))
          end do
        end do



	return
	end

c	**************************************************************

	subroutine kktransf(rhosigma,rlsigma,flag)
	
	integer llim,ulim,flag
	include 'Glob_decl'
	real*8 rhosigma(-Nm:Nm),rlsg(-Nm:Nm),
     .	rlsigma(-Nm:Nm),wo(-Nm:Nm),woj,dspj,resgj,spi,dspi,
     .	dwi,woji1,woji2

c		      oo		
c		       /			
c	rlsigma(w)= -P | dw' rhosigma(w')
c		       /     -----------	
c	             -oo       w' - w


c		KRAMERS-KRONIG

	if(flag.eq.0.or.flag.eq.2) then
	  llim=-1
	  ulim=1		! calculate only ReSig(0)
	else
	  llim=-N+1
	  ulim=N		! calculate ReSig(w) for all w
	end if

	do j=llim,ulim

c	interlacing grid

	wo(j)=.5D0*(w(j-1)+w(j))
	woj=wo(j)
	dspj=rhosigma(j)-rhosigma(j-1)
	resgj=0.0D0

 	
	  do i=-N,N-1

	  if(i.eq.(j-1)) goto 10	! skip the interval (j-1) to j
	  spi=rhosigma(i)
	  dspi=rhosigma(i+1)-rhosigma(i)
	  dwi=w(i+1)-w(i)
	  woji1=w(i)-woj
	  woji2=w(i+1)-woj
	  r1=dlog(woji2/woji1)


	  resgj=resgj-(spi*r1 + dspi )
	  resgj=resgj-(dspi/dwi)*(woj -w(i))*r1
 10	  continue

	  end do

	  resgj=resgj - dspj
	  rlsg(j)=resgj
	 end do
	 rlsg(ulim+1)=rlsg(ulim)
	 rlsg(llim-1)=rlsg(llim)
	
	 do i=llim-1,ulim
	  rlsigma(i)=0.5d0*(rlsg(i)+rlsg(i+1))
	 end do

	 return
	 end


c	**************************************************************

	real*8 function kktsing(rhosigma,wj)
	
	include 'Glob_decl'
	integer sp
	real*8 rhosigma(-Nm:Nm),dwi,wji1,wji2,wj,resgj,spi,dspi

c		 0 
c		 /			
c     kktsing=   | dw' rhosigma(w')	rhosigma=Im(Pi)
c		 /     -----------	
c	        -oo    w' - wj




	resgj=0.0D0

 	
	  do  99 i=-N,-1

	  spi=rhosigma(i)
	  dspi=rhosigma(i+1)-rhosigma(i)
	  dwi=w(i+1)-w(i)
	  wji1=w(i)-wj
	  if(wji1.eq.0.d0) wji1=1.D-10
	  wji2=w(i+1)-wj
	  if(wji2.eq.0.d0) goto 99
	  r1=dlog(dabs(wji2/wji1))


	  resgj=resgj+(spi-(dspi/dwi)*wji1)*r1+dspi

 99	  continue

	
	kktsing=resgj


	return
	end

c	**************************************************************
	
        subroutine takeoutpole(flag)
 
        include 'Glob_decl'
        integer i1,i2,i3,locate,flag,flag_lp,flag_sp
        real*8 lomp,romp,momp,y(-Nm:Nm),xx(-Nm:Nm),m,wl(2*Nm)
        real*8 wji,norm,poleq,sp,sign
        external locate,poleq
        real*8 xxpl(-Nm:Nm),xxline(-Nm:Nm)
        include 'Glob_cons'
 
c       POLE OUTSIDE THE BAND

        do spin=1,2
 
        if(spin.eq.1) then
c          First down spin
           write(6,*) 'DOWN SPIN'
           write(6,*) '********************'
           sp=-1.d0
           llim=1
           ulim=N
        else
c          Now up spin
           write(6,*) 'UP SPIN'
           write(6,*) '********************'
           sp=1.d0
           llim=-N
           ulim=-1
        end if
                                                                                
        norm=0.d0
        do i=-N,N
          wl(i+N+1)=w(i)
          xx(i)=dreal(Gfscript(spin,i))
          y(i)=-dimag(Gfscript(spin,i))/pi
          norm=norm+y(i)*dw(i)
        end do
                                                                                
        r1=0.d0
        i1=llim
        do i=llim+1,ulim
          if(y(i).gt.r1.and.dabs(w(i)).gt.b1) then
            i1=i
            r1=y(i1)
          end if
        end do
        write(6,*) 'spin=',sp,' norm=',norm,'  Y-max',y(i1)
                                                                                
        if(dabs(norm-1.d0).lt.2.D-2.and.y(i1).le.100.d0) then
            Ap(spin)=0.d0
            flag_lp=0
        else
            flag_lp=1
        end if
                                                                                
        if(flag_lp.eq.1) then
                                                                                
        omp(spin)=w(i1)
        momp=w(i1)
        Ap(spin)=1.d0/(1.d0+V**2*(1.d0+t**2/(4.d0*momp**2))/
     .          (momp-t**2/(4.d0*momp))**2)
                                                                                
        write(6,*) 'From Numerics spin=',sp,' Ap(1), Pol. val'
        write(6,*) Ap(spin),y(i1)
                                                                                
        i2=i1
        do while(w(i2).gt.(w(i1)-0.1d0))
          i2=i2-1
        end do
        i3=i1
        do while(w(i3).lt.(w(i1)+0.1d0))
          i3=i3+1
        end do
                                                                                
        !i2=i1-5
        !do while(y(i2).gt.y(i2-1))
         ! i2=i2-1
        !end do
        !i3=i1+5
        !do while(y(i3).gt.y(i3+1))
         ! i3=i3+1
        !end do
                                                                                
        m=(y(i3)-y(i2))/(w(i3)-w(i2))
        do i=i2+1,i3
           r2=y(i2)+m*(w(i)-w(i2))
           y(i)=r2
        end do
        end if          !  flag_lp=1
                                                                                
c       Small pole
 304    i=locate(2*N+1,wl,1,2*N+1,w_m)-N-1
        r1=y(-i)
        i1=-i
        do j=-i+1,i
          if(y(j).gt.r1) then
             i1=j
             r1=y(j)
          end if
        end do
                                                                                
        if(r1.gt.5.d0) then     ! Cut out the small pole
          flag_sp=1
        else
          flag_sp=0
        end if
                                                                                
        if(flag_sp.eq.1) then
          i2=i1-1
          do while(y(i2).gt.y(i2-1))
            i2=i2-1
          end do
          i3=i1+1
          do while(y(i3).gt.y(i3+1))
            i3=i3+1
          end do
                                                                                
          r1=0.d0
          do i=i2,i3
            r1=r1+y(i)*dw(i)
          end do
          write(6,*) 'Small pole weight=',r1
                                                                                
          m=(y(i3)-y(i2))/(w(i3)-w(i2))
          do i=i2+1,i3
             r2=y(i2)+m*(w(i)-w(i2))
             y(i)=r2
          end do
                                                                                
        end if
                                                                                
                                                                                
        if(flag_lp.eq.1) then
        r1=0.d0
        do i=-N,N
          r1=r1+y(i)*dw(i)
        end do
        Ap(spin)=1.d0-r1
                                                                                
        write(6,*) 'Poles Numerics spin=',sp,' omp,Ap'
        write(6,*) omp(spin),Ap(spin)
                                                                                
        end if          ! flag_lp=1
                                                                                
                                                                                
        if(flag.ne.0.and.(flag_lp.eq.1.or.flag_sp.eq.1))
     .  call kktransf(y,xx,1)
                                                                                
        do i=-N,N
          Gfsc_lp(spin,i)=dcmplx(xx(i),-pi*y(i))
        end do
                                                                                
        end do          !  spin=1,2
        return
        end


c	************************************************************
	
	real*8 function poleq(sp,xx)
	real*8 sp,xx
	include 'Glob_decl'


	poleq=xx**3+xx**2*(sp*x-ei-ec)-xx*(t**2/4.d0+ec*
     .	(sp*x-ei)+V**2)-t**2*(sp*x-ei)/4.d0

     	return
	end

c	************************************************************
	subroutine quadratic(x,y,z)
	real*8 x(3)
	complex*16 y(3),z(3),a,b,c

	integer i

	c=((x(2)-x(3))*(y(1)-y(2))-(x(1)-x(2))*(y(2)-y(3)))/
     .	  ((x(1)-x(2))*(x(2)-x(3))*(x(1)-x(3)))
        b=(y(2)-y(3))/(x(2)-x(3))-c*(x(2)+x(3))
	a=y(1)-b*x(1)-c*x(1)**2

c	write(6,*) 'a=',a,' b=',b,' c=',c
c	write(6,*) y(1),a+b*x(1)+c*x(1)**2
c	write(6,*) y(2),a+b*x(2)+c*x(2)**2
c	write(6,*) y(3),a+b*x(3)+c*x(3)**2

c	x0=(-b+dsqrt(b**2-4.0*a*c))/(2.0*c)
c	write(6,*) x0
c	x0=(-b-dsqrt(b**2-4.0*a*c))/(2.0*c)
c	write(6,*) x0

	z(1)=a
	z(2)=b
	z(3)=c

	return
	end
c	************************************************************


	subroutine readinput(flag)
	include 'Glob_decl'
	integer N1,locate,sp,init,flag
	real*8 fermic,solveq,n0up,n0dn
	real*8 wl(2*Nm),dwl(2*Nm),wi,sgn
	complex*16 pi0l(2*Nm),sigl(2,2*Nm)
	character*100 str

	external fermic,solveq,locate,sgn


	include 'Glob_cons'


	  if(flag.eq.1.or.temp.gt.0.d0) then
	  write(6,*) 'Enter U'
	  read(5,*) U
	  write(6,*) 'U=',U
	  end if


	  call makegrid(1)

	  r1=solveq(1,1)

	  do i=-N,N
	    Gc(i)=Gc0(i)
	    Gf(i)=0.5d0*(Gfscript(1,i)+Gfscript(2,i))
	  end do

          !call Grid_Gfscript

          !call makegrid(1)

          r1=solveq(1,1)

	  n0up=0.d0
	  n0dn=0.d0
	  do i=-N,0
	    n0up=n0up+(-dimag(Gfscript(2,i)))*dw(i)*ferm(i)/pi
	    n0dn=n0dn+(-dimag(Gfscript(1,i)))*dw(i)*ferm(i)/pi
          end do
	  n0=n0up+n0dn
	  mu0=n0up-n0dn

	  do i=-N,N
c            write(13,*) w(i),-dimag(Gfscript(1,i))/pi,
c     .      -dimag(Gfscript(2,i))/pi
	    Gc(i)=Gc0(i)
	    Gf(i)=0.5d0*(Gfscript(1,i)+Gfscript(2,i))
	    do sp=1,2
	       Gfsc_lp(sp,i)=Gfscript(sp,i)
	    end do
	  end do
c	  close(13)

          call evalpi0(flag)
	  if(flag.eq.0) then
	  open(unit=20,file='Pi0.dat',status='unknown')
	  do i=-N,N
	    write(20,*) w(i),dimag(pi0(i)),dreal(pi0(i))
	  end do
	  close(20)
	  end if

	  i=1
	  open(unit=23,file='Sigdn.dat',status='unknown')
 9	  read(23,*,end=10) r1,r2,r3
 	  wl(i)=r1
	  sigl(1,i)=dcmplx(r3,sgn(r1)*r2)
	  sigl(1,i)=dcmplx(r3,r2)
	  i=i+1
	  goto 9
 10	  close(23)
 	  N1=i-1

	  i=1
	  open(unit=24,file='Sigup.dat',status='unknown')
 11	  read(24,*,end=12) r1,r2,r3
 	  wl(i)=r1
	  sigl(2,i)=dcmplx(r3,sgn(r1)*r2)
	  sigl(2,i)=dcmplx(r3,r2)
	  i=i+1
	  goto 11
 12	  close(24)
 	  N1=i-1

	  do i=-N,N
	    wi=w(i)
	    j=locate(N1,wl,1,N1,wi)
	    do sp=1,2
	      sigma(sp,i)=sigl(sp,j)+(sigl(sp,j+1)-sigl(sp,j))*
     .		        ((wi-wl(j))/(wl(j+1)-wl(j)))
     	    end do
     	  end do

	  write(6,*) 'Reading done'
	return
	end


c	************************************************************

	subroutine writeoutput
	include 'Glob_decl'
	real*8 sgn

	external sgn
	include 'Glob_cons'

	open(unit=10,file='Gc.dat',status='unknown')
	do i=-N,N
	  write(10,*) w(i),-dimag(gc(i))/pi,dreal(gc(i))
	end do
	close(10)
	open(unit=11,file='Gf.dat',status='unknown')
	do i=-N,N
	  write(11,*) w(i),-dimag(gf(i))/pi,dreal(gf(i))
	end do
	close(11)
	open(unit=12,file='SSig.dat',status='unknown')
	do i=-N,N
	  write(12,*) w(i),dimag(ssig(i)),dreal(ssig(i))
	end do
	close(12)
	open(unit=13,file='Pi0.dat',status='unknown')
	do i=-Npi0,Npi0
	  write(13,*) w_pi(i),dimag(pi0_G(i)),
     .	  dreal(pi0_G(i))
	end do
	close(13)
	open(unit=14,file='Sigdn.dat',status='unknown')
	do i=-N,N
	  write(14,*) w(i),dimag(sigma(1,i)),
     .	  dreal(sigma(1,i))
	end do
	close(14)
	open(unit=15,file='Sigup.dat',status='unknown')
	do i=-N,N
	  write(15,*) w(i),dimag(sigma(2,i)),
     .	  dreal(sigma(2,i))
	end do
	close(15)
	open(unit=16,file='Gfscup.dat',status='unknown')
	do i=-N,N
	  write(16,*) w(i),-dimag(Gfscript(2,i))/pi,
     .	  dreal(Gfscript(2,i))
	end do
	close(16)
	open(unit=17,file='Gfscdn.dat',status='unknown')
	do i=-N,N
	  write(17,*) w(i),-dimag(Gfscript(1,i))/pi,
     .	  dreal(Gfscript(1,i))
	end do
	close(17)
        open(unit=18,file='inputen.dat',status='unknown')
        do i=-N,N
        write(18,*)float(i)*8.d0/dfloat(N)
        enddo
        close(18) 
	return
	end

c	************************************************************
	subroutine writeGfsclp
	include 'Glob_decl'
	include 'Glob_cons'

	
	  open(unit=38,file='Rhosc_lp.dat',status='unknown')
	  open(unit=39,file='Re_lp.dat',status='unknown')
          do i=-N,N
            write(38,*) w(i),-dimag(Gfsc_lp(1,i))/pi,
     .      -dimag(Gfsc_lp(2,i))/pi
            write(39,*) w(i),dreal(Gfsc_lp(1,i)),
     .      dreal(Gfsc_lp(2,i))
          end do
	  close(38)
	  close(39)
	  open(unit=42,file='Polwt.dat',status='unknown')
          write(42,*) Ap(1),Ap(2),omp(1),omp(2)
          write(42,*) 'Ap(1),Ap(2),omp(1),omp(2)'
          close(42)

	return
	end

c	************************************************************
	subroutine writeGfscript
	include 'Glob_decl'
	integer u1,u2
	include 'Glob_cons'

	open(unit=16,file='Gfscup.dat',status='unknown')
	do i=-N,N
	  write(16,*) w(i),-dimag(Gfscript(2,i))/pi,
     .	  dreal(Gfscript(2,i))
	end do
	close(16)
	open(unit=17,file='Gfscdn.dat',status='unknown')
	do i=-N,N
	  write(17,*) w(i),-dimag(Gfscript(1,i))/pi,
     .	  dreal(Gfscript(1,i))
	end do
	close(17)

	return
	end

c	************************************************************
	subroutine writeSigGf
	include 'Glob_decl'
	include 'Glob_cons'

	open(unit=23,file='Sigdn.dat',status='unknown')
	open(unit=24,file='Sigup.dat',status='unknown')
	do i=-N,N
  	  write(23,*) w(i),dimag(sigma(1,i)),dreal(sigma(1,i))
  	  write(24,*) w(i),dimag(sigma(2,i)),dreal(sigma(2,i))
	end do
	close(23)
	close(24)

	open(unit=30,file='Gf.dat',status='unknown')
	open(unit=31,file='Gc.dat',status='unknown')
	do i=-N,N
	  write(30,*) w(i),-dimag(gf(i))/pi,dreal(gf(i))
	  write(31,*) w(i),-dimag(gc(i))/pi,dreal(gc(i))
	end do
	close(30)
	close(31)

	return
	end

c	************************************************************
	subroutine initialise(flag)
	include 'Glob_decl'
	real*8 solveq,n0up,n0dn
	integer flag,sp

	external solveq
	include 'Glob_cons'
	  
	  r1=solveq(1,1)

	  do i=-N,N
	     Gf(i)=0.5d0*(Gfscript(1,i)+Gfscript(2,i))
	  end do

          call Grid_Gf(0)

          call makegrid(1)

          r1=solveq(1,1)

	  do i=-N,N
	     Gf(i)=0.5d0*(Gfscript(1,i)+Gfscript(2,i))
	  end do

	  n0up=0.d0
	  n0dn=0.d0
	  do i=-N,0
	     n0up=n0up+(-dimag(Gfscript(2,i)))*dw(i)*ferm(i)/pi
	     n0dn=n0dn+(-dimag(Gfscript(1,i)))*dw(i)*ferm(i)/pi
          end do
	  n0=n0up+n0dn
	  mu0=n0up-n0dn

	  do i=-N,N
            Gc(i)=Gc0(i)
	    Gf(i)=0.5d0*(Gfscript(1,i)+Gfscript(2,i))
	    pi0(i)=zero
	    do sp=1,2
	       Gfsc_lp(sp,i)=Gfscript(sp,i)
	       sigma(sp,i)=zero
	    end do
	  end do
	  Ap(1)=0.d0
	  Ap(2)=0.d0

          call evalPi0(flag)
	  open(unit=20,file='Pi0.dat',status='unknown')
	  do i=-N,N
	    write(20,*) w(i),dimag(pi0(i)),dreal(pi0(i))
	  end do
	  close(20)

	  return
	  end

c	************************************************************

        subroutine imagpi0_ft(rho1,rho2,impi0)

        include 'Glob_decl'
        integer ljvec(-Nm:Nm),llim,ulim,lj,flag,locate
        integer jmin,jmax,ierr
        real*8 wji,sgn,fji,wjmin,wjmax
        real*8 rho1(-Nm:Nm),rho2(-Nm:Nm),d2(-Nm:Nm),rho2ji,
     .         impi0(-Nm:Nm),wl(2*Nm),fermic
        real time1,time2,time3

c       Finds the imaginary part of Pi0  (Retarded)

        external fermic,sgn,locate
        include 'Glob_cons'

        time1=secnds(0.0)

        do i=-N,N-1
            wl(i+N+1)=w(i)
            d2(i)=(rho2(i+1)-rho2(i))/(w(i+1)-w(i))
        end do
        wl(2*N+1)=w(N)
        d2(N)=d2(N-1)

        do 98 i=-N,N
          llim=-N
          wjmin=max(w(i)-w(N),-w(N))
          wjmax=min(w(i)+w(N),w(N))
          jmin=locate(2*N+1,wl,1,2*N+1,wjmin)-N-1+1
          jmax=locate(2*N+1,wl,1,2*N+1,wjmax)-N-1
          r1=0.d0
          llim=-N
          do j=jmin,jmax
            wji=w(j)-w(i)
            do while(wji.gt.w(llim+1))
              llim=min(llim+1,N)
            end do
            lj=llim
c           This implies that  ==> w(lj)<= w_j-w_i <=w(lj+1)
            rho2ji=rho2(lj)+d2(lj)*(wji-w(lj))
	    fji=fermic(wji,beta)
            r1=r1+rho1(j)*rho2ji*(ferm(j)-fji)*dw(j)
          end do
          impi0(i)=-pi*r1
 98     continue
        close(15)

        time2=secnds(time1)
        write(6,*) 'convol',time2

        return
        end

c	**************************************************************

        subroutine convolsg_ft(rho1,rho2,totsig,flag,spin)

        include 'Glob_decl'
        integer ljvec(-Nm:Nm),llim,ulim,lj,flag,spin,locate,
     .  li,ui,sp,spbar,jmin,jmax,unit
        real*8 wi,wji,sgn,theta,sign,piji,spi,fom,fermic,spj
        real*8 rho1(-Nm:Nm),rho2(-Nm:Nm),d2(-Nm:Nm),rho2ji,
     .         wl(2*Nm),resig(-Nm:Nm),wjmin,wjmax,
     .         rhosigma(-Nm:Nm),fji,fj,fi
        complex*16 totsig(-Nm:Nm)

        external sgn,theta,fermic
        external locate
        include 'Glob_cons'

c       Retarded self energies

        do i=-N,N-1
            wl(i+N+1)=w(i)
            d2(i)=(rho2(i+1)-rho2(i))/(w(i+1)-w(i))
        end do
        wl(N+N+1)=w(N)
        d2(N)=d2(N-1)

        do 198 i=-N,N
          wjmin=max(-w(N),-w(N)-w(i))
          wjmax=min(w(N)-w(i),w(N))
          jmin=locate(2*N+1,wl,1,2*N+1,wjmin)-N-1+1
          jmax=locate(2*N+1,wl,1,2*N+1,wjmax)-N-1
          llim=-N
          r1=0.d0
          do j=jmin,jmax
            wji=w(i)+w(j)
            do while(wji.gt.w(llim+1))
              llim=min(llim+1,N)
            end do
            lj=llim
c           This implies that  ==> w(lj)<= w_j+w_i <=w(lj+1)
            rho2ji=rho2(lj)+(rho2(lj+1)-rho2(lj))*((wji-w(lj))/
     .                                          (w(lj+1)-w(lj)))
	    fji=fermic(wji,beta)
            r1=r1+rho1(j)*rho2ji*dw(j)*dabs(theta(-w(j))-
     .            fji)
          end do
          rhosigma(i)=r1/pi
          r1=rhosigma(i)
          if(r1.lt.0.d0) then
             write(6,*) 'rhosigma negative ',i,w(i),r1
             stop
          end if
 198    continue


c       Large pole contribution
        if(spin.eq.-1) then
          sp=1
          spbar=2
        else
          sp=2
          spbar=1
        end if

        if(Ap(spbar).gt.0.d0) then

        fom=fermic(omp(spbar),beta)
        do i=-N,N
          wi=(omp(spbar)-w(i))
          if(dabs(wi).lt.w(N).and.wi.ne.0.d0) then
          j=locate(2*N+1,wl,1,2*N+1,wi)-N-1
          r1=rho1(j)+(rho1(j+1)-rho1(j))*((wi-w(j))/
     .    (w(j+1)-w(j)))
          rhosigma(i)=rhosigma(i)+Ap(spbar)*r1
     .    *dabs(theta(-wi)-fom)/pi
          end if
          r1=rhosigma(i)
          if(r1.lt.0.d0) then
             write(6,*) 'Pole contrib rhosigma negative ',i,w(i),r1
             stop
          end if
        end do

	end if

 401    write(6,*) 'CONVOLUTION DONE, NOW KKTRANSF'
        call kktransf(rhosigma,resig,flag)

        do i=-N,N
           totsig(i)=resig(i)-pi*ii*rhosigma(i)
        end do

        return
        end

c	**************************************************************

                          
