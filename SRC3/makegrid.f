c	************************************************************


	subroutine makegrid(flag)

	include 'Glob_decl'
	integer grid,M1,M2,M3,M4_l,M5_l,M5_r,
     .	M4_r,M6,M7,lN,rN,i1,sp,i2,i3,j1,j2,M8,M9,flag
	real*8 k,kp,wmax,wmax1,tep,uw,fermic
        real*8 wl(-Nm:Nm),dwl(-Nm:Nm),b_l,db_l,db_r,b_r

	external fermic,grid

	wmax=3.0d0*max(t,x)
	wmax1=15.0d0*max(t,x)
c	if flag=0  Nonverbose mode
c	       =1  Verbose mode

	if(frac.ge.100)then
        N=5000
	ep0=0.001d0
        do i=-N,N
	  w(i)=dfloat(i)*ep0
         dw(i)=ep0
          ferm(i)=fermic(w(i),beta)
       	end do
       	return
	end if
c       SET UP THE FREQUENCY GRID

        w(0)=0.0d0
	dw(0)=0.d0

	ep1=ep0*(b1-ep2)/(ep2-ep0)
	k=dlog((b1-ep0)/(b1-ep2))
	kp=(1.d0-dexp(-k))
	M1=aint(dlog(b1/ep1+1.d0)/k)
	if((M1/2*2-M1).ne.0) M1=M1+1
        do i = 1,M1
           w(i)=ep1*(dexp(i*k)-1.d0)
           w(-i)=-w(i)
           dw(i)=kp*(w(i)+ep1)
           dw(-i)=dw(i)
        enddo
	N=M1
	
	do i=-N,N
	     wl(i)=w(i)
	     dwl(i)=dw(i)
	end do

	i=0
	do while(w(i).lt.wm_l)
	    i=i+1
	end do
	i1=i

	do while(w(i).lt.wm_r)
	    i=i+1
	end do
	i2=i

	tep=min(dw(i1),(w(i2)-w(i1))*ep_wm)
	j=i1
	do while(wl(j).le.wm_r)
	    j=j+1
	    wl(j)=wl(j-1)+tep
	    dwl(j)=tep
	end do
	j1=j-1

	do i=i2+1,N
	    wl(j1+i-i2)=w(i)
	    dwl(j1+i-i2)=dw(i)
	end do
	N=j1+N-i2

 501	do i=0,N
	    w(i)=wl(i)
	    dw(i)=dwl(i)
	    w(-i)=-wl(i)
	    dw(-i)=dwl(i)
	end do

	M1=N

	uw=max(w(N),b2-db2)
	M2=grid(w(N),uw,ep2,1.d0)
	N=N+M2

	M3=grid(w(N),2.d0*db2+w(N),ep3,1.d0)	! ImSigma
	N=N+M3

c	Given b3_p, db3_p, b3_n, db3_n

	if(b3_p.lt.b3_n) then
	  b_l=b3_p
	  db_l=db3_p
	  b_r=b3_n
	  db_r=db3_n
	else
	  b_l=b3_n
	  db_l=db3_n
	  b_r=b3_p
	  db_r=db3_p
	end if
	

	uw=max(w(N),b_l-db_l)
	M4_l=grid(w(N),uw,ep2,1.d0)
	N=N+M4_l

	M5_l=grid(w(N),b_l+db_l,ep3,1.d0)	! Gf, Gfscript
	N=N+M5_l

	uw=max(w(N),b_r-db_r)
	M4_r=grid(w(N),uw,ep2,1.d0)
	N=N+M4_r

	M5_r=grid(w(N),b_r+db_r,ep3,1.d0)	! Gf, Gfscript
	N=N+M5_r

	uw=max(w(N),b4-db4)
	M6=grid(w(N),uw,ep2,1.d0)
	N=N+M6

	M7=grid(w(N),2.d0*db4+w(N),ep4,1.d0)	! ImPi0
	N=N+M7

	uw=max(w(N),wmax)
	M8=grid(w(N),wmax,ep2,1.d0)
	N=N+M8

	M9=grid(w(N),wmax1,ep5,1.d0)
	N=N+M9

	do i=-N,0
	  w(i)=-w(-i)
	  dw(i)=dw(-i)
	end do

	if(flag.eq.1) then
	  write(6,*) M1,M2,M3,M4_l,M5_l,M4_r,M5_r
	  write(6,*) M6,M7,M8,M9
	  write(6,*) 'N=',N
	end if
	do i=-N,N
	   write(12,*) w(i),dw(i)
	end do
	close(12)

	if(N.gt.Nm) then
	  write(6,*) 'N > Nm: segmentation fault'
	  stop
	end if

c	Defining the fermi function as an array.
 400	do i=-N,N
	     ferm(i)=fermic(w(i),beta)
	end do

 401	return
	end

c	************************************************************
	
	integer function grid(aa,bb,ep,sign)
	include 'Glob_decl'

	integer M1
	real*8 aa,bb,ep,sign
	
	j=dint(sign)
	M1=aint(sign*(bb-aa)/ep)
	if((M1/2*2-M1).ne.0) M1=M1+1
	do i=N+1,N+M1
	    w(j*i)=aa+sign*ep*dfloat(i-N)
	    dw(j*i)=ep
	end do
	
	grid=M1

	return
	end

