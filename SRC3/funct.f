c       **********************************************
        double precision function theta(x)
        double precision x

        if(x.lt.0.D0) then
          theta=0.D0
        else if(x.eq.0.d0) then
          theta=0.5d0
        else
          theta=1.d0
        end if

        return
        end

c       **********************************************

	double precision function sgn(x)
	double precision x

	if(x.lt.0.D0) then
	  sgn=-1.D0
	else
	  sgn=1.D0
	end if

	return
	end

c	**********************************************	  

	double precision function fermic(x,beta)
	double precision x,beta
c	fermi function which avoids underflow/overflow.  If beta=0,
c	it is assumed that T=0 is meant!

	if(beta.eq.0.0D0) then
	  if(x.lt.0.0D0) then	  
	    fermic=1.0D0
	  else if(x.eq.0.0D0) then
	    fermic=0.5D0
	  else
	    fermic=0.0D0
	  end if
	else
	  if(x.lt.0.0D0) then	  
	    fermic=1.0D0/(dexp(beta*x)+1.0D0)
	  else
	    fermic=dexp(-beta*x)/(dexp(-beta*x)+1.0D0)
	  end if
	end if
	return
	end

c	**********************************************	  
	double precision function dfermic(x,beta)
	double precision x,beta
c
c	This function returns the minus derivative of the fermi function
c
c                               1
c	fermic(x,beta)= ----------------
c                       exp(-beta x) + 1
c
c       d fermic(x,beta)     -beta exp(-beta x)
c	----------------  = --------------------- = -dfermic(x,beta)
c              d x           (exp(-beta x) + 1)^2
c
c
	if(x.lt.0.0D0) then	  
	  dfermic=beta*dexp(beta*x)/(dexp(beta*x)+1.0D0)**2
        else
	  dfermic=beta*dexp(-beta*x)/(dexp(-beta*x)+1.0D0)**2
	end if

	return
	end

c	**********************************************	  
c
	INTEGER FUNCTION LOCATE(N,xx,ll,ul,x)
	integer N,ll,ul
	real*8 xx(1:N),x
c
c	given an array xx(1:n),and given a value x, returns a value j 
c 	such that x is between xx(j) and xx(j+1). xx(1:n) must be
c 	monotonic.

	integer jl,jm,ju


	if(x.lt.xx(ll)) then
	  locate=ll
c	  pause'input out of range in locate,left'
	  return
	end if

	if(x.gt.xx(ul)) then
	  locate=ul-1
c	  pause'input out of range in locate,right'
	  return
	end if

	jl=ll
	ju=ul
        do while((ju-jl) .gt. 1) 
	  jm=(ju+jl)/2
	  if (x .ge. xx(jm)) then
	    jl=jm
	  else
	    ju=jm
	 endif
	end do 
	locate=jl
c
	return
	end

c	****************************************************************
	double complex function gauss(z)
c	This block calculates -i*sqrt(pi)*w(z)
	double complex z,ii
	double precision pi2,x,y
	logical flag
	parameter(pi2=1.7724539D0)
	ii=dcmplx(0.D0,1.D0)
	if(dimag(z).lt.0.0D0) then
	  call wofz(dreal(z),-dimag(z),x,y,flag)
	  gauss=-pi2*ii*dcmplx(x,-y)
	else
	  call wofz(dreal(z),dimag(z),x,y,flag)
	  gauss=-pi2*ii*dcmplx(x,y)
	end if

	if(flag) write(6,*) 'error in cerfjd'
	return
	end


c	**************************************************
c		WS(Z)
c	**************************************************

	  double complex function ws(z)
	  double complex z
	  integer nnorm,idos,iflag1,iflag2,levmax
	  double precision dband,epsvh,zr,zi,acc,dnorm
	  double precision tdnorm,error1,area1,rho1,rho2,rho0,
     &    y1min,yplus,yminus,y1max,y2min,y2max,wsr,wsi,error2,
     &    area2
          double complex sem,gauss,flat

	  common/dos/dband,epsvh,zr,zi,acc,dnorm,nnorm,idos
	  parameter(pi=3.1415927D0)

	  external smsn
	  external rho1,rho2,rho0
	  external sem
	  external gauss
	  external flat

	   
	   zr=dreal(z)
	   zi=dimag(z)
	   z=dcmplx(zr,zi)

	  if(idos.eq.2) then
 	  	ws=sem(z)
	        goto 354
	  else if(idos.eq.1) then
	        ws=gauss(z)
		goto 354
	  else if(idos.eq.0) then
		ws=flat(z)
		goto 354
	  end if	
	  if (nnorm .eq. 1) then
	   tdnorm=0.0D0
	   call smsn(rho0,-dband,dband,acc,tdnorm,error1,area1,iflag1,levmax)
	   dnorm=tdnorm
	   nnorm=0
	  end if




	      y1min=0.0D0
	      yplus=(dband-zr)
	      yminus=(dband+zr)
	      y1max=dmax1(yplus,yminus)

	   if (zr .gt. dband) y1min =zr-dband
	   if (zr .lt. -dband) y1min = -dband-zr

	  wsr=0.0D0
          call smsn(rho1,y1min,y1max,acc,wsr,error1,area1,iflag1,levmax)
	  
	  
	  if (zi .eq. 0.0D0) then 
	   wsi=-pi*rho0(zr)
	  else
	   y2min=datan((-dband-zr)/zi)
	   y2max=datan((dband-zr)/zi)
	   wsi=0.0D0

	   call smsn(rho2,y2min,y2max,acc,wsi,error2,area2,iflag2,levmax)
	  end if


	  ws=dcmplx(wsr,wsi)


 354	  return
	  end


c	  ***********************************************************
        complex*16 function sem(z)
c       This block calculates Hilbert transform for a semicircular 
c       DOS of bandwidth=2*dband
c       In general sem=8*(z-sqrt(z**2-D**2/4))/D**2
        double complex z,z1,z2,ii,sem1,sem2,sz2
        double precision dband,pi

	dband=1.0d0
        include 'Glob_cons'

        z1=dcmplx(dband,0.D0)
					             
        z2=z**2-z1**2
				                  
        sz2=zsqrt(z2)
												               
        sem1=2.D0/(z+sz2)
												                    
        sem2=2.D0/(z-sz2)
        
   
        if(dimag(sem1).le.0.D0) then
             sem=sem1
        else if(dimag(sem2).le.0.D0) then
             sem=sem2
        else 
             write(6,*) 'no causal root found'
             stop 
        end if
    
        return
        end  


c	*************************************************************
	
	  double complex function flat(z)
	  double complex z,ii
	  double precision rews,imws,rez,imz,pi
	  double precision dband,epsvh,zr,zi,acc,dnorm,rho0

	  common/dos/dband,epsvh,zr,zi,acc,dnorm,nnorm,idos

	
c	  Flat Band from -D/2 to D/2 where 1.D0
	  
	  include 'Glob_cons'

	  rez=dreal(z)
	  imz=dimag(z)
	  
	  if(dabs(rez).eq.dband/2.D0) rez=rez+1.D-6
	  if(dabs(imz).le.1.D-6) then
	    rews=dlog(dabs((rez+dband/2.D0)/(rez-dband/2.D0)))
	    if(dabs(rez).ge.dband/2.D0) then
	      imws=0.D0
	    else
	    imws=-pi
	    end if
	  else
	    rews=0.5D0*dlog(((rez+dband/2.D0)**2+imz**2)/
     &	         	   ((rez-dband/2.D0)**2+imz**2))

	    imws=-(datan((dband/2.D0-rez)/imz)+
     &		  datan((dband/2.D0+rez)/imz))
	  end if

 355	  flat=dcmplx(rews,imws)

 	  return
	  end

c	*************************************************************
	  double precision function rho1(y)
	  double precision y
	  integer nnorm,idos
	  double precision epy1,emy1
	  double precision dband,epsvh,zr,zi,acc,dnorm,rho0
	  common/dos/dband,epsvh,zr,zi,acc,dnorm,nnorm,idos
	  
	  external rho0
	  
	  rho1=0.0D0


	  epy1=zr+y
	  emy1=zr-y
	  
	  rho1=-(rho0(epy1)-rho0(emy1))*y/(y*y+zi*zi)
	  return
	  end

c	  ***********************************************************
	  double precision function rho2(y)
	  double precision y

	  integer nnorm,idos
	  double precision dband,epsvh,zr,zi,acc,dnorm,rho0
	  double precision ey2
	  common/dos/dband,epsvh,zr,zi,acc,dnorm,nnorm,idos
	  
	  external rho0
	  
	  ey2=zr+dtan(y)*zi

	  rho2=-rho0(ey2)
	  return
	  end

c
c	  ***********************************************************
c	  ***********************************************************
	  double precision function rho0(eps)
	  double precision eps
	  integer nnorm,idos
	  double precision dband,epsvh,zr,zi,acc,dnorm
	  double precision a0,a1,a2,a3,a4
	  common/dos/dband,epsvh,zr,zi,acc,dnorm,nnorm,idos

	  rho0=0.0D0
	  if (dabs(eps) .gt. dband) return
	  if(idos.eq.0) then
	    rho0=1.0D0
	  else if(idos.eq.1) then
            rho0=dexp(-eps**2)
          else if(idos.eq.2) then
           rho0=dsqrt(1.0D0-(eps/dband)**2)
          else if(idos.eq.3) then
	   rho0=dlog(2.0D0*dband/dsqrt((eps-epsvh)**2+1.d-8))
	   rho0=rho0*dsqrt(1.0D0-(eps/dband)**2)
	  else if(idos.eq.4) then
	   a0=-0.5236D0
	   a1=-0.6014D0
	   a2=+0.2516D0
	   a3=-0.0953D0
	   a4=+0.1964D0
	   rho0=a0*log10(dsqrt((eps-a1)**2+0.000001D0))+a2+a3*eps+a4*eps**2
	  else if(idos.eq.5) then
	   a0=-0.461112D0
	   a1=0.290368D0
	   a2=0.028830D0
	   a3=0.000394D0
	   rho0 = a0*log10(dsqrt(eps**2+a3*a3))+a1+a2*eps**2
	  else if(idos.eq.6) then
	   a0=-0.478585D0
	   a1=-0.201103D0
	   a2=+0.278186D0
	   a3=-0.027951D0
	   a4=+0.056084D0
	   rho0=a0*log10(dsqrt((eps-a1)**2+0.000001D0))+a2+a3*eps+a4*eps**2
	  end if
	  rho0=rho0/dnorm
	  return
	  end

c	****************************************************************
c
c 
c
      SUBROUTINE SMSN(F,A,B,ACC,ANS,ERROR,AREA,IFLAG,LEVMAX)
c      IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT INTEGER*4 (I-N)	
C     THIS ROUTINE IS COPIED FROM NUMERICAL COMPUTING BY SHAMPINE
C     AND ALLEN, W.B.SAUNDERS&CO.,PHILA.,PA. 1973
C
C          IS AN ADAPTIVE, ITERATIVE CODE BASED ON SIMPSON'S RULE
C     IT IS DESIGNED TO EVALUATE THE INDEFINITE INTEGRAL OF A CON-
C     TINUOUS FUNCTION WITH FINITE LIMITS OF INTEGRATION.
C
C
C     F - NAME OF FUNCTION WHOSE INTEGRAL IS DESIRED. THE FUNCTION
C         NAME F MUST APPEAR IN AN EXTERNAL STATEMENT IN THE CALLING
C         PROGRAM.
C     A,B - LOWER AND UPPER LIMITS OF INTEGRATION.
C     ANS - APPROXIMATE VALUE OF THE INTEGRAL OF F(X) FROM A TO B
C     AREA - APPROXIMATE VALUE OF THE INTEGRAL OF DABS(F(X)) FROM
C            A TO B
C     ERROR - ESTIMATED ERROR OF ANS. USER MAY WISH TO EXTRA-
C             POLATE BY FORMING ANS+ERROR TO GET WHAT IS OFTEN BUT NOT
C             ALWAYS A MORE ACCURATE RESULT.
C     ACC - DESIRED ACCURACY OF ANS. CODE TRIES TO MAKE
C           DABS(ERROR).LE.ACC*DABS(AREA).
C
C     IFLAG = 1 FOR NORMAL RETURN
C             2 IF IT NECESSARY TO GO TO 500 LEVELS OR USE A
C               SUBINTERVAL TOO SMALL FOR MACHINE WORD LENGTH.
C               ERROR MAY BE UNRELIABLE IN THIS CASE.
C             3 IF MORE THAN 2000 FUNCTION EVALUATIONS ARE USED.
C               ROUGH APPROXIMATIONS ARE USED TO COMPLETE THE
C               COMPUTATIONS AND ERROR IS USUALLY UNRELIABLE.
C
C
      INTEGER IFLAG,LEVMAX
      DOUBLE PRECISION FV(5),LORR(500),F1T(500),F2T(500),F3T(500),
     .DAT(500),ARESTT(500),ESTT(500),EPST(500),PSUM(500),DIFF,AREA,
     .DX,ALPHA,F,A,B,ACC,ANS,ERROR
      EXTERNAL F
C
C     SET U TO APPROXIMATELY THE UNIT ROUND OFF OF SPECIFIC MACHINE
C
      U=1.0d-16
C
C     INITIALIZE
C
      FOURU=4.0D0*U
      IFLAG=1
      EPS=ACC
      ERROR=0.0D0
      LVL=1
      LEVMAX=1
      LORR(LVL)=1
      PSUM(LVL)=0.0D0
      ALPHA=A
      DA=B-A
      AREA=0.0D0
      AREST=0.0D0
      FV(1)=F(ALPHA)
      FV(3)=F(ALPHA+0.5D0*DA)
      FV(5)=F(ALPHA+DA)
      KOUNT=3
      WT=DA/6.0D0
      EST=WT*(FV(1)+4.0D0*FV(3)+FV(5))
C
C     BASIC STEP: HAVE ESTIMATE EST OF INTEGRAL ON (ALPHA,ALPHA+DA).
C     BISECT AND COMPUTE ESTIMATES OF LEFT AND RIGHT HALF INTERVALS.
C     SIMILARLY TREAT INTEGRAL OF DABS(F(X))
C     SUM IS BETTER VALUE FOR INTEGRAL AND DIFF/15.0 IS APPROX. ERROR.
C
1     DX=0.5D0*DA
      FV(2)=F(ALPHA+0.5D0*DX)
      FV(4)=F(ALPHA+1.5D0*DX)
      KOUNT=KOUNT+2
      WT=DX/6.0D0
      ESTL=WT*(FV(1)+4.0D0*FV(2)+FV(3))
      ESTR=WT*(FV(3)+4.0D0*FV(4)+FV(5))
      SUM=ESTL+ESTR
      ARESTL=WT*(DABS(FV(1))+DABS(4.0D0*FV(2))+DABS(FV(3)))
      ARESTR=WT*(DABS(FV(3))+DABS(4.0D0*FV(4))+DABS(FV(5)))
      AREA=AREA+((ARESTL+ARESTR)-AREST)
      DIFF=EST-SUM
C
C     IF ERROR IS READABLE, GO TO 2. IF THE INTERVAL IS TOO SMALL OR
C     TOO MANY LEVELS OR TOO MANY FUNCTION EVALUATIONS, SET
C     A FLAG AND GO TO 2 ANYWAY.
C
      IF(DABS(DIFF).LE.EPS*DABS(AREA))GO TO 2
      IF(DABS(DX).LE.FOURU*DABS(ALPHA))GO TO 5
      IF(LVL.GE.500)GO TO 5
      IF(KOUNT.GE.2000)GO TO 6
C
C     HERE TO RAISE LEVEL. STORE INFORMATION TO PROCESS RIGHT HALF
C     INTERVAL LATER. INITIALIZE FOR 'BASIC STEP' SO AS TO TREAT
C     LEFT HALF INTERVAL.
C
      LVL=LVL+1
      LORR(LVL)=0
      F1T(LVL)=FV(3)
      F2T(LVL)=FV(4)
      F3T(LVL)=FV(5)
      DA=DX
      DAT(LVL)=DX
      AREST=ARESTL
      ARESTT(LVL)=ARESTR
      EST=ESTL
      ESTT(LVL)=ESTR
      EPS=EPS/1.4D0
      EPST(LVL)=EPS
      FV(5)=FV(3)
      FV(3)=FV(2)
      GO TO 1
C
C     READ APPROXIMATE INTEGRAL SUM. IF IT WAS ON THE LEFT INTERVAL,
C     GO TO 'MOVE RIGHT'. IF A RIGHT INTERVAL, ADD RESULT TO FINISH
c
C     AT THIS LEVEL. ARRAY LORR(MNEMONIC FOR LEFT OR RIGHT)
C     TELLS WHETHER LEFT OR RIGHT INTERVAL AT EACH LEVEL.
C
2     ERROR=ERROR+DIFF/15.0D0
      IF(LVL.GT.LEVMAX)LEVMAX=LVL
3     IF(LORR(LVL).EQ.0)GO TO 4
      SUM=PSUM(LVL)+SUM
      LVL=LVL-1
      IF(LVL.GT.1)GO TO 3
      ANS=SUM
      RETURN
C
C     'MOVE RIGHT' RESTORE SAVED INFORMATION TO PROCESS RIGHT HALF
C     INTERVAL.
C
4     PSUM(LVL)=SUM
      LORR(LVL)=1
      ALPHA=ALPHA+DA
      DA=DAT(LVL)
      FV(1)=F1T(LVL)
      FV(3)=F2T(LVL)
      FV(5)=F3T(LVL)
      AREST=ARESTT(LVL)
      EST=ESTT(LVL)
      EPS=EPST(LVL)
      GO TO 1
C
C
C     READ POOR VALUE. SET APPROPRIATE FLAGS.
C
5     IFLAG=2
      GO TO 2
6     IFLAG=3
      GO TO 2
      END
c
c	***************************************************************
C  (C) Copr. 1986-92 Numerical Recipes Software 3ki-i]3|-9.
c
c
C      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 1, PP. 47.
      SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
C
C  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
C  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
C  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
C  MEANS SQRT(-1).
C  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
C  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
C  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
C  OF THE FUNCTION.
C  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
C
C
C  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
C     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
C                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
C                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
C                FLOATING-POINT ARITHMETIC
C     RMAXEXP  = LN(RMAX) - LN(2)
C     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
C                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
C  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
C  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
C
C
C  PARAMETER LIST
C     XI     = REAL      PART OF Z
C     YI     = IMAGINARY PART OF Z
C     U      = REAL      PART OF W(Z)
C     V      = IMAGINARY PART OF W(Z)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  XI, YI      ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
C
C  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
C  PUT TO 0 UPON UNDERFLOW;
C
C  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
C  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
C
*
*
*
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      LOGICAL A, B, FLAG
      PARAMETER (FACTOR   = 1.12837916709551257388D0,
     *           RMAXREAL = 0.5D+154,
     *           RMAXEXP  = 708.503061461606D0,
     *           RMAXGONI = 3.53711887601422D+15)
*
      FLAG = .FALSE.
*
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
*
C
C     THE FOLLOWING IF-STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
*
      QRHO = X**2 + Y**2
*
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
*
      A     = QRHO.LT.0.085264D0
*
      IF (A) THEN
C
C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
C  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
        QRHO  = (1-0.85*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0/J
        YSUM  = 0.0D0
        DO 10 I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0/J
 10     CONTINUE
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
*
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
*
      ELSE
C
C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
C  CONTINUED FRACTION
C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
C  TO OBTAIN THE REQUIRED ACCURACY
C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
C
*
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = IDINT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
*
        B = (H.GT.0.0)
*
        IF (B) QLAMBDA = H2**KAPN
*
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
*
        DO 11 N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
 11     CONTINUE
*
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
*
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
*
      END IF
*
*
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C
*
      IF (YI.LT.0.0) THEN
*
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
*
C
C         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
C         AGAINST OVERFLOW
C
          IF ((YQUAD.GT.RMAXGONI).OR.
     *        (XQUAD.GT.RMAXEXP)) GOTO 100
*
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
*
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
*
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END

c	************************************************************

      SUBROUTINE LINT(XA,YA,Y2A,N,KLO,X,Y) 
      INTEGER N,KLO,KHI,K
      DOUBLE PRECISION XA(N),YA(N),Y2A(N),X,Y,H,B
      
      IF(KLO.EQ.0) THEN		! MEANS UNSET IN MAIN PROGRAM
      KLO=1 
      END IF
      KHI=N 
1     IF (KHI-KLO.GT.1) THEN 
        K=(KHI+KLO)/2 
        IF(XA(K).GT.X)THEN 
          KHI=K 
        ELSE 
          KLO=K 
        ENDIF 
      GOTO 1 
      ENDIF 
      H=XA(KHI)-XA(KLO) 
      IF (H.EQ.0.D0) PAUSE 'Bad XA input.' 
      B=(X-XA(KLO))
      Y=YA(KLO)+Y2A(KLO)*B
      RETURN 
      END 

c	**********************************************	  
      SUBROUTINE LINE(X,Y,N,Y2) 
      PARAMETER (NMAX=5000) 
      INTEGER N,I
      DOUBLE PRECISION X(N),Y(N),Y2(N)

      DO I=1,N-1
        Y2(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
      END DO
      Y2(N)=Y2(N-1)
        
      RETURN 
      END 

c	**********************************************	  

	subroutine quadfit(n,a,b,np,c)
	integer n,nmax,i,j,k,np,npp
	parameter(nmax=50)
	real*8 a(n),b(n),c(n),s(0:nmax,0:nmax),
     .	x(nmax,nmax),y(nmax,1)


	if(n.gt.nmax) then
	  write(6,*) 'In quadfit, array limits exceeded'
	  return
	endif
	
	do i=0,2*np
	  do j=0,2*np
	    s(i,j)=0.d0
	  end do
	end do

	do i=0,2*np
	  do j=0,1
	     do k=1,n
	       s(i,j)=s(i,j)+a(k)**i*b(k)**j
	     end do
	  end do
	end do

	npp=np+1
	do i=npp,1,-1
	  y(i,1)=s(npp-i,1)
	end do

	do i=npp,1,-1
	  do j=npp,1,-1
	    x(i,j)=s(npp-j+npp-i,0)
	  end do
	end do

	call gaussj(x,npp,nmax,y,1,1)

	do i=1,npp
c	  write(6,*) 'x^',npp-i,' ',  y(i,1)
	  c(i)=y(i,1)
	end do

	return
	end

c	************************************************************

	SUBROUTINE gaussj(a,n,np,b,m,mp)
	INTEGER m,mp,n,np,NMAX
	REAL*8 a(np,np),b(np,mp)
	PARAMETER (NMAX=50)
c	    Linear  equation solution by Gauss-Jordan elimination
c	    a(1:n,1:n) is an input matrix stored in an array of
c	    physical dimensions np by np. b(1:n,1:m) is an input 
c	    matrix containing the right-hand side vectors, stored
c	    in an array of physical dimensions np by mp. On output,
c	    a(1:n,1:n) is replaced by its matrix inverse, and 
c	    b(1:n,1:m) is replaced by the corresponding set of 
c	    solution vectors.
c	    Parameter: NMAX is the largest anticipated value of n.
	INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),
     .		ipiv(NMAX)	! The integer arrays ipiv, indxr
     				! and indxc are used for bookkeeping
	REAL*8 big,pivinv,dum	! on the pivoting.
	
	do 11 j=1,n
	  ipiv(j)=0
 11	enddo 
	do 22 i=1,n
	  big=0.d0
	  do 13 j=1,n
	     if(ipiv(j).ne.1) then
	       do 12 k=1,n
	          if(ipiv(k).eq.0) then
		    if(dabs(a(j,k)).ge.big)then
		      big=dabs(a(j,k))
		      irow=j
		      icol=k
		    endif
		  else if (ipiv(k).gt.1) then
		    pause 'singular matrix in gaussj'
		  endif
 12	       enddo 
	     endif
 13	  enddo 
	  ipiv(icol)=ipiv(icol)+1
c	  We now have the pivotal element, so we interchane rows, 
c	  if needed, to put the pivotal element on the diagonal. The 
c	  columns are not physically interchanged, only relabelled;
c	  indxc(i), the column of the ith pivot element, is the ith
c	  column that is reduced, while indxr(i) is the row in which 
c	  that pivot element was originally located. If indxr(i).ne.
c	  indxc(i), there is an implied column interchane. With this 
c	  of bookkeeping, the solution b's will end up in the correct 
c	  order, and the inverse matrix will be scrambled by columns.

	  if(irow.ne.icol)then
	      do 14 l=1,n
	         dum=a(irow,l)
		 a(irow,l)=a(icol,l)
		 a(icol,l)=dum
 14	      enddo 
	      do 15 l=1,m
	         dum=b(irow,l)
		 b(irow,l)=b(icol,l)
		 b(icol,l)=dum
 15	      enddo 
	  endif
	  indxr(i)=irow		! We are now ready to divide the pivot
	  indxc(i)=icol		! row by the pivot element located at 
	  			! irow and icol.
	  if(a(icol,icol).eq.0.d0) pause 'singular matrix in gaussj'
	  pivinv=1.d0/a(icol,icol)
	  a(icol,icol)=1.d0
	  do 16 l=1,n
	       a(icol,l)=a(icol,l)*pivinv
 16	  enddo 
	  do 17 l=1,m
	       b(icol,l)=b(icol,l)*pivinv
 17	  enddo 
	  
	  do 21 ll=1,n			! Next, we reduce the rows ..
	      if(ll.ne.icol) then	! .. except for the pivot one.
	         dum=a(ll,icol)
		 a(ll,icol)=0.d0
		 do 18 l=1,n
		      a(ll,l)=a(ll,l)-a(icol,l)*dum
 18		 enddo 
		 do 19 l=1,m
		      b(ll,l)=b(ll,l)-b(icol,l)*dum
 19		 enddo 
	      endif
 21	  enddo 
 22	enddo 	! This is the end of the main loop over
			! columns of the reduction.
	do 24 l=n,1,-1	! It only remains to unscramble the solution
  	  if(indxr(l).ne.indxc(l))then	! in view of the column interchanges.
	    do 23 k=1,n			! We do this by interchanging pairs
	      dum=a(k,indxr(l))		! of columns in the reverse order that 
	      a(k,indxr(l))=a(k,indxc(l))  ! the permutation was built up.
	      a(k,indxc(l))=dum
 23	    enddo 
	  endif
 24	enddo 
	return
	END
		
c	**********************************************	  
