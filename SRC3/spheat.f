	Subroutine intener() 

	include 'Glob_decl'
	real*8 H11,H22,H33
	real*8 H1,H2
        complex*16 gamma,gauss,sem
	
	external fermic,gauss,sem
      
	include 'Glob_cons'


        H11=0.d0
	H1=0.d0
	H2=0.d0
        

         
        do i=-N,N
         r2=-2.d0*w(i)*dimag(Gfcpa(i)
     .       +Gccpa(i))/pi
	 r3=dimag(ssig(i)*Gfcpa(i))/pi
         H1=H1+dw(i)*r2*ferm(i)
	 H2=H2+dw(i)*r3*ferm(i)
        end do
         H11=H1+H2  	
	 
	open(unit=31,file='intener.dat',status='unknown')
        write(31,*) frac,H11
	close(31)
        write(6,*) frac,H11


	return
	
	end

c	*************************************************************
