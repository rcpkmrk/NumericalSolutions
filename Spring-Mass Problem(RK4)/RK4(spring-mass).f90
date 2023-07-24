!----------------------------------------------------------------
!..A RK4 SOLVER for Spring-mass System Problem
!----------------------------------------------------------------
Module data
  integer,parameter   :: neq=4
  real, parameter :: k = -5., m1 = 1., m2 = 2. !spring constant and masses
  real,dimension(neq) :: x0=(/1.,0.,1.,0./)  !..IC
End module

Program RK4
  use data
  real,dimension(neq) :: x
  character*40 fname
!..Read the stepsize
  write(*,'(/,(a))',advance='no')' Enter TimeStep and FinalTime :> '
        read(*,*) dt,tmax
        write(*,'(/,(a))',advance='no')' Enter Output File Name [solution.dat] :> '
        read(*,'(a)') fname
        if( fname .eq. ' ') fname = 'solution.dat'
        !..write the header and the initial conditions into output file
         open(1,file=fname,form='formatted')

!..Set the Initial Conditions and output them
   time      = 0.
   x       = x0
   write(1,'(5E15.5)') time, (x(n),n=1,neq)
   !..solution loop
          DO WHILE (time .lt. tmax)
           call SRK4(dt,time,x)
           time = time + dt
           write(1,'(5E15.5)') time, (x(n),n=1,neq)
          ENDDO
   !..Close the output file
   close(1)
stop
End
!----------------------------------------------------------------------
       Subroutine SRK4(dt,time,x)
       use data
       real,dimension(neq) :: x,xtmp,k1,k2,k3,k4,kf

       dt2 = 0.5*dt
       call ODES(time,x,k1)
       xtmp(:) = x(:) + dt2*k1(:)
       call ODES(time,xtmp,k2)
       xtmp(:) = xtmp(:) + dt2*k2(:)
       call ODES(time,xtmp,k3)
       xtmp(:) = xtmp(:) + dt*k3(:)
       call ODES(time,xtmp,k4)
       kf(:) = (k1(:)+2.*k2(:)+2.*k3(:)+k4(:))/6.
       x(:) = x(:) + dt*kf(:)

       return
       End
!--------------------------------------------------
Subroutine ODES(time,x,f)         !..Define the ODEs
  use data
  real,dimension(neq) :: x,f
  f(1) = x(2)
  f(2) = (2.*k/m1)*x(1)-(k/m1)*x(3)
  f(3) = x(4)
  f(4) = -(k/m2)*x(1)+2*(k/m2)*x(3)
return
End
