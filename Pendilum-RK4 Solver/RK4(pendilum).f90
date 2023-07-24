!----------------------------------------------------------------
!..A RK4 SOLVER for Pendilum Problem
!----------------------------------------------------------------
Module data
  integer,parameter   :: neq=2
  real, parameter :: gl = 2., kml = 3.
  real,dimension(neq) :: theta0=(/75.,0./)  !..IC
End module

Program RK4
  use data
  real,dimension(neq) :: theta
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
   theta     = theta0
   write(1,'(3e15.7)') time, (theta(n),n=1,neq)
   !..solution loop
          DO WHILE (time .lt. tmax)
           call SRK4(time,dt,theta)
           time = time + dt
           write(1,'(3e15.7)') time, (theta(n),n=1,neq)
          ENDDO

   !..Close the output file
   close(1)
stop
End
!----------------------------------------------------------------------
       Subroutine SRK4(time,dt,theta)
       use data
       real,dimension(neq) :: theta,t_tmp,k1,k2,k3,k4,k

       dt2 = 0.5*dt
       call ODES(time,theta,k1)
       t_tmp(:) = theta(:) + dt2*k1(:)
       call ODES(time,t_tmp,k2)
       t_tmp(:) = t_tmp(:) + dt2*k2(:)
       call ODES(time,t_tmp,k3)
       t_tmp(:) = t_tmp(:) + dt*k3(:)
       call ODES(time,t_tmp,k4)
       k(:) = (k1(:)+2.*k2(:)+2.*k3(:)+k4(:))/6.
       theta(:) = theta(:) + dt*k(:)

       return
       End
!--------------------------------------------------
Subroutine ODES(time,theta,f)         !..Define the ODEs
  use data
  real,dimension(neq) :: theta,f
  f(1) = theta(2)
 ! f(2) = -gl*sin(theta(1))-kml*theta(2)  !Non-linear
  f(2) = -gl*(theta(1))-kml*theta(2)      !Linearized
return
End
