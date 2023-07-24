!------------------------------------------------------------------------
!  RK4 SOLVER for a system of ODEs                                      |
!------------------------------------------------------------------------
Module data
   integer,parameter   :: neq=4
   real,parameter      :: k_z=20E3, c_z=300., k_alpha= 15E3, rho_inf= 1.225, pi= 4.d0*atan(1.d0), &
   			   S= 15., M= 150., i= 50., l= .15, u_inf= 40., c_alpha= 150.
   real,dimension(neq) :: q0=(/0.,0.,8.*pi/180.,0./)     !..IC
End module

PROGRAM sysRK4
   use data   
   real,dimension(neq) :: q, q_degree
   character*40 fname
 
!..read the input data
   write(*,'(/,(a))',advance='no')' Enter TimeStep and FinalTime :> '
   read(*,*) dt,tmax
   write(*,'(/,(a))',advance='no')' Enter Output File Name [solution.dat] :> '
   read(*,'(a)') fname
   if( fname .eq. ' ') fname = 'solution.dat'
   open(1,file=fname,form='formatted')

   time = 0.
   q    = q0
   call CL(q,c_l)
   !..Write out ICs
   write(1,'(6E15.5)') time,(q(n),n=1,2),q(3)*180/pi,q(4),c_l
   !..solution loop
   DO WHILE (time .lt. tmax)
     call SRK4(time,dt,q)
     q_degree(:) = q(:)*180/pi
     time = time + dt 
     call CL(q,c_l)
     write(1,'(6E15.5)') time,(q(n),n=1,2),(q_degree(k),k=3,4),c_l
   ENDDO
   close(1)

stop
End
!----------------------------------------------------------------------
Subroutine SRK4(time,dt,q)
   use data   
   real,dimension(neq) :: q,qtmp,k1,k2,k3,k4
       
   dt2 = 0.5*dt
   call ODES(time,q,k1)
   qtmp(:) = q(:) + dt2*k1(:)
   call ODES(time,qtmp,k2)
   qtmp(:) = qtmp(:) + dt2*k2(:)
   call ODES(time,qtmp,k3)
   qtmp(:) = qtmp(:) + dt*k3(:)
   call ODES(time,qtmp,k4)
   q(:) = q(:) + dt*((k1(:)+2.*k2(:)+2.*k3(:)+k4(:))/6.)
return
End
!------------------------------------------------------------------------
Subroutine CL(q,c_l)              !..define c_l
   use data
   real,dimension(neq) :: q
   
   if(q(3).lt.12.*pi/180.) then
   	c_l = 2*pi*q(3)
   else if(q(3).lt.16.*pi/180.) then
   	c_l = 2*pi*12.*pi/180.
   else if(q(3).lt.19.*pi/180.) then
   	c_l = 2*pi*12.*pi/180.+2*pi*(16.*pi/180.-q(3))
   else
   	stop
   end if
return
End
!------------------------------------------------------------------------
Subroutine ODEs(time,q,f)              !..define the ODE's 
   use data   
   real,dimension(neq) :: q,f
   
   call CL(q,c_l)

   f(1) = q(2)
   f(2) = (((0.5*rho_inf*u_inf*S*c_l*sqrt(u_inf**2.+q(2)**2.))) -  &
   c_z*q(2) - k_z*q(1))/M
   f(3) = q(4)
   f(4) = (((0.5*rho_inf*u_inf*S*c_l*l*sqrt(u_inf**2.+q(2)**2.)**.5)) &
   -  c_alpha*q(4) - k_alpha*q(3))/i
return
End 
