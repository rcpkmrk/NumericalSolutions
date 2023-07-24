!------------------------------------------------------------------------
!  RK4 SOLVER for a system of ODEs                                      |
!------------------------------------------------------------------------
Module data
   integer,parameter  :: neq=2
   real,parameter     :: x=0.1, L=1.,&
                         Vinf=20., Tinf=300., Twall=500., &
                         nu=5.2e-5, Pr=0.72, &
                         Re  = L*Vinf/nu, &
                         Rex = (x * Vinf) / nu, &
                         Cf  = (0.664) * sqrt(Rex), &
                         Beta= Cf * Pr * sqrt(Re) / (2**(1.5))
   real,dimension(neq):: T0=(/Twall, 0./)         !..ICs
End module

PROGRAM sysRK4
   use data   
   real,dimension(neq) :: T
   character*40 fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                  Constant Stepsize Algorithm                     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !..read the input data
   write(*,'(/,(a))',advance='no')' Enter d_eta and  eta_max :> '
   read(*,*) deta,etamax
   write(*,'(/,(a))',advance='no')' Enter DT/d_eta(0)  :> '
   read(*,*) T0(2)
   write(*,'(/,(a))',advance='no')' Enter Output File Name [solution.dat] :> '
   read(*,'(a)') fname
   if( fname .eq. ' ') fname = 'solution.dat'
!..write the header and the initial conditions into output file
   open(1,file=fname,form='formatted')

   eta = 0.
   T   = T0
   write(1,'(5E15.7)') eta,(T(n),n=1,neq)
!..solution loop

   DO WHILE (eta .lt. etamax)
     call SRK4(eta,deta,T)
     eta = eta + deta
     ! y = eta * sqrt((2 * nu * x)/Vinf)   for 2nd question
     write(1,'(6E15.5)') eta,(T(n),n=1,neq)
   ENDDO

  close(1)
stop
End
!----------------------------------------------------------------------
Subroutine SRK4(eta,deta,T)
   use data   
   real,dimension(neq) :: T,Ttmp,k1,k2,k3,k4,k
   de2 = 0.5*deta
   call ODE(eta,T,k1)
   Ttmp(:) = T(:) + de2*k1(:)
   call ODE(eta,Ttmp,k2)
   Ttmp(:) = Ttmp(:) + de2*k2(:)
   call ODE(eta,Ttmp,k3)
   Ttmp(:) = Ttmp(:) + deta*k3(:)
   call ODE(eta,Ttmp,k4)
   k(:) = (k1(:)+2.*k2(:)+2.*k3(:)+k4(:))/6.
   T(:) = T(:) + deta * k(:)
return
End

!------------------------------------------------------------------------
Subroutine ODE(eta,T,f)              !..define the ODE's 
   use data   
   real,dimension(neq) :: T,f
   f(1) = T(2)
   f(2) = (-1) * Beta * T(2) *  eta**2
return
End 
