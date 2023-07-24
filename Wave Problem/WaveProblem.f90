!-------------------------------------------------------------------------
!..AN EXPLICIT FD SOLVER FOR THE 1-D DIFFUSION/LINEAR CONVECTION EQUATION  
!-------------------------------------------------------------------------
module data
  integer, parameter :: imax=201, ntout=1
  integer :: ntmax
  real, dimension(imax) :: wave_n, wave_np1, xg
  real :: dnum,sigma,diff_number,sigma2, dx=0.1, x1, pi=ACOS(-1.),dt=1,alpha=1.0e-3,V=1

end module data

program EXPLICIT_FDE
   use data

   call INIT             !..Read the input data, and initialize the wave
   DO nt = 1,ntmax       !..Start the solution loop 
      do i = 2,imax-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !..diffusion
     wave_np1(i) = wave_n(i) + dnum*(wave_n(i+1) - 2.*wave_n(i) + wave_n(i-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !..convection(FTCS)
!     wave_np1(i) = wave_n(i) - 0.5*sigma*(wave_n(i+1) - wave_n(i-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !..convection(FTBS)
!     wave_np1(i) = wave_n(i) - sigma*(wave_n(i) - wave_n(i-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !..Convection-Diffusion(FTCS)
!      wave_np1(i) = wave_n(i)-(sigma/2)*(wave_n(i+1)-wave_n(i-1))+&
!      diff_number*(wave_n(i+1)-2*wave_n(i)+wave_n(i-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !..Convection-Diffusion(FTCS+FTBS(only for convective term))
!      wave_np1(i) = wave_n(i)-sigma*(wave_n(i)-wave_n(i-1))+diff_number*&
!     (wave_n(i+1)-2*wave_n(i)+wave_n(i-1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo
      wave_n = wave_np1
!..Output intermediate solutions
      if( MOD(nt,ntout) .eq. 0 .or. nt .eq. ntmax ) call IO(nt)
   ENDDO                        

stop
end program EXPLICIT_FDE
!------------------------------------------------------------------------
subroutine INIT
  use data

  write(*,'(/(a))',advance='no')'  Enter dnum/sigma and ntmax : '
  read(*,*) dnum, ntmax
  diff_number = alpha*dt/(dx**2)     !values for question3
  sigma2 = V*dt/dx                   !values for question3
  sigma = dnum
  x = -(imax-1)*dx/2.
  do i = 1,imax             !..Initialize the wave 
    xg(i) = x
  ! if( x .gt. -1. .and. x .lt. 1. ) then
   !   wave_n(i) = 2*SIN(pi*x)                    !..sin wave
   !endif
    if( x .gt. -1. .and. x .lt. 0. ) then
       wave_n(i) = 1. + x                       !..triangular wave
    else if( x .gt. 0. .and. x .lt. 1. ) then
       wave_n(i) = 1. - x
    else
       wave_n(i) = 0.
    endif
    x = x + dx
  enddo
  call IO(0)

return 
end subroutine INIT
!-------------------------------------------------------------------
subroutine IO(nt)
   use data
   character :: fname*32,string*6,ext*3
   write(string,'(f5.3)') float(nt)/1000
   read(string,'(2x,a3)') ext
   fname = 'wave-'//ext//'.txt'
   open(1,file=fname,form='formatted')
   write(1,'(2e14.6)') (xg(i), wave_n(i), i=1,imax)
   close(1)
return 
end subroutine IO

