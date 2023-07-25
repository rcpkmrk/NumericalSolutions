Module data
  integer, parameter :: imax=201,jmax=101,kres=10, ksol=1000
  integer :: ibs,ibe,jbs,jbe,kmax
  real, parameter :: xl=10., yl=5., omega=1., rl2allow=-5.
  real, dimension(imax,jmax) :: phi_k,phi_kp1,u,v
  real :: dx,dy,beta2,dx2i,dy2i,uinf, rl2=1., x(imax), y(jmax)
  logical, dimension(imax,jmax) :: flow
End module

!------------------------------------------------------------------------------|
!..A POINT/LINE ITERATIVE SOLVER FOR ELLIPTIC PDEs                             |
!------------------------------------------------------------------------------|
 program ELLIPTIC
   use data
!..Read the input data, generate the grid data and initialize the solution
  k = 1
  call INIT
!..Start the iterative solution loop
  DO WHILE (k .lt. kmax .and. rl2 .gt. rl2allow )
     k = k+1
     call POINT_ITERATE           !..Point iterative solutions
!    call LINE_ITERATE            !..Line iterative solutions
     call BC                      !..Apply BCs 
!..Update phi_kp1 array and evaluate the L2 norm of the residual
     rl2 = SQRT(SUM((phi_kp1 - phi_k)**2))
     if(k .eq. 2) rl2_1=rl2
     rl2 = ALOG10(rl2/rl2_1)
     if( MOD(k,kres).eq.0 .or. k .eq. 2) print*,k,rl2  !       'Residual @ k =',
     phi_k = phi_kp1
!..Output intermediate solutions
     if( MOD(k,ksol) .eq. 0 .and. k .ne. kmax) call QOUT(k)
  ENDDO
  !print*, 'Residual @ k =',k,rl2
  read*,a
  call QOUT(k)

stop
end
!------------------------------------------------------------------------
subroutine INIT
 use data

  write(*,'(a)',advance='no') 'Enter kmax: '
  read(*,*) kmax

  dx = xl/(imax-1)    !dx=0.05
  dy = yl/(jmax-1)    !dy=0.05
  beta2 = (dx/dy)**2
  dx2i = 1./(2*dx)
  dy2i = 1./(2*dy)
  uinf = 1.

!..Set wall/object boundary locations and body grid points
  ibs = 4/dx + 1       !ibs=81
  ibe = 6/dx + 1       !ibe=121
  jbs = 2/dy + 1       !jbs=41
  jbe = 3/dy + 1       !jbe=61
  flow = .true.
  flow( ibs:ibe, jbs:jbe ) = .false.

!..Grid generation 
  do i=1,imax
     x(i)= dx*(i-1) 
  enddo
  do j=1,jmax
     y(j)= dy*(j-1) 
  enddo

  phi_k   = 1.       !..Initial guess, a bad one
  phi_kp1 = 1.
  call BC
  call QOUT(1)

return 
end
!-----------------------------------------------------------------
subroutine POINT_ITERATE
!..Implement  Point Jacobi, Gauss-Seidel and SOR methods
  use data
  cm  = 0.5/(1.+beta2)
  do j = 2,jmax-1
  do i = 2,imax-1
     !..solve for phi_kp1 in the flow domain
!     phi_kp1(i,j) = cm*(phi_k(i-1,j)+phi_k(i+1,j)+beta2*(phi_k(i,j-1)+phi_k(i,j+1)))  !Point Jacobi
     
 !     phi_kp1(i,j) = cm*(phi_kp1(i-1,j)+phi_k(i+1,j)+beta2*(phi_kp1(i,j-1)+phi_k(i,j+1)))  !Gauss-Seidel

     phi_kp1(i,j) = (1-omega)*phi_k(i,j) + (omega*cm)*(phi_kp1(i-1,j)+phi_k(i+1,j)+beta2*(phi_kp1(i,j-1)+phi_k(i,j+1)))  !SOR
  enddo
  enddo

return 
end
!-------------------------------------------------------------------
subroutine BC
 use data
  do j = 2,jmax-1  !..Set left/right farfield BC
     phi_k  (1,j)    = phi_k(2,j) - uinf*dx
     phi_kp1(1,j)    = phi_k(1,j)
  enddo

  do i = 1,imax    !..Set bottom/top farfield BC
     phi_k(i,1)      = phi_k(i,2)    
     phi_k(i,jmax)   = phi_k(i,jmax-1) 
     phi_kp1(i,1)    = phi_k(i,1)
     phi_kp1(i,jmax) = phi_k(i,jmax)
  enddo

!..Set wall BCs around the object
  do i = ibs,ibe
     phi_k(i,jbs)   = phi_k(i,jbs-1)
     phi_k(i,jbe)   = phi_k(i,jbe+1)
     phi_kp1(i,jbs) = phi_k(i,jbs)
     phi_kp1(i,jbe) = phi_k(i,jbe)
  enddo
  do j = jbs,jbe
     phi_k(ibs,j)   = phi_k(ibs-1,j)
     phi_k(ibe,j)   = phi_k(ibe+1,j)
     phi_kp1(ibs,j) = phi_k(ibs,j)
     phi_kp1(ibe,j) = phi_k(ibe,j)
  enddo

return 
end
!-------------------------------------------------------------------
subroutine VELOCITY
  use data
!..Compute u and v velocity components & enforce the FS/wall velocities
  u = 0.
  v = 0.
  u(1,:) = uinf
  u(imax,:) = (phi_k(imax,:)- phi_k(imax-1,:))/dx
  do i = 2,imax-1
  u(i,1)    = (phi_k(i+1,1) - phi_k(i-1,1))*dx2i
  u(i,jmax) = (phi_k(i+1,jmax) - phi_k(i-1,jmax))*dx2i
  do j = 2,jmax-1
     if( flow(i,j) ) then
       u(i,j) = (phi_k(i+1,j)- phi_k(i-1,j))*dx2i
       v(i,j) = (phi_k(i,j+1)- phi_k(i,j-1))*dy2i
     endif
  enddo
  enddo
return 
end
!-------------------------------------------------------------------
subroutine QOUT(k)          !..Output the solution
  use data
  character fname*32,string*9,ext*6
  write(string,'(f9.6)') float(k)/1000000
  read(string,'(3x,a6)') ext
  fname = 'q-'//ext//'.plt'
  open(1,file=fname,form='formatted')
  write(1,*) ' variables="x","y","phi","u","v","k","r12"'
  write(1,*) ' zone i=',imax, ', j=',jmax
  call VELOCITY()
  do j=1,jmax
  do i=1,imax
     write(1,*) x(i),y(j),phi_k(i,j),u(i,j),v(i,j),k,r12
  enddo
  enddo
  close(1)
return
end
!-------------------------------------------------------------------
subroutine THOMAS(is,ie, a,b,c,f)
!
!  Solution of a tridiagonal system of n equations of the form
!  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = F(i)  for k=is,ie
!  the solution X(i) is stored in F(i)
!  A(is-1) and C(ie+1) are not used.
!  Solution is returned in array F
!
   use data, only : imax
   real, dimension(imax) :: a,b,c,f
   do i=is+1,ie
      w    = a(i)/b(i-1)
      b(i) = b(i) - w*c(i-1)
      f(i) = f(i) - w*f(i-1)
   enddo
   f(ie) = f(ie)/b(ie)
   do i=ie-1,is,-1
      f(i) = (f(i) - c(i)*f(i+1))/b(i)
   enddo
return
end
