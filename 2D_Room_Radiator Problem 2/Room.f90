Module data
  integer, parameter :: imax=201, jmax=121, outfreq=1000
  real, parameter    :: xl=10., yl=6., omega=1.5, rl2tol=-5., &
  k_coeff=0.025
  real :: dx,dy,beta2,c1,rl2=1.
  real :: x(imax),y(jmax),qk(imax,jmax),qkp1(imax,jmax), &
  qx(imax,jmax), qy(imax,jmax)
  integer :: irs,ire,jrs,jre, kmax, k=1
  logical :: domain(imax,jmax)
End module

!------------------------------------------------------------------------------|
!..A ITERATIVE/DIRECT SOLVER FOR ELLIPTIC PDEs                                 |
!------------------------------------------------------------------------------|
program ELLIPTIC
   use data

!..Read the input data, generate the grid data and initialize the solution
   call INIT 
!..Start the iterative solution loop 
   DO while ( k .lt. kmax .and. rl2 .gt. rl2tol )
      k = k+1
!     call GS         !..Point iterative solutions; Jacobi, GS, SOR
      call LGS        !..Line iterative solution; LGS,SLOR, ADI
      call BC         !..Apply BC
!..Evaluate L2 norm of the the residual
      rl2 = SQRT(SUM((qkp1 - qk)**2))
      if(k .eq. 2) rl2_1=rl2
      rl2 = ALOG10(rl2/rl2_1)
      print*, ' Norm. Log of residual @ k =',k,rl2
      write(9,*) k, rl2
      
!..Output intermediate solutions
      if( MOD(k,outfreq).eq.0 .and. k.ne.kmax) call QOUT
      qk = qkp1
   ENDDO
!question 3 outputs
      do i=1,imax
         do j=1,jmax
            if (x(i) .eq. 5.) then
               write(12,*) y(j),qk(i,j)
            endif
         enddo
      enddo
      
      do j=1,jmax
         do i=1,imax
            if (y(j) .eq. 3.) then
               write(15,*) x(i),qk(i,j)
            endif
         enddo
      enddo
   
   call QOUT
return 
end

!------------------------------------------------------------------------
subroutine INIT
 use data

  write(*,'(a)',advance='no') 'Enter kmax: '
  read(*,*) kmax

  dx = xl/(imax-1)
  dy = yl/(jmax-1)
  beta2 = (dx/dy)**2
  c1    = 0.5/(1.+beta2)

!..Set the radiator location and solution domain
  irs = 2/dx + 1
  ire = 8/dx + 1
  jrs = 5.6/dy + 1
  jre = 5.8/dy + 1
  domain = .true.
  domain( irs:ire, jrs:jre ) = .false.

!..Grid generation 
  do i=1,imax
     x(i)= dx*(i-1) 
  enddo
  do j=1,jmax
     y(j)= dy*(j-1) 
  enddo

  q_k = 0.       !..Initial guess
  call BC
  call QOUT
  open(9,file='rl2.dat',form='formatted')
  open(12,file='T vs. y(x=5).dat',form='formatted')
  open(15,file='T vs. x(y=3).dat',form='formatted')
  

return 
end

!-------------------------------------------------------------------
subroutine BC
  use data
!..Set the BCs
  qk(:,1)      = 20.
  qkp1(:,1)    = 20.
  qk  (imax,:) = 20.
  qkp1(imax,:) = 20.
  qk(:,jmax)   = 0.
  qkp1(:,jmax) = 0.
  qk  (1,:)    = qk(2,:)
  qkp1(1,:)    = qk(2,:)

!..Set BCs for the radiator
 qk( irs:ire, jrs:jre )   = 50
 qkp1( irs:ire, jrs:jre ) = 50

return 
end

!-------------------------------------------------------------------
subroutine GS   !..Jacobi->GS->SOR
   use data
   do j = 2,jmax-1
   do i = 2,imax-1
!..Solve for q^k+1
     if( domain(i,j) ) then 
        !qkp1(i,j) = c1*(qkp1(i-1,j)+qk(i+1,j)+beta2*(qkp1(i,j-1)+qk(i,j+1)))  !GS
        
        !qkp1(i,j) = c1*(qk(i-1,j)+qk(i+1,j)+beta2*(qk(i,j-1)+qk(i,j+1)))  !Point Jacobi
        
        qkp1(i,j) = (1-omega)*qk(i,j) + (omega*c1)*(qkp1(i-1,j)+qk(i+1,j)+beta2*(qkp1(i,j-1)+qk(i,j+1))) !SOR
 
     endif
   enddo
   enddo
return 
end

!-------------------------------------------------------------------
subroutine LGS
   use data
   implicit none
   integer :: i, j

  !..Line Iterative Method
   real, dimension(imax) :: cl, cm, cu, rhs
  do j =2, jmax-1
     do i = 2, imax-1
        cl(i) = 1.
        cm(i) =  -2. * (1. + beta2)
        cu(i) = 1.
        rhs(i) = -beta2*(qk(i,j+1) + qkp1(i,j-1))
     enddo
     rhs(2)      = rhs(2) -  qkp1(1,j)
     rhs(imax-1) = rhs(imax-1) -  qkp1(imax,j)
     
     call THOMAS(2,imax-1, cu,cm,cl,rhs)   !calling thomas with proper values

     qkp1(2:imax-1,j) = rhs(2:imax-1)    !the new value which is equal to the one obtained by thomas

  enddo
   return
  end

!-------------------------------------------------------------------
subroutine FLUX()
  use data
!..Compute the heat flux vectors
  do j = 2,jmax-1
  do i = 2,imax-1
       qx(i,j) = -k_coeff*((qk(i+1,j)-qk(i-1,j))/(2*dx))
       qy(i,j) = -k_coeff*((qk(i,j+1)-qk(i,j-1))/(2*dy))
  enddo
  enddo
 return
end

!-------------------------------------------------------------------
subroutine QOUT          !..Output the solution in tecplot format
   use data
   character :: fname*32,string*8,ext*5
   write(string,'(f8.5)') float(k)/100000
   read(string,'(3x,a5)') ext
   fname = 'q-'//ext//'.plt' 
   call FLUX()
   open(1,file=fname,form='formatted')
   write(1,*) ' variables="x","y","q","qx","qy" '
   write(1,*) ' zone i=',imax, ', j=',jmax
   write(1,'(5e12.4)')((x(i),y(j),qk(i,j),qx(i,j),qy(i,j),i=1,imax),j=1,jmax)
   close(1)
return
end

!-------------------------------------------------------------------
subroutine THOMAS(is,ie, a,b,c,f)
   use data, only : imax
   implicit none
   integer :: is, ie, isp1, i, ii, iepis
   real :: z
   real, dimension(imax) ::  a,b,c,f,x
!
!  Solution of a tridiagonal system of n equations of the form
!  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = F(i)  for k=is,ie
!  the solution X(i) is stored in F(i)
!  A(is-1) and C(ie+1) are not used.
!  Solution is returned in array F
!
   x(is)=c(is)/b(is)
   f(is)=f(is)/b(is)
   isp1 = is+1
   do i=isp1,ie
      z   =1./(b(i)-a(i)*x(i-1))
      x(i)=z*c(i)
      f(i)=(f(i)-a(i)*f(i-1))*z
   enddo
   iepis=ie+is
   do ii=isp1,ie
      i=iepis-ii
      f(i)=f(i)-x(i)*f(i+1)
   enddo
return
end

!-------------------------------------------------------------------
subroutine GAUSS(n,A,B)
   real :: a(n,n),b(n)
!..Convert to upper triangular form
   do k = 1,n-1
      if (abs(a(k,k)).gt.1.e-6) theN
      do i = k+1, n
         x = a(i,k)/a(k,k)
         do j = k+1, n
            a(i,j) = a(i,j) -a(k,j)*x
         enddo
         b(i) = b(i) - b(k)*x
      enddo
      else
        write (6,*) 'zero pivot found in line:',k
        stop
      endif
   enddo
!..Back substitution
   do i = n,1,-1
      rsum = b(i)
      if (i.lt.n) then
         do j= i+1,n
            rsum = rsum - a(i,j)*b(j)
         enddO
      endif
      b(i) = rsum/a(i,i)
   enddo
return
end
