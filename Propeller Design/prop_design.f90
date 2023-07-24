Module data
  real, parameter :: pi=4*atan(1.), rho=1.225 !density[kg/m^3]
  integer :: n
  real :: T0, n_b, U_inf, U_l, C_l, alpha_cld, Cd_Cl, n_i
  real :: R, omega, beta, U_s, U, fi, c, dT, dM, n_l, b, x, y, z, root1, root2, U_R, M, sigma
End module

program Design
   use data
   
   character*40 fname
   
   print*,"Please enter the values by order as: || T, n_b, U_inf, U_l, C_l, alpha_cld, Cd_Cl, n_i, n ||"
   read*,T, n_b, U_inf, U_l, C_l, alpha_cld, Cd_Cl, n_i, n
   
   write(*,'(/,(a))',advance='no')' Enter Output File Name [solution.dat] :> '
   read(*,'(a)') fname
   if( fname .eq. ' ') fname = 'solution.dat'
   open(1,file=fname,form='formatted')

   U = U_inf/n_i  !eqn1
   U_s = ((2/n_i)-1)*U_inf  !eqn2
   R = sqrt(((2*T)/(rho*pi*((U_inf)**2)))*((n_i**2)/(4-4*n_i)))  !eqn3
   omega = (sqrt(U_l**2-U**2))/R  !eqn4
   T=T0
   dr = R/n
   M = 0.
   write(1,100)
   do i=3,n      !first two values are omitted
      r1 = i*dr
      beta = U/(omega*r1)
      ! finding roots of xb^2+yb+z=0
        x=1
        y=beta*Cd_Cl*n_i-1
        z=beta*(beta+Cd_Cl)*(1-n_i)

        if (y**2-4*x*z .ge. 0.) then
        root1 = (-y+sqrt(y**2-4*x*z))/(2*x)
        root2 = (-y-sqrt(y**2-4*x*z))/(2*x)
        else
        continue
        endif
        
        b=root1
        if (root2 .lt. root1) then
        b=root2
        endif
      fi = asin(beta/(sqrt(((1-b)**2)+beta**2)))  !eqn8
      theta = (alpha_cld + fi)*180/pi  !eqn7
      U_R = U/sin(fi)
      c = ((8*pi*r1)/(n_b*C_l))*(1-(U_inf/U))*((sin(fi)**2)/(cos(fi+Cd_Cl)))  !eqn5
      sigma = (n_b*c)/(2*pi*r1)
      dT = pi*r1*sigma*rho*(U_R**2)*C_l*cos(fi+Cd_Cl)*dr  !eqn9
      dM = pi*(r1**2)*sigma*rho*(U_R**2)*C_l*sin(fi+Cd_Cl)*dr  !eqn10
      n_l = (U_inf/(omega*r1))*(cos(fi+Cd_Cl)/sin(fi+Cd_Cl))  !eqn11
      
      T = T + dT
      M = M + dM
      
    Tf = T - T0
    P_h = omega*M  !eqn13
    P_f = Tf*U_inf  !eqn14
    eff = P_f/P_h  !eqn15

   write(1,101) r1,c,theta,dT,dM,n_l,sigma

   enddo
   
   write(1,102)
   write(1,103) Tf,M,P_h,P_f,eff
   
   close(1)
   
  100 format ('     R(m)          C(m)      Theta(deg)     dT(N)        dM(Nm)       ETA_L       Sigma')
  101 format (7(1x,e12.5))
  102 format ('       T            M           PH           PF           ETA')
  103 format (5(1x,e12.5))
   

 stop
end
!     2500 2 80 320 0.4 0 0.015 0.96 20
!     200 2 50 150 0.4 0 0.05 0.99 20
!     110 2 20 102.087 1.58 6.5 0.01187147 0.96 20
