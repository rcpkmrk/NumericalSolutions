Module ac_data
  implicit none
  real*8, parameter :: mass = 15000, &       ! mass of a/c in kg
                       wing_area = 60, &     ! wing planform area of a/c in m^2
                       CL0 = 0.1896, &           ! zero-AOA lift coeff.
                       CD0 = 0.025d0, &         ! zero-lift drag coeff.
                       dCL_dalpha = 5.52d0, &   ! lift coeff. curve slope
                       K2 = 0.05d0              ! parasite drag coeff. in drag polar:
                                                ! CD = CD0 + K2 * CL^2
End module

Module flight_data
  implicit none
  real*8, parameter :: Vground = 150, &      ! ground speed in m/s
                       h_alt0 = 5000.d0,  &      ! initial altitude
                       grav = 9.81d0, &         ! gravitiy
                       Vamp_turb = 10, &     ! ghust amplitude in m/s
                       kturb = 0.01            ! wavno of ghust velocity m^-1

  real*8 :: alpha, &    ! instantaneous AOA in radians
            CL, &       !            "  lift coef.
            CD, &       !            "  drag coef.
            rho_inf, &  !            "  air density
            Vinf, &     !            "  free stream vel.
            Vy_ac, &    !            "  vertical vel. of aircraft
            Vx_ac, &    !            "  horizontal vel. of aircraft
            y_ac, &     !            "  altitude
            x_ac, &     !            "  ground position
            Lift, &     !            "  lift on a/c
            Drag, &     !            "  drag on a/c
            Vy_turb, &  !            "  vertical ghust vel.
            time,    &  !         time
            t_max       !   final time
End module

Program ac_in_turb
use ac_data
 use flight_data
 implicit none

 real*8 ::  dt, acc_y, y_i, y_avg_1, y_avg_2, pi = 4.d0 * atan( 1.d0 )
 real*8 ::  y_a, y_a_1, y_a_2, x_a, x_a_1, x_a_2, x_a_3, x_avg, i
 real*8 ::  acc_y_i_1, acc_y_i_2, acc_y_i_3, k1, k_2, k3 ,k4 ,k5 ,k6 ,k7 ,k8, k9, k10, k11, k12
 real*8 ::  y_i_1, y_i_2, y_i_3, y_i_4, y_i_5, y_i_6, Vinf_i_1, Vinf_i_2, Vinf_i_3, Lift_i_1, Lift_i_2, Lift_i_3

print*,"Choose a method:"
print*,"1)Euler's Method"
print*,"2)Heun's Method"
print*,"3)RK4 Method"
read*,i
print*,"Enter the step size:"
read*,dt
if(i .eq. 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!            !!!!!!!!!!!!!!!!
        !Euler's Method!            !Euler's Method!
        !!!!!!!!!!!!!!!!            !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call init_flight()

 open( 20, file = 'results.dat', form = 'formatted' )

 t_max = 25.d0
 acc_y = 0.d0

 write( 20, *)'" time [s]" "Vy_turb [m/2]" "Vinf [m/s]" "alpha [deg]" "a_y [m/s2]" "V_y [m/s]" "y [m]"'
 write( 20, '(7(e15.7,2x))') time, Vy_turb, Vinf, alpha * 180.d0/pi, acc_y, Vy_ac, y_ac

 do while ( time < t_max )
    ! turb vertical velocity
    Vy_turb = Vamp_turb * sin(kturb * (x_ac + Vx_ac * time))
    ! AOA and lift coef.
    alpha = atan((Vy_turb - Vy_ac)/(Vx_ac))
    CL = CL0 + alpha * dCL_dalpha
    ! air density
    ! ignore changes
    
    ! freestream vel.
    Vinf = sqrt( Vx_ac**2 + ( Vy_turb - Vy_ac )**2 )
    ! lift on a/c
    Lift = 0.5 * rho_inf * wing_area * (Vinf**2) * CL
    ! vertical acceleration
    acc_y = (Lift - (mass * grav)) / mass
    ! new time
    time = time + dt
    ! new vertical velocity
    Vy_ac = Vy_ac + dt * acc_y
    ! new vertical position
    y_ac = y_ac + dt * Vy_ac
    ! new horizonal ground vel. -- maintained / no change
    Vx_ac = Vground
    ! new a/c position
    x_ac = x_ac +  Vx_ac * dt
    write( 20, '(7(e15.7,2x))') time, Vy_turb, Vinf, alpha * 180.d0/pi, acc_y, Vy_ac, y_ac
 enddo

 close( 20 )
 stop
else if (i .eq. 2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!             !!!!!!!!!!!!!!!
        !Heun's Method!             !Heun's Method!
        !!!!!!!!!!!!!!!             !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call init_flight()
 open( 20, file = 'results.dat', form = 'formatted' )
 t_max = 25.d0
 acc_y = 0.d0
 write( 20, *)'" time [s]" "Vy_turb [m/2]" "Vinf [m/s]" "alpha [deg]" "a_y [m/s2]" "V_y [m/s]" "y [m]"'
 write( 20, '(7(e15.7,2x))') time, Vy_turb, Vinf, alpha * 180.d0/pi, acc_y, Vy_ac, y_ac
 
 do while ( time < t_max )
    ! turb vertical velocity
    Vy_turb = Vamp_turb * sin(kturb * (x_ac + Vx_ac * time))
    ! AOA and lift coef.
    alpha = atan((Vy_turb - Vy_ac)/(Vx_ac))
    CL = CL0 + alpha * dCL_dalpha
    ! air density
    ! ignore changes

    ! freestream vel.
    Vinf = sqrt( Vx_ac**2 + ( Vy_turb - Vy_ac )**2 )
    ! lift on a/c
    Lift = 0.5 * rho_inf * wing_area * (Vinf**2) * CL
    ! vertical acceleration
    acc_y = (Lift - (mass * grav)) / mass
    ! new time
    time = time + dt
    ! new vertical velocity
    y_i = acc_y
    y_i_1 = Vy_ac + dt * y_i        !predictor
    Vinf_i_1 = sqrt( Vx_ac**2 + ( Vy_turb - y_i_1 )**2 )            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Lift_i_1 = 0.5 * rho_inf * wing_area * (Vinf_i_1**2) * CL       !! NOTE: Vinf_i_1, Vinf_i_2, Vinf_i_3, Lift_i_1, Lift_i_2 and   !!
    acc_y_i_1 = (Lift_i_1 - (mass * grav)) / mass                   !! Lift_i_3 variables are just inter-variables to get the       !!
    y_i_2 = acc_y_i_1                                               !! variable acc_y_i_1                                           !!
    y_avg_1 = 0.5 * (y_i + y_i_2)                                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Vy_ac = Vy_ac + dt * y_avg_1         !corrector
    ! new vertical position
    y_a = Vy_ac
    y_a_1 = y_ac + dt * y_a       !predictor
    y_a_2 = Vy_ac
    y_avg_2 = (Vy_ac + Vy_ac)/2
    y_ac = y_ac + dt * y_avg_2       !corrector
    ! new horizonal ground vel. -- maintained / no change
    Vx_ac = Vground
    ! new a/c position
    x_a = Vx_ac
    x_a_1 = x_ac + dt * x_a        !predictor
    x_a_2 = Vx_ac
    x_avg = (Vx_ac + Vx_ac)/2
    x_ac = x_ac +  x_avg * dt        !corrector
    write( 20, '(7(e15.7,2x))') time, Vy_turb, Vinf, alpha * 180.d0/pi, acc_y, Vy_ac, y_ac
 enddo

 close( 20 )
 stop
 else if (i .eq. 3) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!             !!!!!!!!!!!!
        !RK4 Method!             !RK4 Method!
        !!!!!!!!!!!!             !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call init_flight()

 open( 20, file = 'results.dat', form = 'formatted' )
 t_max = 25.d0
 acc_y = 0.d0
 write( 20, *)'" time [s]" "Vy_turb [m/2]" "Vinf [m/s]" "alpha [deg]" "a_y [m/s2]" "V_y [m/s]" "y [m]"'
 write( 20, '(7(e15.7,2x))') time, Vy_turb, Vinf, alpha * 180.d0/pi, acc_y, Vy_ac, y_ac
 
 do while ( time < t_max )
    ! turb vertical velocity
    Vy_turb = Vamp_turb * sin(kturb * (x_ac + Vx_ac * time))
    ! AOA and lift coef.
    alpha = atan((Vy_turb - Vy_ac)/(Vx_ac))
    CL = CL0 + alpha * dCL_dalpha
    ! air density
    ! ignore changes

    ! freestream vel.
    Vinf = sqrt( Vx_ac**2 + ( Vy_turb - Vy_ac )**2 )
    ! lift on a/c
    Lift = 0.5 * rho_inf * wing_area * (Vinf**2) * CL
    ! vertical acceleration
    acc_y = (Lift - (mass * grav)) / mass
    ! new time
    time = time + dt
    ! new vertical velocity
    k1 = acc_y
    y_i_1 = Vy_ac + k1 * (0.5) * dt
    Vinf_i_1 = sqrt( Vx_ac**2 + ( Vy_turb - y_i_1 )**2 )              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Lift_i_1 = 0.5 * rho_inf * wing_area * (Vinf_i_1**2) * CL         !! NOTE: Vinf_i_1, Vinf_i_2, Vinf_i_3, Lift_i_1, Lift_i_2 and   !!
    acc_y_i_1 = (Lift_i_1 - (mass * grav)) / mass                     !! Lift_i_3 variables are just inter-variables to get the       !!
    k_2 = acc_y_i_1   !k2 was used as drag coefficient                !! variable acc_y_i_1                                           !!
                                                                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    y_i_2 = Vy_ac + k_2 * (0.5) * dt
    Vinf_i_2 = sqrt( Vx_ac**2 + ( Vy_turb - y_i_2 )**2 )
    Lift_i_2 = 0.5 * rho_inf * wing_area * (Vinf_i_2**2) * CL
    acc_y_i_2 = (Lift_i_2 - (mass * grav)) / mass
    k3 = acc_y_i_2

    y_i_3 = Vy_ac + k3 * dt
    Vinf_i_3 = sqrt( Vx_ac**2 + ( Vy_turb - y_i_3 )**2 )
    Lift_i_3 = 0.5 * rho_inf * wing_area * (Vinf_i_3**2) * CL
    acc_y_i_3 = (Lift_i_3 - (mass * grav)) / mass
    k4 = acc_y_i_3

    Vy_ac = Vy_ac + ((k1 + 2*k_2 + 2*k3 + k4) / 6) * dt   ! new vertical velocity

    ! new vertical position
    k5 = Vy_ac
    y_i_4 = y_ac + k5 * (0.5) * dt
    k6 = Vy_ac
    y_i_5 = y_ac + k6 * (0.5) * dt
    k7 = Vy_ac
    y_i_6 = y_ac + k7 * (0.5) * dt
    k8 = Vy_ac

    y_ac = y_ac + ((k5 + 2*k6 + 2*k7+ k8) / 6) * dt   ! new vertical position

    ! new horizonal ground vel. -- maintained / no change
    Vx_ac = Vground
    ! new a/c position
    k9 = Vx_ac
    x_a_1 = x_ac + k9 * (0.5) * dt
    k10 = Vx_ac
    x_a_2 = x_ac + k10 * (0.5) * dt
    k11 = Vx_ac
    x_a_3 = x_ac + k11 * (0.5) * dt
    k12 = Vx_ac
    
    x_ac = x_ac + ((k9 + 2*k10 + 2*k11+ k12) / 6) * dt

    write( 20, '(7(e15.7,2x))') time, Vy_turb, Vinf, alpha * 180.d0/pi, acc_y, Vy_ac, y_ac
 enddo

 close( 20 )
  endif

 stop
End

Subroutine init_flight()
 use ac_data
 use flight_data
 implicit none
 real*8 :: qinf, air_density

 time = 0.d0
 y_ac = h_alt0
 rho_inf = air_density( y_ac )
 Vinf =  Vground
 Vy_ac = 0.d0
 Vx_ac = Vground
 qinf =  0.5d0 * rho_inf * Vinf**2
 CL = ( mass * grav ) / ( qinf * wing_area )
 alpha = ( CL - CL0 ) / dCL_dalpha
 CD = CD0 + K2 * CL**2
 x_ac = 0.d0
 return
End subroutine

Function air_density( h )
 implicit none
 real*8 :: air_density, h
 air_density = 1.15d0
 return
End
