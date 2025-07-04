!**********************************************************************************************
!
! Implementation of the plume dynamics/dilution
!
! Implemented the dilution model used by Karcher et al 2015
! Equations are solved numerically, using forward euler 
! Note: There is a parameter Beta which a senstivity study can be conducted on
!
! By Gorakh Adhikari 27/06/2025
!
!**********************************************************************************************

subroutine set_mixing_timescale()

!**********************************************************************************************
!
! subroutine determines the mixing time scale, tau_m
!
! By Gorakh Adhikari 26/06/2025
!
!**********************************************************************************************

use con_mod 

implicit none

epsilon = 0.0285

x_m = r_0 * sqrt(2 / epsilon)
tau_m = x_m / initial_velocity

write(*,*) 'tau_m (s):', tau_m
end subroutine set_mixing_timescale

!**********************************************************************************************

subroutine update_plume_variables(time)

!**********************************************************************************************
! 
! subroutine updates plume_variables:
!   Temperature, inital plume area,
!
! It follows the dilution model described in:
!   Karcher et al 2015
!
! Note that the Beta parameter contrails the intesity 
! of turbulent mixing (can conduct sensitivity study on it)
! 
! By Gorakh Adhikari 26/06/2025
! 
!**********************************************************************************************

use con_mod

implicit none

double precision dilution 
double precision beta
double precision, intent(in) :: time

beta = 0.9 

if (time > tau_m) then

    dilution = (tau_m/time) ** beta
    T = amb_temp + (initial_temp - amb_temp) * dilution

else

    T = initial_temp 
    dilution = 1
    
end if

p_wsat = sat_pressure_w_func(T)
p_isat = sat_pressure_i_func(T)

amb_pw = amb_Si * sat_pressure_i_func(amb_temp)

smw = ((amb_pw + G * (T - amb_temp)) / p_wsat) - 1

sw = smw 

end subroutine update_plume_variables

!**********************************************************************************************

subroutine droplet_activation(nsoot, nwater)

!**********************************************************************************************
!
! Subroutine determines the soot particles that can activate
! Following from the activation equation outlined in Karcher et al 2015
!
! TODO: implement a newton-rhapson method for the activation
!
! By Gorakh Adhikari 03/07/25
!
!**********************************************************************************************

use pbe_mod, only: m, v_m     
use con_mod

implicit none

double precision, dimension(m), intent(inout) :: nsoot
double precision, dimension(m), intent(inout) :: nwater 

double precision :: activation_radius, kelvin_radius, activation_vol
double precision :: kappa

integer :: i

kelvin_radius = 1 ! [nm]
kappa = 0.005
activation_radius = (kelvin_radius / (54.d0 * kappa)**(1.d0/3.d0)) * (sw ** (-2.d0/3.d0))
activation_vol = ((4.d0 * pi) / 3 ) * activation_radius**3

do i = 1, m

    if ((v_m(i) >= activation_vol).AND.(sw >= 0)) then
        nwater(i) = nwater(i) + nsoot(i)
        nsoot(i) = 0.d0 
    end if

end do 

end subroutine droplet_activation