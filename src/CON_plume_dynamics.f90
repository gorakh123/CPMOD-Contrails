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

! print *, 'taum = ', tau_m
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
    
end subroutine update_plume_variables