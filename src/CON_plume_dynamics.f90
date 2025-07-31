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

subroutine update_plume_variables(nsoot, namb, nwater, nice, time, dt)

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
use pbe_mod, only: m
use con_mod

implicit none

double precision dilution 
double precision beta
double precision, intent(in) :: time
double precision, intent(in) :: dt
double precision, dimension(m), intent(inout) :: nwater
double precision, dimension(m), intent(inout) :: nice
double precision, dimension(m), intent(inout) :: nsoot
double precision, dimension(m), intent(inout) :: namb
double precision :: R = 287.d0 !Specific gas constant of air J/kgK

amb_pvap = amb_Si * sat_pressure_i_func(amb_temp)

beta = 0.9d0 

if (time > tau_m) then

    dilution = (tau_m/time) ** beta
    T = amb_temp + (initial_temp - amb_temp) * dilution
    dTdt = -(beta / time) * (T - amb_temp)

else

    T = initial_temp 
    dilution = 1.d0
    dTdt = 0.d0
    
end if

p_wsat = sat_pressure_w_func(T)
p_isat = sat_pressure_i_func(T)

amb_pw = amb_Si * sat_pressure_i_func(amb_temp)

smw = ((amb_pw + G * (T - amb_temp)) / p_wsat) - 1.d0
pvap_m = amb_pw + G * (T - amb_temp)

amb_rho = amb_p / (T * R)

if (time == 0.d0) then
    sw = smw
    si = pvap_m / p_isat - 1
    pvap = pvap_m
    Dilution_prev = dilution
    amb_rhoprev = amb_rho
    si_prev = si
    dpvapdt_m = 0.d0
else
    dpvapdt_m = (pvap_m - pvap_mprev) / dt
end if 

nsoot = (nsoot * amb_rho * dilution) / (amb_rhoprev * Dilution_prev)
namb = (namb * amb_rho * dilution) / (amb_rhoprev * Dilution_prev)
nwater = (nwater * amb_rho * dilution) / (amb_rhoprev * Dilution_prev)
nice = (nice * amb_rho * dilution) / (amb_rhoprev * Dilution_prev)



smw_prev = smw
pvap_mprev = pvap_m
amb_rhoprev = amb_rho
Dilution_prev = dilution
pvap = pvap + (dpvapdt_m * dt)

sw = pvap / p_wsat - 1.d0
si = pvap / p_isat - 1.d0
dsidt = (si - si_prev) / dt

if ((dsidt < 0).AND.(amb_Si<1.d0).AND.(pvap < amb_pvap)) then
    pvap = amb_pvap

else if ((dsidt < 0).AND.(amb_Si>=1.d0).AND.(pvap < p_isat)) then
    pvap = p_isat
    
end if

sw = pvap / p_wsat - 1
si = pvap / p_isat - 1

si_prev = si


end subroutine update_plume_variables

!**********************************************************************************************

subroutine droplet_activation(nsoot, namb, nwater)

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
double precision, dimension(m), intent(inout) :: namb
double precision, dimension(m), intent(inout) :: nwater 

double precision :: activation_radius_soot, activation_radius_amb, kelvin_radius
double precision :: kappa_soot, kappa_amb

integer :: i

kelvin_radius = 1 ! [nm]
kappa_soot = 0.005
kappa_amb = 0.5
activation_radius_soot = (kelvin_radius / (54.d0 * kappa_soot)**(1.d0/3.d0)) * (sw ** (-2.d0/3.d0))
activation_radius_amb = (kelvin_radius / (54.d0 * kappa_amb)**(1.d0/3.d0)) * (sw ** (-2.d0/3.d0))

do i = 1, m

    if ((v_m(i) >= activation_radius_amb).AND.(sw >= 0)) then
        nwater(i) = nwater(i) + namb(i)
        namb(i) = 0.d0
    end if


    if ((v_m(i) >= activation_radius_soot).AND.(sw >= 0)) then
        nwater(i) = nwater(i) + nsoot(i)
        nsoot(i) = 0.d0
    end if

end do 

end subroutine droplet_activation

!**********************************************************************************************


subroutine droplet_freeze(nwater, nice)

!**********************************************************************************************
!
! Subroutine determines the droplets that freeze
! Following from the activation equation outlined in Karcher et al 2015
!
! 
!
! By Gorakh Adhikari 14/07/25
!
!**********************************************************************************************

    use pbe_mod, only: v_m, m
    use con_mod

    implicit none

    double precision, dimension(m), intent(inout) :: nwater
    double precision, dimension(m), intent(inout) :: nice
    double precision :: LWV
    double precision :: tau_frz
    double precision :: J_ice
    double precision :: lambda

    integer index
    tau_frz = 1.d0/( -3.574d0 * dTdt)
    J_ice = 10.d0**6 * (exp(-3.574d0*T + 858.719d0))
    do index = 1, m
        LWV = ((v_m(index)**3) * (4.d0*pi)/3.d0) * 10.d-27
        ! LWV = ((v_m(index)**3 - mean_radius**3) * (4.d0*pi)/3.d0) * 10.d-27
        lambda = LWV * J_ice * tau_frz
        !write(*,*) 'lambda: ', lambda, 'LWV: ', LWV, 'J_ice: ', J_ice, 'tau_frz: ', tau_frz,'Temp: ', T
        

        if ((lambda>=1.d0).AND.(nwater(index)>=0.d0)) then
            !write(*,*) 'Checkpoint'
            nice(index) = nice(index) + nwater(index)
            nwater(index) = 0.d0
        end if
    end do

end subroutine droplet_freeze