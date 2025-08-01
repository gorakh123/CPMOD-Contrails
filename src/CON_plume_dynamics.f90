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

module activation_funcs

    real(kind=8), save :: tmp_r_dry, tmp_T, tmp_kappa
    real(kind=8), parameter :: Mw = 18.016_8/1.0e3_8 !! molar weight of water [kg/mol]
    real(kind=8), parameter :: rho_w = 1000.0_8 !! density of water [kg/m^3]

    contains
    function fminbound(f, a, b, tol) result(x_min)
        !! Find the minimizer of the function f(x) on the interval [a,b].
        !! Uses the golden-section search algorithm.
        !!
        !! Arguments:
        !!   f   : external real function of one variable
        !!   a,b : real, the lower and upper bounds
        !!   tol : real, the convergence tolerance
        !!
        !! Returns:
        !!   x_min : real, approximate minimizer of f on [a,b].
        !! Obtained from https://github.com/lboscagli/Population-Balance-of-Particles-in-Flows/blob/main/src_mp/module_thermo.f90
        implicit none
        interface
        function f(x) result(fval)
            real(kind=8), intent(in) :: x
            real(kind=8) :: fval
        end function f
        end interface
        real(kind=8), intent(in) :: a, b, tol
        real(kind=8) :: x_min
        real(kind=8) :: gr, c, d, fc, fd, left, right

        gr = (sqrt(5.0_8)-1.0_8)/2.0_8
        left = a; right = b
        c = right - (right-left)*gr
        d = left + (right-left)*gr
        fc = f(c); fd = f(d)

        do while (abs(right-left) > tol)
        if (fc < fd) then
            right = d
            d = c
            fd = fc
            c = right - (right-left)*gr
            fc = f(c)
        else
            left = c
            c = d
            fc = fd
            d = left + (right-left)*gr
            fd = f(d)
        end if
        end do
        x_min = (left+right)/2.0_8
    end function fminbound

    function sigma_w(T) result(sigma)
        !! Compute surface tension of water at temperature T [K].
        !! Returns sigma [J/m^2].
        real(kind=8), intent(in) :: T
        real(kind=8) :: sigma
        sigma = 0.0761_8-1.55e-4_8*(T-273.15_8)
    end function sigma_w

    function Seq(r, r_dry, T, kappa) result(S_eq)
        !! Compute equilibrium supersaturation over a particle per κ-Köhler theory.
        !! r,r_dry [m]; T [K]; kappa [-].
        !! Returns S_eq [-].
        real(kind=8), intent(in) :: r, r_dry, T, kappa
        real(kind=8) :: S_eq, A, B
        A = (2.0_8*Mw*sigma_w(T))/(8.314_8*T*rho_w*r)
        B = (r**3-r_dry**3)/(r**3-(r_dry**3)*(1.0_8-kappa))
        S_eq = exp(A)*B - 1.0_8
    end function Seq

    function neg_Seq_fixed(r) result(fval)
        real(kind=8), intent(in) :: r
        real(kind=8) :: fval
        fval = -1.0_8*Seq(r, tmp_r_dry, tmp_T, tmp_kappa)
    end function neg_Seq_fixed


    subroutine kohler_crit(T, r_dry, kappa, r_crit, s_crit)
    !! Compute Köhler critical radius and supersaturation for dry particle.
    !! T [K], r_dry [m], kappa [-].
    !! approx: logical - use approximate eqn.
    !! r_crit, s_crit are returned.  
        implicit none
        real(kind=8), intent(in) :: T, r_dry, kappa
        real(kind=8), intent(out) :: r_crit, s_crit
            tmp_r_dry = r_dry
        tmp_T = T
        tmp_kappa = kappa
        ! Call minimizer on neg_Seq_fixed
        r_crit = fminbound(neg_Seq_fixed, r_dry, r_dry*1.0e4_8, 1.0e-10_8)
        s_crit = Seq(r_crit, r_dry, T, kappa)  
    end subroutine kohler_crit

    subroutine critical_curve(T, kappa, rs, rcrits, scrits)
    !! Generate Köhler critical radii and supersaturations across dry sizes.
    !! T [K], r_a,r_b [m] dry radius bounds; kappa [-]; approx logical.
    !! rs(n), rcrits(n), scrits(n) returned.
    
        use pbe_mod, only: m
        real(kind=8), intent(in) :: T, kappa, rs(m)
        real(kind=8), intent(out) :: rcrits(m), scrits(m)
        integer :: i
        real(kind=8) :: r_crit, s_crit
        do i = 1,m
            call kohler_crit(T, rs(i), kappa, r_crit, s_crit)
            rcrits(i) = r_crit
            scrits(i) = s_crit
        end do
    end subroutine critical_curve

end module activation_funcs



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

subroutine update_plume_variables(nsoot, namb, namb_mixed, nvola, nwater, nice, time, dt)

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
double precision, dimension(m), intent(inout) :: nvola
double precision, dimension(m), intent(inout) :: namb_mixed
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
    T_prev = T
    si_prev = si
    dpvapdt_m = 0.d0
else
    dpvapdt_m = (pvap_m - pvap_mprev) / dt
end if 

nsoot = (nsoot * amb_rho * dilution) / (amb_rhoprev * Dilution_prev)
namb_mixed = namb * (amb_temp / T) * (1 - dilution)
nvola = (nvola * amb_rho * dilution) / (amb_rhoprev * Dilution_prev)
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

subroutine droplet_activation(nsoot, namb, nvola, nwater)

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
double precision, dimension(m), intent(inout) :: nvola

double precision :: activation_radius_soot, activation_radius_amb, activation_radius_vola, kelvin_radius
double precision :: kappa_soot, kappa_amb, kappa_vola

integer :: i

kelvin_radius = 1 ! [nm]
kappa_soot = 0.005d0
kappa_amb = 0.5d0
kappa_vola = 0.5d0
activation_radius_soot = (kelvin_radius / (54.d0 * kappa_soot)**(1.d0/3.d0)) * (sw ** (-2.d0/3.d0))
activation_radius_amb = (kelvin_radius / (54.d0 * kappa_amb)**(1.d0/3.d0)) * (sw ** (-2.d0/3.d0))
activation_radius_vola = (kelvin_radius / (54.d0 * kappa_vola)**(1.d0/3.d0)) * (sw ** (-2.d0/3.d0))

do i = 1, m

    if ((v_m(i) >= activation_radius_amb).AND.(sw >= 0)) then
        nwater(i) = nwater(i) + namb(i)
        namb(i) = 0.d0
    end if


    if ((v_m(i) >= activation_radius_soot).AND.(sw >= 0)) then
        nwater(i) = nwater(i) + nsoot(i)
        nsoot(i) = 0.d0
    end if

    if ((v_m(i) >= activation_radius_vola).AND.(sw >= 0)) then
        nwater(i) = nwater(i) + nvola(i)
        nvola(i) = 0.d0
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


subroutine droplet_activation2(nsoot, namb, nvola, nwater)

    use pbe_mod, only: m, v_m   
    use con_mod
    use activation_funcs

    implicit none

    double precision, dimension(m), intent(inout) :: nsoot
    double precision, dimension(m), intent(inout) :: namb
    double precision, dimension(m), intent(inout) :: nwater
    double precision, dimension(m), intent(inout) :: nvola
    double precision, dimension(m) :: radius_m

    double precision, dimension(m) :: scrit_soot, scrit_amb, scrit_vola
    double precision, dimension(m) :: rcrit_soot, rcrit_amb, rcrit_vola

    double precision :: kappa_amb = 0.5
    double precision :: kappa_soot = 0.005
    double precision :: kappa_vola = 0.5 ! Duble check

    integer :: index

    radius_m = v_m * 10.d-9

    call critical_curve(T, kappa_soot, radius_m, rcrit_soot, scrit_soot)
    call critical_curve(T, kappa_amb, radius_m, rcrit_amb, scrit_amb)
    call critical_curve(T, kappa_vola, radius_m, rcrit_vola, scrit_vola)

    do index = 1, m

        if ((scrit_soot(index) <= sw).AND.(sw >= 0)) then
            nwater(index) = nwater(index) + nsoot(index)
            nsoot(index) = 0
        end if

        if ((scrit_amb(index) <= sw).AND.(sw >= 0)) then
            nwater(index) = nwater(index) + namb(index)
            namb(index) = 0
        end if

        if ((scrit_vola(index) <= sw).AND.(sw >= 0)) then
            nwater(index) = nwater(index) + nvola(index)
            nvola(index) = 0
        end if

    end do

end subroutine droplet_activation2