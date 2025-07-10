!**********************************************************************************************
!
! PBE finite volume discretisation for growth
!
!**********************************************************************************************



!**********************************************************************************************

subroutine growth_tvd(ni,index,growth_source, w_or_i, Jw_l, Jw_r, n_satw)

!**********************************************************************************************
!
! Growth for finite volume method
! Size can be diameter, volume etc.
!
! Stelios Rigopoulos, Fabian Sewerin, Binxuan Sun
!
!**********************************************************************************************

use pbe_mod

use con_mod, Only: T, amb_p, pi, sat_pressure_w_func, sw

implicit none

double precision, dimension(m), intent(in) :: ni
integer, intent(in)                        :: index
double precision, intent(out)              :: growth_source

double precision :: phi
double precision :: g_terml
double precision :: g_termr

double precision :: gnl,gnr           !< (G*N) at left surface and right surface
double precision :: nl                !< Number density in left cell
double precision :: nll               !< Number density in left-left cell
double precision :: nc                !< Number density in this cell
double precision :: nr                !< Number density in right cell
double precision :: nrr               !< Number density in right-right cell
double precision :: eps               !< Tolerance for upwind ratio
double precision :: rl,rr             !< r+ at left and right surface

parameter(eps = 1.D1*epsilon(1.D0))

!**********************************************************************************************
!Contrail variables

integer, intent(in) :: w_or_i         ! Parameter determines wether water (1) or ice (0) growth is used

double precision, intent(out) :: Jw_r
double precision, intent(out) :: Jw_l
double precision :: volH20 = 0.03d0 ! in nm**3 
double precision :: alpha_w = 1.d0
double precision :: radius_water_r, radius_water_l
double precision :: thermal_speed
double precision, intent(out) :: n_satw
double precision :: diffusion_coef
double precision :: kB = 1.380649d0 * 10.d0**(-23)
double precision :: massH20 = 2.99d0 * 10.d0**(-26)

!**********************************************************************************************

radius_water_r = (v(index) * (3.d0 / (4.d0 * pi))) ** (1.d0 / 3.d0 ) * 10.d0**(-9) !m
radius_water_l = (v(index-1) * (3.d0 / (4.d0 * pi))) ** (1.d0 / 3.d0 ) * 10.d0**(-9)

thermal_speed = sqrt( (8 * kB * T) / (massH20)) ! m/s might need to change to nm/s
n_satw = sat_pressure_w_func(T) / (kB * T)
diffusion_coef = (2.66d0 * 10.d0 ** (-5)) * (T / 298.15d0)**(1.81d0) * (1.01325 * 10**5) / (amb_p) ! units m^2/s

Jw_r = (pi * radius_water_r**2 * alpha_w * thermal_speed * n_satw * sw) / (1 + (alpha_w * thermal_speed * radius_water_r) / (4 * diffusion_coef))
Jw_l = (pi * radius_water_l**2 * alpha_w * thermal_speed * n_satw * sw) / (1 + (alpha_w * thermal_speed * radius_water_l) / (4 * diffusion_coef))
!Only growth to the right present at nucleation interval

if (growth_function==1) then

  !Size-independent growth
  g_termr = g_coeff1
  g_terml = g_coeff1

else if (growth_function==2) then

  !Power-law growth
  g_termr = g_coeff1*(v(index)**g_coeff2)
  g_terml = g_coeff1*(v(index-1)**g_coeff2)

else if (growth_function==3.AND.w_or_i==1) then

  g_termr = Jw_r * volH20 ! units [nm**3/s]
  g_terml = Jw_l * volH20

end if

!----------------------------------------------------------------------------------------------
!TVD scheme ref: S.Qamar et al 2006: A comparative study of high resolution schemes for solving
!                population balances in crystallization
!----------------------------------------------------------------------------------------------

if ((sw>0.d0).AND.(g_termr >=0).AND.(g_terml>=0)) then

  ! growth rate is along positive direction
  if (index==1) then

    gnl = 0.0D0
    gnr = g_termr * 0.5 * (ni(1)+ni(2))

  else if (index==m) then

    rl = (ni(m) - ni(m-1) + eps) / (ni(m-1) - ni(m-2) + eps)
    phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
    gnl = g_terml * (ni(m-1) + 0.5 * phi * (ni(m-1) - ni(m-2)))
    gnr = g_termr * (ni(m) + 0.5*(ni(m) - ni(m-1)))

  else

    ! Fluxes at cell right surface
    nl = ni(index-1)
    nc = ni(index)
    nr = ni(index+1)
    rr = (nr - nc + eps) / (nc - nl + eps)
    phi = max(0.0d0, min(2.0d0 * rr, min((1.0d0 + 2.0d0 * rr) / 3.0d0, 2.0d0)))
    gnr = g_termr * (nc + 0.5 * phi * (nc - nl))

    ! Fluxes at cell left surface
    if (index==2 ) then
      gnl = g_terml * 0.5 * (ni(1)+ni(2))
    else
      nl = ni(index-2)
      nc = ni(index-1)
      nr = ni(index)
      rl = (nr - nc + eps) / (nc - nl + eps)
      phi = max(0.0d0, min(2.0d0 * rl, min((1.0d0 + 2.0d0 * rl) / 3.0d0, 2.0d0)))
      gnl = g_terml * (nc + 0.5 * phi * (nc - nl))
    end if
  end if
end if

if (i_gm==1) then
  ! For mass-conservative growth scheme, apply it after the first interval
  if (index>1) then
    growth_source = (v(index-1)*gnl - v(index)*gnr) / (0.5*(v(index)**2-v(index-1)**2))
  else
    growth_source = (gnl - gnr) / dv(index)
  end if
else
  growth_source = (gnl - gnr) / dv(index)
end if

end subroutine growth_tvd

!**********************************************************************************************