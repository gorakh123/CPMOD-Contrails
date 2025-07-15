!**********************************************************************************************
!
! PBE solver subroutines
! Stelios Rigopoulos
!
!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_integ(ni,dt)

!**********************************************************************************************
!
! Temporal integration
!
! Stelios Rigopoulos (02/06/2019)
!
! Modified 14/06/06
! Modified 19/12/2017
! Modified 13/05/2018
! Modified 02/06/2019
! Modified 25/06/2020
!
!**********************************************************************************************

use pbe_mod

implicit none

integer index 

double precision, dimension(m), intent(inout) :: ni
double precision, intent(in)                  :: dt

double precision niprime(m),nitemp(m)
double precision k1(m),k2(m),k3(m),k4(m)
!----------------------------------------------------------------------------------------------

if (solver_pbe == 1) then

  !Euler explicit
  call pbe_ydot(ni,niprime,dt)
  ni = ni + niprime * dt

else if (solver_pbe == 2) then

  !Runge-Kutta 2nd order
  call pbe_ydot(ni,niprime, dt)
  nitemp = ni + 0.5D0 * niprime * dt
  call pbe_ydot(nitemp,niprime,dt)
  ni = ni + niprime * dt

else if (solver_pbe == 3) then

  !Runge-Kutta 4th order
  call pbe_ydot(ni,niprime,dt)
  k1 = niprime * dt
  nitemp = ni + 0.5D0 * k1
  call pbe_ydot(nitemp,niprime,dt)
  k2 = niprime * dt
  nitemp = ni + 0.5D0 * k2
  call pbe_ydot(nitemp,niprime,dt)
  k3 = niprime * dt
  nitemp = ni + k3
  call pbe_ydot(nitemp,niprime,dt)
  k4 = niprime * dt
  ni = ni + (1.D0 / 6.D0) * k1 + (1.D0 / 3.D0) * k2 + (1.D0 / 3.D0) * k3 + (1.D0 / 6.D0) * k4

end if

! ensure n_i does not fall below 0

do index = 1,m 
  if (ni(index) < 0.d0) then

    if (index==m) then
      ni(index) = 0.d0
    else if ((ni(index + 1) - ni(index)) > 0 ) then
      ni(index + 1) = ni(index + 1) + ni(index) ! attempt to ensure conservation of number density
      ni(index) = 0.d0
    else if (((ni(index - 1) - ni(index) )> 0).AND.(m>1)) then
      ni(index-1) = ni(index-1) + ni(index)
    else
      ni(index) = 0
    end if

  end if
end do

end subroutine pbe_integ

!**********************************************************************************************



!**********************************************************************************************

subroutine pbe_ydot(ni,niprime,timestep)

!**********************************************************************************************
!
! Calculates the right hand side of the ODEs to be integrated
!
! By Stelios Rigopoulos
! 14/01/2002
! Modified 04/05/2017
! Modified 23/06/2020
! Modified 25/06/2020
! Modified 05/07/2020
!
!**********************************************************************************************

use pbe_mod
use con_mod, only: sw, smw, P_w, pi, pvap, pvap_m, dpvapdt_m, p_wsat, T

implicit none

double precision, dimension(m), intent(in)  :: ni
double precision, dimension(m), intent(out) :: niprime

double precision dn(m)
double precision, intent(in) :: timestep

double precision growth_source,growth_mass_source,params(1)
double precision :: Lw_index, dpvapsink_index
double precision :: gterml, gtermr
double precision :: Lw_total, dpvapsink_total 
double precision :: nwsat
double precision :: volH20 = 0.03d0 !nm**3

double precision :: gnl, gnr, nil, nir

integer index

double precision :: kB = 1.380649d0 * 10.d0**(-23)

!----------------------------------------------------------------------------------------------

niprime = 0.
params(1) = 0.

!Nucleation
if (max_nuc>0) then
! Note: expressions for i_gm==0 and i_gm==1 are equivalent but written separately for clarity
  if (i_gm==0) then
    do index = 1,max_nuc
      niprime(index) = nuc(index)/dv(index)
    end do
  else if (i_gm==1) then
    do index = 1,max_nuc
      niprime(index) = nuc(index)*v_m(index)/(0.5*(v(index)**2-v(index-1)**2))
    end do
  end if
end if

!Growth
if ((sw>0.d0)) then
  Lw_total = 0
  dpvapsink_total = 0
  !write(*,*) 'T ', 'sw ', 'thermal speed ', 'n_w,sat ', 'v(index) ', 'g_termr '
  do index = 1,m
    call growth_tvd(ni,index,growth_source,1,gterml,gtermr,nwsat, gnl, gnr)
    niprime(index) = niprime(index) + growth_source
    nwsat = nwsat * 10.d-9 ! ensure units are consistant
    Lw_index = ((4 * pi) / (volH20 * nwsat)) * (ni(index) *dv(index)) * (v_m(index)**2) * ((gterml + gtermr)*0.5D0)
    dpvapsink_index = ((4 * pi * p_wsat) / (volH20 * nwsat)) * (ni(index) *dv(index)) * (v_m(index)**2) * ((gterml + gtermr)*0.5D0) !tweak nwsat
    !write(*,*) Lw_index

    if ((Lw_index == Lw_index).AND.(Lw_index>=0.d0)) then
      Lw_total = Lw_total + Lw_index
      dpvapsink_total = dpvapsink_total + dpvapsink_index
    end if

    nil = (gnl / dv(index)) * timestep
    nir = (gnr / dv(index)) * timestep

    ! Attmept to track dry core radius

    ! if (index==1) then
    !   rcore_array(index) = rcore_array(index)
    ! else
    !   if (((ni(index) + nil - nir) > 0).AND.(ni(index)>nir)) then
    !     rcore_array(index) =(((ni(index)-nir) * rcore_array(index)) + (nil * rcore_array(index-1))) / (ni(index) + nil - nir)
    !   else if (((ni(index) + nil -nir) > 0).AND.(ni(index)<=nir)) then
    !     rcore_array(index) = rcore_array(index-1)
    !   else
    !     rcore_array(index) = 0


    ! ensure ni(index) does not turn negative

    !   if (ni(index) > nir) then
    !     rcore_array(index) =(((ni(index)-nir) * rcore_array(index) ) + (ni * rcore_array(index-1))) / (ni(index) + nil - nir)
    !   else if (nir >= ni(index) + nil) then
    !     rcore_array(index) = 0
    !   else if (nir >= ni(index)) then
    !     rcore_array(index) = rcore_array(index-1)
    !   end if
    !  end if
 ! end if

  end do
  !write(*,*) 'Lw total: ', Lw_total, 'Pw :', P_w , 'nwsat: ', nwsat, 'volH20: ', volH20, 'last gterml: ', gterml
  !sw = sw + (-Lw_total * timestep)
  pvap = pvap + (-dpvapsink_total * timestep)
  !sw = pvap / p_wsat + 1
end if

!Aggregation
if (agg_kernel>0) then
  ! CFV formulation of Liu and Rigopoulos (2019)
  ! Note 1: current value of niprime is augmented within pbe_agg_cfv
  ! Note 2: contracting grid is not implemented
  call pbe_agg_cfv(dv,v_m,ni,niprime)
end if

!Fragmentation
if (break_const>0.) then
  call pbe_breakage_cfv(niprime,ni)
end if

end subroutine pbe_ydot

!**********************************************************************************************