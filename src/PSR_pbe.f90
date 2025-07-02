!**********************************************************************************************

subroutine psr_pbe()

!**********************************************************************************************
!
! Perfectly stirred reactor for the homogeneous PBE
! Stelios Rigopoulos (06/05/2017)
! Modified 25/06/2020
!
!**********************************************************************************************

implicit none

double precision, allocatable :: ni(:)

double precision moment(0:1)
double precision int_time,tin,meansize,dt,time

integer i,i_step,n_steps,iflag,flowflag,nin,i_write,n_write,i_writesp
integer agg_kernel_update,n_pbe_grid

! TEST variables
! double precision Temperature, Pw, RH, rho, time

!**********************************************************************************************

! Initialisation

! Initialise PBE
call pbe_read(n_pbe_grid)
call contrail_read()

allocate(ni(n_pbe_grid))
call pbe_grid()
call pbe_init(ni)
! Read PSR input data
open(30,file='psr/psr.in')
do i=1,2
  read(30,*)
end do
read(30,*) int_time
read(30,*) dt
read(30,*) agg_kernel_update
read(30,*) i_writesp
read(30,*) n_write
close(30)

! Read contrail input data


! Initialise PSR integration
n_steps = int_time/dt
i_write = 0
time = 0
!----------------------------------------------------------------------------------------------

open(22, file='pbe/plume_variables.out')

! Integration

! Write initial moments
call PBE_moments(ni,moment,meansize)

! initialise tau_m
call set_mixing_timescale()

call con_output(ni)

do i_step = 1,n_steps
  
  ! update contrail plume variables
  call update_plume_variables(time)
  call plume_var_output(22, time)
  ! The following should be done if the kernel should be updated at each time step due to e.g. 
  ! temperature dependency
  if (agg_kernel_update==1) then
    ! Insert here the expression for updating the kernel
    ! agg_kernel_const = 
    call PBE_agg_beta(2)
  end if

  ! Integrate
    call pbe_integ(ni,dt)

  ! Calculate moments
  call pbe_moments(ni,moment,meansize)

  ! Write moments
  ! Write PSD
  if ((i_write==n_write).or.(i_step==n_steps)) then
    i_write = 0
    call pbe_output(ni,i_writesp)
  end if
  i_write = i_write + 1

  time = time+dt
end do

!----------------------------------------------------------------------------------------------
close(22)

! Deallocate arrays
deallocate(ni)
call PBE_deallocate()

end subroutine psr_pbe

!**********************************************************************************************