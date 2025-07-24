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

double precision, allocatable :: nsoot(:) 
double precision, allocatable :: nwater(:) 
double precision, allocatable :: nice(:) 

double precision moment(0:1)
double precision int_time,tin,meansize,dt,time

integer i,i_step,n_steps,iflag,flowflag,nin,i_write,n_write,i_writesp,nout_dt
integer agg_kernel_update,n_pbe_grid
double precision :: output_timestep
! TEST variables
! double precision Temperature, Pw, RH, rho, time

!**********************************************************************************************

! Initialisation
! Initialise PBE
call pbe_read(n_pbe_grid)
call contrail_read()

allocate(nsoot(n_pbe_grid))
allocate(nwater(n_pbe_grid))
allocate(nice(n_pbe_grid))

nwater = 0.0
nice= 0.0
call pbe_grid()
call pbe_init(nsoot)

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
read(30,*) output_timestep
close(30)


! Initialise PSR integration
n_steps = int(int_time/dt)
nout_dt = int(output_timestep/dt)

i_write = 0
time = 0


!----------------------------------------------------------------------------------------------

open(22, file='pbe/plume_variables.csv', status='replace', action='write')
write(22,*) 'Temperature (K),smw,sw,si,Ambient Density (kg/m^3),time'
open(33, file='pbe/distribution_data.csv')
write(33,*) 'radius (nm),nsoot,nwater,nice,cumulative n,sw,si,time'
open(55, file='pbe/moments.csv')
write(55,*) 'time(s),n_soot(#/m^3),n_water(#/m^3),n_ice(#/m^3),r_i(nm),r_act(nm)'
! Integration

! Write initial moments
call PBE_moments(nsoot,moment,meansize)

! initialise tau_m
call set_mixing_timescale()

do i_step = 1,n_steps
  
  ! update contrail plume variables
  call update_plume_variables(nsoot, nwater, nice, time, dt)
  ! The following should be done if the kernel should be updated at each time step due to e.g. 
  ! temperature dependency
  
  if (agg_kernel_update==1) then
    ! Insert here the expression for updating the kernel
    ! agg_kernel_const = 
    call PBE_agg_beta(2)
  end if

  !Droplet activation and freezing

  call droplet_freeze(nwater, nice)
  call droplet_activation(nsoot, nwater)  

  ! Integrate
  call pbe_integ(nwater,dt,1)
  call pbe_integ(nice,dt,0)

  ! Calculate moments
  call pbe_moments(nsoot,moment,meansize)

  call con_output(nsoot, nwater, nice, 33, time,i_step,nout_dt)
  call plume_var_output(22, time,i_step, nout_dt)
  call moments_output(nsoot, nwater, nice, 55, time, i_step, nout_dt)

  ! Write moments
  ! Write PSD
  if ((i_write==n_write).or.(i_step==n_steps)) then
    i_write = 0
    call pbe_output(nsoot,i_writesp)
  end if
  i_write = i_write + 1

  time = time+dt
end do

!----------------------------------------------------------------------------------------------
close(22)
close(33)
close(55)

! Deallocate arrays
deallocate(nsoot)
deallocate(nwater)
deallocate(nice)

call PBE_deallocate()
write(*,*) 'Check'
end subroutine psr_pbe

!**********************************************************************************************