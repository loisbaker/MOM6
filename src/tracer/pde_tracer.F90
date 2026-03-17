module pde_tracer

use MOM_diag_mediator, only : diag_ctrl, post_data, register_diag_field
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : vardesc, var_desc
use MOM_restart, only : MOM_restart_CS
use MOM_time_manager, only : time_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_pde_tracer, initialize_pde_tracer
public pde_tracer_column_physics, pde_tracer_end
public register_pde_state_pointers

type, public :: pde_tracer_CS ; private
   integer :: window = 0 !< Width of advection window [s]
   integer :: position = 0 !< Position inside the advection window [s]

   integer :: filter_degree !< Degree (*2) of the Butterworth filter
   real :: filter_cutoff !< Cutoff frequency for the filter

   real, pointer :: u_ptr(:,:,:) => NULL()
   real, pointer :: v_ptr(:,:,:) => NULL()

   integer :: ntr = 0
   type(tracer_registry_type), pointer :: tr_Reg => NULL()
   real, pointer :: tr(:,:,:,:) => NULL()
   type(diag_ctrl), pointer :: diag => NULL()
   type(MOM_restart_CS), pointer :: restart_CSp => NULL()

   integer :: id_tr_u_filt, id_tr_v_filt, id_tr_u_map, id_tr_v_map
end type pde_tracer_CS

contains

function register_pde_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type), intent(in) :: HI
  type(verticalGrid_type), intent(in) :: GV
  type(param_file_type), intent(in) :: param_file
  type(pde_tracer_CS), pointer :: CS
  type(tracer_registry_type), pointer :: tr_Reg
  type(MOM_restart_CS), target, intent(inout) :: restart_CS

  logical :: register_pde_tracer
#include "version_variable.h"
  character(len=40) :: mdl = "pde_tracer"
  character(len=13) :: var_name
  type(vardesc) :: tr_desc
  real, pointer :: tr_ptr(:,:,:) => NULL()

  integer :: isd, ied, jsd, jed, nz
  integer :: m

  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_pde_tracer called with an "// &
         "associated control structure.")
  end if
  allocate(CS)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, MDL, "PDE_TRACER_WINDOW", CS%window, &
       "The width of the advection window for PDE tracers (in seconds). " // &
       "This determines the length of Lagrangian trajectories.", default=6*3600)

  call get_param(param_file, mdl, "PDE_TRACER_FILTER_DEGREE", CS%filter_degree, &
       "The order of the Butterworth filter to use as an impulse function.", &
       default=2)

  call get_param(param_file, mdl, "PDE_TRACER_FILTER_CUTOFF", CS%filter_cutoff, &
       "The cutoff angular frequency (in [rad s-1]) of the Butterworth filter.", &
       default=5e-5, units="rad s-1")

  CS%ntr = 4 ! u, v filtered velocities and maps
  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr))

  do m = 1, CS%ntr
    write(var_name(1:13), '(a6,i2.2)') 'tracer_filt', m
    tr_desc = var_desc(var_name, "1", "Filtered tracer", caller=mdl)
    tr_ptr => CS%tr(:,:,:,m)

    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, tr_desc=tr_desc, &
         registry_diags=.false., restart_CS=restart_CS, mandatory=.false.)
  end do

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS

  register_pde_tracer = .true.
end function register_pde_tracer

subroutine register_pde_state_pointers(CS, u, v)
  type(pde_tracer_CS), pointer :: CS
  real, dimension(:,:,:), target :: u, v

  if (.not. associated(CS)) return

  CS%u_ptr => u
  CS%v_ptr => v

end subroutine register_pde_state_pointers

subroutine initialize_pde_tracer(restart, day, G, GV, h, diag, CS)
  logical, intent(in) :: restart
  type(time_type), target, intent(in) :: day
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h
  type(diag_ctrl), target, intent(in) :: diag
  type(pde_tracer_CS), pointer :: CS

  if (.not. associated(CS)) return

  CS%tr(:,:,:,:) = 0.
  CS%diag => diag

  CS%id_tr_u_filt = register_diag_field("ocean_model", "u_filt", diag%axesTL, &
       day, "Filtered u", "m s-1")
  CS%id_tr_v_filt = register_diag_field("ocean_model", "v_filt", diag%axesTL, &
       day, "Filtered v", "m s-1")
  CS%id_tr_u_map = register_diag_field("ocean_model", "u_map", diag%axesTL, &
       day, "Map for u", "m")
  CS%id_tr_v_map = register_diag_field("ocean_model", "v_map", diag%axesTL, &
       day, "Map for v", "m")
end subroutine initialize_pde_tracer

function filter_impulse_response(degree, cutoff, window, t)
  integer, intent(in) :: degree, window, t
  real, intent(in) :: cutoff
  real :: filter_impulse_response, an, bn, cn, dn, norm_correction
  integer :: n

  real, parameter :: pi = 4.0 * atan(1.0)

  filter_impulse_response = 0.
  norm_correction = 0.

  do n = 1, degree / 2
    cn = cutoff * sin(pi / (2 * degree) * (2*n - 1))
    dn = cutoff * cos(pi / (2 * degree) * (2*n - 1))
    an = cn / degree
    bn = dn / degree

    filter_impulse_response = filter_impulse_response &
      + exp(-cn * abs(t)) * (an * cos(dn * abs(t)) + bn * sin(dn * abs(t)))

    norm_correction = norm_correction + 2 * exp( -cn * window / 2) * (( an * cn + bn * dn)/( cn**2 + dn**2) * cos(dn * window / 2) + ( bn * cn - an * dn )/( cn**2 + dn**2) * sin(dn * window / 2))
  end do
  filter_impulse_response = filter_impulse_response / (1.0 - norm_correction)
end function filter_impulse_response

function filter_integrated_impulse_response(degree, cutoff, window, t)
  ! Finds the integral from -window/2 to t of the impulse response
  integer, intent(in) :: degree, window, t
  real, intent(in) :: cutoff
  real :: filter_integrated_impulse_response, an, bn, cn, dn, norm_correction
  integer :: n

  real, parameter :: pi = 4.0 * atan(1.0)

  filter_integrated_impulse_response = 0.
  norm_correction = 0.

  do n = 1, degree / 2
    cn = cutoff * sin(pi / (2 * degree) * (2*n - 1))
    dn = cutoff * cos(pi / (2 * degree) * (2*n - 1))
    an = cn / degree
    bn = dn / degree

    factor1 = (an * cn + bn * dn)/( cn**2 + dn**2)
    factor2 = ( bn * cn - an * dn )/( cn**2 + dn**2)

    if (t <= 0) then
      filter_integrated_impulse_response = filter_integrated_impulse_response &
        + factor1 * (exp(cn * t) * cos(dn * t) - exp(-cn * window / 2) * cos(dn * window / 2)) &
        + factor2 * (- exp(cn * t) * sin(dn * t) - exp(-cn * window / 2) * sin(dn * window / 2))
    else
      filter_integrated_impulse_response = filter_integrated_impulse_response &
        - 2 * exp( -cn * t) * (factor1 * cos(dn * t) + factor2 * sin(dn * t))
    end if

    norm_correction = norm_correction + 2 * exp( -cn * window / 2) * (factor1 * cos(dn * window / 2) + factor2 * sin(dn * window / 2))
  end do

  if (t <= 0) then
    filter_integrated_impulse_response = filter_integrated_impulse_response / (1.0 - norm_correction)

  else
    filter_integrated_impulse_response = (1 - filter_integrated_impulse_response) / (1.0 - norm_correction)/2 + 0.5
  end if

end function filter_integrated_impulse_response

function heaviside(t)
  integer, intent(in) :: t
  real :: heaviside

  if (t >= 0) then
    heaviside = 1.0
  else
    heaviside = 0.0
  end if
end function heaviside

subroutine pde_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_old, h_new, ea, eb
  type(forcing), intent(in) :: fluxes
  real, intent(in) :: dt
  type(unit_scale_type), intent(in) :: US
  type(pde_tracer_CS), pointer :: CS
  real, optional, intent(in) :: evap_CFL_limit, minimum_forcing_depth

  integer :: i, j, k, is, ie, js, je, nz, m
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work
  real :: u_on_h, v_on_h, impulse_response, midpoint_mask, heaviside_factor, integrated_impulse_response

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) return

  if (.not. (associated(CS%u_ptr) .and. associated(CS%v_ptr))) call MOM_error(FATAL, "velocity pointers must be associated")

  print *, "pde_tracer_column_physics, dt:", dt, ", position:", CS%position

  ! XXX check that dt divides CS%window evenly
  CS%position = CS%position + dt
  !if (CS%position == CS%window) then
  if (abs(CS%position - CS%window) < 0.5 * dt) then
    CS%position = 0
    print *, "At window endpoint, resetting filter tracers"

    ! reset and post data (temporarily move this outside of the if statement so we can see fields at all timesteps)
     if (CS%id_tr_u_filt > 0) call post_data(CS%id_tr_u_filt, CS%tr(:,:,:,1), CS%diag)
     if (CS%id_tr_v_filt > 0) call post_data(CS%id_tr_v_filt, CS%tr(:,:,:,2), CS%diag)
     if (CS%id_tr_u_map > 0) call post_data(CS%id_tr_u_map, CS%tr(:,:,:,3), CS%diag)
     if (CS%id_tr_v_map > 0) call post_data(CS%id_tr_v_map, CS%tr(:,:,:,4), CS%diag)

    CS%tr(:,:,:,:) = 0.
  end if
  
  ! If we want to output these fields every timestep, we do that here
  ! if (CS%id_tr_u_filt > 0) call post_data(CS%id_tr_u_filt, CS%tr(:,:,:,1), CS%diag)
  ! if (CS%id_tr_v_filt > 0) call post_data(CS%id_tr_v_filt, CS%tr(:,:,:,2), CS%diag)
  ! if (CS%id_tr_u_map > 0) call post_data(CS%id_tr_u_map, CS%tr(:,:,:,3), CS%diag)
  ! if (CS%id_tr_v_map > 0) call post_data(CS%id_tr_v_map, CS%tr(:,:,:,4), CS%diag)
  
  ! Check if we are within the "midpoint" timestep
  if (abs(CS%position - (real(CS%window) / 2.0)) < (0.5 * dt)) then
    midpoint_mask = 1.0
    print *, "midpoint mask is 1, position:", CS%position
  else
    midpoint_mask = 0.0
  end if

  ! every timestep, add the impulse response
  do k = 1,nz ; do j = js,je ; do i = is,ie
    ! interpolations onto tracer point
    u_on_h = (CS%u_ptr(I-1,j,k) + CS%u_ptr(I,j,k)) / 2
    v_on_h = (CS%v_ptr(i,J-1,k) + CS%v_ptr(i,J,k)) / 2

    ! impulse response at this point in the window (might want to centre this on the timestep)
    impulse_response = filter_impulse_response(CS%filter_degree, CS%filter_cutoff, CS%window, CS%window / 2 - CS%position)
    heaviside_factor = heaviside(CS%position - CS%window / 2)
    !integrated_impulse_response = filter_integrated_impulse_response(CS%filter_degree, CS%filter_cutoff, CS%window, -CS%window / 2 + CS%position)

    ! Start with just the low pass (i.e. don't hit with the full velocity at the midpoint)
    CS%tr(i,j,k,1) = CS%tr(i,j,k,1) - dt * impulse_response * u_on_h + midpoint_mask * u_on_h ! This definition finds the wave component
    CS%tr(i,j,k,2) = CS%tr(i,j,k,2) - dt * impulse_response * v_on_h + midpoint_mask * v_on_h ! This definition finds the wave component
    ! CS%tr(i,j,k,1) = CS%tr(i,j,k,1) + dt * impulse_response * u_on_h ! This definition finds the mean component
    ! CS%tr(i,j,k,2) = CS%tr(i,j,k,2) + dt * impulse_response * v_on_h ! This definition finds the mean component
    CS%tr(i,j,k,3) = CS%tr(i,j,k,3) - dt * heaviside_factor * u_on_h ! These maps remap to midpoint position
    CS%tr(i,j,k,4) = CS%tr(i,j,k,4) - dt * heaviside_factor * v_on_h ! These maps remap to midpoint position
    ! CS%tr(i,j,k,3) = CS%tr(i,j,k,3) - dt * integrated_impulse_response * u_on_h ! These maps remap to mean position
    ! CS%tr(i,j,k,4) = CS%tr(i,j,k,4) - dt * integrated_impulse_response * v_on_h ! These maps remap to mean position
  enddo; enddo ; enddo
end subroutine pde_tracer_column_physics

subroutine pde_tracer_end(CS)
  type(pde_tracer_CS), pointer :: CS

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  end if
end subroutine pde_tracer_end

end module pde_tracer
