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

   integer :: id_tr_u, id_tr_v
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
  character(len=8) :: var_name
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

  CS%ntr = 2 ! u, v velocity components
  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr))

  do m = 1, CS%ntr
    write(var_name(1:8), '(a6,i2.2)') 'pdemap', m
    tr_desc = var_desc(var_name, "1", "PDE Map", caller=mdl)
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

  CS%id_tr_u = register_diag_field("ocean_model", "pdemap_u", diag%axesTL, &
       day, "PDE map for u", "m s-1")
  CS%id_tr_v = register_diag_field("ocean_model", "pdemap_v", diag%axesTL, &
       day, "PDE map for v", "m s-1")
end subroutine initialize_pde_tracer

function filter_impulse(degree, cutoff, t)
  integer, intent(in) :: degree, t
  real, intent(in) :: cutoff
  real :: filter_impulse, an, bn, cn, dn
  integer :: n

  real, parameter :: pi = 4.0 * atan(1.0)

  filter_impulse = 0.

  do n = 1, degree / 2
    cn = cutoff * sin(pi / (2 * degree) * (2*n - 1))
    dn = cutoff * cos(pi / (2 * degree) * (2*n - 1))
    an = cn / degree
    bn = dn / degree

    filter_impulse = filter_impulse + exp(-cn * t) * (an * cos(dn * t) + bn * sin(dn * t))
  end do
end function filter_impulse

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

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) return

  if (.not. (associated(CS%u_ptr) .and. associated(CS%v_ptr))) call MOM_error(FATAL, "velocity pointers must be associated")

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m = 1, CS%ntr
      do k = 1,nz ; do j = js,je ; do i = is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo

      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
           evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    end do
  else
    do m = 1, CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    end do
  end if

  print *, "pde_tracer_column_physics, dt:", dt, ", position:", CS%position

  ! XXX check that dt divides CS%window evenly
  CS%position = CS%position + dt
  if (CS%position == CS%window) then
    CS%position = 0

    ! reset and post data
    if (CS%id_tr_u > 0) call post_data(CS%id_tr_u, CS%tr(:,:,:,1), CS%diag)
    if (CS%id_tr_v > 0) call post_data(CS%id_tr_v, CS%tr(:,:,:,2), CS%diag)

    CS%tr(:,:,:,:) = 0.
  end if

  do k = 1,nz ; do j = js,je ; do i = is,ie
    ! XXX does this only happen at the midpoint of the window?
    CS%tr(i,j,k,1) = CS%tr(i,j,k,1) - filter_impulse(CS%filter_degree, CS%filter_cutoff, CS%position) * (CS%u_ptr(I-1,j,k) + CS%u_ptr(I,j,k)) / 2
    CS%tr(i,j,k,2) = CS%tr(i,j,k,2) - filter_impulse(CS%filter_degree, CS%filter_cutoff, CS%position) * (CS%v_ptr(i,J-1,k) + CS%v_ptr(i,J,k)) / 2
  enddo ; enddo ; enddo
end subroutine pde_tracer_column_physics

subroutine pde_tracer_end(CS)
  type(pde_tracer_CS), pointer :: CS

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  end if
end subroutine pde_tracer_end

end module pde_tracer
