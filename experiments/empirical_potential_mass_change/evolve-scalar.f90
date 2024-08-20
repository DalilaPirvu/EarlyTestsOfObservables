!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS

#include "macros.h"
#define SMOOTH 1

program Gross_Pitaevskii_1d
  ! the following modules are incorporated:
  use, intrinsic :: iso_c_binding
  use gaussianRandomField
  use integrator
  use constants
  use eom
  implicit none

  real(dl), pointer :: time
  integer :: nTime = 10*nLat, spec = nyq
  integer :: step = 1, alph = 8
  integer :: sim, minSim = 0, nSims = 50
  real(dl) :: temp = 1._dl
  real(dl) :: phi0 = 1._dl
  real(dl) :: lambda = 1._dl !0.423_dl !1._dl
  real(dl) :: m2bare = 1._dl
  real(dl) :: m2eff = 0.4_dl

  integer, parameter :: inFile = 70, cpFile = 71
  real(dl), dimension(:,:), pointer :: fld

  fld(1:nLat,1:2) => yvec(1:nVar-1) ! store the field in yvec
  time => yvec(nVar) ! last position stores the time?
  call setup(nVar)


  call initialize_rand(93286123,12)  

  do sim = 0, nSims ! run nSims simulations for the parameters, each with its output files
      call initialise_fields(fld, nyq, phi0, spec, temp, m2bare)

      if (sim >= minSim) then
          call time_evolve(sim, step, temp, lambda, m2bare)
          print*, "Simulation ", sim, " in ", nSims , "with", temp, " done!"
      endif

  enddo
  print*, 'm2 in V = ', m2bare, 'm2 in field = ', m2eff, 'len = ', len, 'dt = ', dx/alph, 'dx = ', dx, 'dk = ', dk, 'phi0 = ', phi0, 'alph = ', alph, 'spec = ', spec

contains

  subroutine initialise_fields(fld, kmax, phi, klat, temp, m2inphi)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: kmax, klat
    real(dl), intent(in) :: phi, temp, m2inphi

    fld(:,1) = 0._dl ! free field
    fld(:,2) = 0._dl
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here

    call initialize_thermal_fluctuations(fld, len, m2inphi, kmax, phi, klat, temp)
  end subroutine initialise_fields


  subroutine initialize_thermal_fluctuations(fld, len, m2inphi, km, phi, kc, temp)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: len, m2inphi, phi, temp
    integer, intent(in), optional :: km, kc

    real(dl), dimension(1:size(fld(:,1)/2+1)) :: spec, w2eff
    real(dl), dimension(1:size(fld(:,1))) :: df

    integer :: i
    real(dl) :: norm

    norm = 1._dl / phi / sqrt(2._dl * len)

    do i = 1, nyq
       w2eff(i) = m2inphi + dk**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:km) = norm / w2eff(2:km)**0.25 * sqrt( 2._dl/(exp(w2eff(2:km)**0.5/temp)-1._dl))
    call generate_1dGRF(df,spec(:kc))
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(:kc))
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_thermal_fluctuations



  subroutine time_evolve(sim, step, temp, lamb, m2inV)
    real(dl) :: dt, dtout, temp, lamb
    real(dl), intent(in) :: m2inV
    integer :: i, j, sim, step, outsize
    
    dt = dx/alph
    if (dt > dx) print*, "Warning, violating Courant condition" !i.e. alph > 1
    outsize = alph !* nLat/nTime
    dtout = dt*outsize
    
    do i = 1, nTime ! this loops over time slices
       call output_fields(fld, dt, dtout, sim, temp, lamb, m2inV)
       do j = 1, outsize ! how many integrations are needed to evolve by one time slice; depends on how many crosses
          call gl10(yvec, dt, lamb, m2inV)
       enddo
   enddo
  end subroutine time_evolve
  

  subroutine setup(nVar)
    integer, intent(in) :: nVar
    call init_integrator(nVar)
    call initialize_transform_1d(tPair,nLat)
  end subroutine setup

  character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str

  character(len=20) function real_str(k)
    real(dl), intent(in) :: k
    write (real_str, '(f12.4)') k
    real_str = adjustl(real_str)
  end function real_str


  subroutine output_fields(fld, dt, dtout, sim, temp, lamb, m2inV)
    real(dl), dimension(1:nLat, 1:2) :: fld
    real(dl) :: dt, dtout, temp, lamb, m2inV
    logical :: o
    integer :: m, sim
    integer, parameter :: oFile = 98
    
    inquire(file='/gpfs/dpirvu/thermal_bubbles/phi4_empirical_potential_null_m2'//trim(real_str(m2inV))//'_t'//trim(str(nTime))//'_x'//trim(str(nLat))//'_temp'//trim(real_str(temp))//'_phi0'//trim(real_str(phi0))//'_lamb'//trim(real_str(lamb))//'_sim'//trim(str(sim))//'_fields.dat', opened=o)
    
    if (.not.o) then
       open(unit=oFile,file='/gpfs/dpirvu/thermal_bubbles/phi4_empirical_potential_null_m2'//trim(real_str(m2inV))//'_t'//trim(str(nTime))//'_x'//trim(str(nLat))//'_temp'//trim(real_str(temp))//'_phi0'//trim(real_str(phi0))//'_lamb'//trim(real_str(lamb))//'_sim'//trim(str(sim))//'_fields.dat')
    
       write(oFile,*) "# Lattice Parameters dx = ", dx, 'len = ', len
       write(oFile,*) "# Time Stepping parameters dt = ", dt, "dt_out = ", dtout
       write(oFile,*) "# Other Parameters 'm2 in V = ", m2inV, "dk = ", dk
    endif

    do m = 1, nLat
       write(oFile,*) fld(m,:)
    enddo
  end subroutine output_fields



end program Gross_Pitaevskii_1d

