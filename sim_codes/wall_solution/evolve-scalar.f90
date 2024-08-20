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
  integer :: nTime = nLat
  integer :: step = 1, alph = 8
  real(dl) :: a1 = 6._dl**0.5, a2 = 0.06905335_dl, a3 = -6._dl**0.5 !a2=0.13810679_dl for nLat=512
  integer, parameter :: inFile = 70, cpFile = 71
  real(dl), dimension(:,:), pointer :: fld

  fld(1:nLat,1:2) => yvec(1:nVar-1) ! store the field in yvec
  time => yvec(nVar) ! last position stores the time?
  call setup(nVar)

  call initialize_rand(93286123,12)
  call initialise_fields(fld, a1, a2, a3)
  call time_evolve(step, lambda, m2eff, gam)
  print*, 'm2eff = ', m2eff, ' len = ', len, 'dt = ', dx/alph, 'dx = ', dx, 'dk = ', dk, 'alph = ', alph

contains

  subroutine initialise_fields(fld, a1, a2, a3)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: a1, a2, a3
    integer :: ii

    do ii = 1, nLat
        !fld(ii,1) = a1 * TANH(a2 * ii) + a3
        fld(ii,1) = a1*(TANH(a2 * (ii-nLat/4._dl)) - TANH(a2*(ii-3._dl*nLat/4._dl))) + a3
    enddo

    fld(:,2) = 0._dl
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here

  end subroutine initialise_fields


  subroutine time_evolve(step, lambda, m2eff, gam)
    real(dl) :: dt, dtout, lambda, m2eff, gam
    integer :: i, j, step, outsize
    
    dt = dx/alph
    if (dt > dx) print*, "Warning, violating Courant condition" !i.e. alph > 1
    outsize = alph !* nLat/nTime
    dtout = dt*outsize
    
    do i = 1, nTime ! this loops over time slices
       do j = 1, outsize ! how many integrations are needed to evolve by one time slice; depends on how many crosses
          call gl10(yvec, dt, lambda, m2eff, gam)
       enddo

       if (i > nTime-step*nLat) then
           if (mod(i, step) == 0) then
               call output_fields(fld, dt, dtout)
           endif
       endif
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


  subroutine output_fields(fld, dt, dtout)
    real(dl), dimension(1:nLat, 1:2) :: fld
    real(dl) :: dt, dtout
    logical :: o
    integer :: m
    integer, parameter :: oFile = 98
    inquire(file='/gpfs/dpirvu/thermal_bubbles/phi4_gam'//trim(real_str(gam))//'_t'//trim(str(nTime))//'_x'//trim(str(nLat))//'_fields.dat', opened=o)
    if (.not.o) then
       open(unit=oFile,file='/gpfs/dpirvu/thermal_bubbles/phi4_gam'//trim(real_str(gam))//'_t'//trim(str(nTime))//'_x'//trim(str(nLat))//'_fields.dat')
       write(oFile,*) "# Lattice Parameters dx = ", dx, 'len = ', len
       write(oFile,*) "# Time Stepping parameters dt = ", dt, "dt_out = ", dtout
       write(oFile,*) "# Other Parameters m2eff = ", m2eff, "dk = ", dk
    endif

    do m = 1, nLat
       write(oFile,*) fld(m,1)
    enddo
  end subroutine output_fields

end program Gross_Pitaevskii_1d

