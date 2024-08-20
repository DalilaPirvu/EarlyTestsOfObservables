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
  integer :: nTimeMax = 30000
  integer :: sim, nSims = 50, minSim = 0
  integer :: spec = nyq
  real(dl), parameter :: alph = 8._dl
  real(dl), parameter :: phi0 = twopi / 7.5_dl
  integer, parameter :: inFile = 70, cpFile = 71
  real(dl), dimension(:,:), pointer :: fld

  fld(1:nLat,1:2) => yvec(1:nVar-1) ! store the field in yvec
  time => yvec(nVar) ! last position stores the time?
  call initialize_rand(93286123,12)
  call setup(nVar)

  do sim = 0, nSims-1 ! run nSims simulations for the parameters, each with its output files
      call initialise_fields(fld)
      if (sim >= minSim) then !(any(sim == goodSims)) then ! (sim >= minSim) then
          call time_evolve(sim)
          print*, "Simulation ", sim, " in ", nSims , " done!"
      endif    
  end do

contains

  subroutine initialise_fields(fld)
    real(dl), dimension(:,:), intent(inout) :: fld

    fld(:,1) = 0.5_dl*twopi
    fld(:,2) = 0._dl
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here

    call initialize_vacuum_fluctuations(fld) ! Change this call as necessary
  end subroutine initialise_fields

!!!!!!!!!!!!!!!!!!
! Time Evolution !
!!!!!!!!!!!!!!!!!!

  subroutine time_evolve(sim) 
    real(dl) :: dt, dtout, sum_cos_fld
    integer :: i, j, k, m, sim
    logical :: bool

    dt = dx/alph
    if (dt > dx) print*, "Warning, violating Courant condition"
    dtout = dt*alph
    call output_fields(fld, dt, dtout, sim)
    
    i = 1
    k = 1
    bool = .False.

    do while ( i <= 2000 )
!       print*,  k

       do j = 1, int(dtout/dt)
          call gl10(yvec,dt)
       end do
       
!       if (k > nTimeMax-1024) then
       call output_fields(fld, dt, dtout, sim)
!       endif
       if ( .not. bool ) then
          sum_cos_fld = 0._dl
          do m = 1, nLat
             sum_cos_fld = sum_cos_fld + cos(fld(m,1))
          end do
          sum_cos_fld = sum_cos_fld/nLat
!          print*,  sum_cos_fld

          if ( k == nTimeMax ) then
             exit
          end if

          if ( sum_cos_fld >= -0.5_dl ) then
             bool = .True.
             print*, k
          end if

          k = k + 1
       else
          i = i + 1
!          print*,  i
       end if
    end do
  end subroutine time_evolve

!!!!!!!!!!!!!!!!
! Fluctuations !
!!!!!!!!!!!!!!!!

  subroutine initialize_vacuum_fluctuations(fld)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), dimension(1:size(fld(:,1)/2+1)) :: spectrum, w2eff
    real(dl), dimension(1:size(fld(:,1))) :: df
    real(dl) :: norm; integer :: i

    norm = (0.5_dl)**0.5 / phi0 / sqrt(2._dl) / sqrt(lenLat)
    do i = 1, nyq
       w2eff(i) = m2eff + dk**2*(i-1)**2
    end do
    spectrum = 0._dl
    spectrum(2:nyq) = norm / w2eff(2:nyq)**0.25
    call generate_1dGRF(df,spectrum(:spec))
    fld(:,1) = fld(:,1) + df(:)

    spectrum = spectrum * w2eff**0.5
    call generate_1dGRF(df,spectrum(:spec))
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_vacuum_fluctuations

!!!!!!!!!!
! Output !
!!!!!!!!!! 

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

  subroutine output_fields(fld, dt, dtout, sim)
    real(dl), dimension(1:nLat, 1:2) :: fld
    real(dl), dimension(1:nLat) :: gsq
    real(dl) :: dt, dtout!, energy
    logical :: o; integer :: m, sim
    integer, parameter :: oFile = 98

    tPair%realSpace(:) = fld(:, 1)
    call gradsquared_1d_wtype(tPair, dk)
    gsq(:) = tPair%realSpace(:)   
   
    inquire(file='/gpfs/dpirvu/sims/x'//trim(str(nLat))//'_phi0'//trim(real_str(phi0))//'_lambda'//trim(real_str(lambda))//'_sim'//trim(str(sim))//'_fields.dat', opened=o)
    if (.not.o) then
       open(unit=oFile,file='/gpfs/dpirvu/sims/x'//trim(str(nLat))//'_phi0'//trim(real_str(phi0))//'_lambda'//trim(real_str(lambda))//'_sim'//trim(str(sim))//'_fields.dat')
       write(oFile,*) "# Lattice Parameters dx = ", dx, "nLat = ", nLat, " lenLat = ", lenLat, "spec = ", spec, "dk = ", dk
       write(oFile,*) "# Time Stepping parameters dt = ", dt, "dt_out = ", dtout
       write(oFile,*) "# Other Parameters m2eff = ", m2eff, "phi0 = ", phi0
    endif

!    energy = 0._dl
    do m = 1, nLat
       write(oFile,*) fld(m,1)!, 0.5_dl*gsq(m), 4._dl*nu*(-cos(fld(m,1)) + 0.5_dl*lambda**2*sin(fld(m,1))**2)
!       write(oFile,*) fld(m,:), 0.5_dl*fld(m,:)**2, 0.5_dl*gsq(m), 4._dl*nu*(-cos(fld(m,1)) + 0.5_dl*lambda**2*sin(fld(m,1))**2)
!       energy = energy + 0.5_dl*fld(m,2)**2._dl + 0.5_dl*gsq(m) + 4._dl*nu*( - cos(fld(m,1)) + 0.5_dl*lambda**2._dl * sin(fld(m,1))**2._dl )
    end do
!    print*, energy
  end subroutine output_fields

end program Gross_Pitaevskii_1d
