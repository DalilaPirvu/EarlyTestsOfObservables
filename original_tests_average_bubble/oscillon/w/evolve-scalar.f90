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
  integer :: nTimeMax = 20000
  integer :: sim, minSim = 0, nSims = 5
  integer :: spec = nyq, param
!  integer :: ns = 5
  real(dl), parameter :: alph = 8._dl

  real(dl) :: phi0! = twopi / 7._dl
  integer, parameter :: inFile = 70, cpFile = 71
  real(dl), dimension(:,:), pointer :: fld

  fld(1:nLat,1:2) => yvec(1:nVar-1) ! store the field in yvec
  time => yvec(nVar) ! last position stores the time?
  call initialize_rand(93286123,12)
  call setup(nVar)

  do param = 11, 12, 2 
     phi0 = param / 10._dl
     do sim = 0, nSims-1 ! run nSims simulations for the parameters, each with its output files
        if (sim >= minSim) then
           call initialise_fields(fld)
!           call time_evolve2222222(sim, ns)
           call time_evolve(sim)
           print*, "Simulation ", sim, " in ", nSims , " done!"
        end if
     end do
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

  subroutine time_evolve2222222(sim, ns) 
    real(dl) :: dt, dtout, mean, mean0, stdev0
    real(dl), dimension(1:nLat) :: cosslice
    integer :: i, j, k, sim, ns
    logical :: bool

    dt = dx/alph
    if (dt > dx) print*, "Warning, violating Courant condition"
    dtout = dt*alph

    j = 1
    bool = .False.
    do while ( .not. bool )

       do i = 1, int(dtout/dt)
          call gl10(yvec,dt)
       end do
       call output_fields(fld, dt, dtout, sim)

       do k = 1, nLat
          cosslice(k) = cos(fld(k,1))
       end do

       if (j == 1) then
          mean0 = sum( cosslice ) / size( cosslice )
          stdev0 = sqrt( sum( cosslice**2 ) / size( cosslice ) - mean0**2 )
       else
          mean = sum( cosslice ) / size( cosslice )       

          if ( mean >= mean0 + ns * stdev0 ) then
             bool = .True.
             print*, j
          end if
       end if
       
       if ( j == nTimeMax ) then
          exit
       end if

       j = j + 1
    end do
  end subroutine time_evolve2222222
  
  subroutine time_evolve(sim) 
    real(dl) :: dt, dtout, sum_cos_fld
    integer :: i, j, k, m, sim
    logical :: bool

    dt = dx/alph
    if (dt > dx) print*, "Warning, violating Courant condition"
    dtout = dt*alph

    i = 1
    k = 1
    bool = .False.

    do while ( i <= 1 )
!       print*,  k

       do j = 1, int(dtout/dt)
          call gl10(yvec,dt)
       end do
       call output_fields(fld, dt, dtout, sim)

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

          if ( sum_cos_fld > -0.7_dl ) then
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
    real(dl) :: dt, dtout
    logical :: o; integer :: m, sim
    integer, parameter :: oFile = 98

    inquire(file='/gpfs/dpirvu/sims/typeP'//trim(str(typeP))//'_len'//trim(real_str(lenLat))//'_phi0'//trim(real_str(phi0))//'_lamb'//trim(real_str(lambda))//'_x'//trim(str(nLat))//'_sim'//trim(str(sim))//'_fields.dat', opened=o)
    if (.not.o) then
    open(unit=oFile,file='/gpfs/dpirvu/sims/typeP'//trim(str(typeP))//'_len'//trim(real_str(lenLat))//'_phi0'//trim(real_str(phi0))//'_lamb'//trim(real_str(lambda))//'_x'//trim(str(nLat))//'_sim'//trim(str(sim))//'_fields.dat')
       write(oFile,*) "# Lattice Parameters dx = ", dx, "nLat = ", nLat, " lenLat = ", lenLat, "spec = ", spec, "dk = ", dk
       write(oFile,*) "# Time Stepping parameters dt = ", dt, "dt_out = ", dtout
       write(oFile,*) "# Other Parameters m2eff = ", m2eff, "phi0 = ", phi0
    endif

    do m = 1, nLat
       write(oFile,*) fld(m,1)
    end do
  end subroutine output_fields

end program Gross_Pitaevskii_1d
