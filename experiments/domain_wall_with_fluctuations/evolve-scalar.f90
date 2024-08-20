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
  integer :: nTime = 128, spec = nyq
  integer :: alph = 8
  integer :: sim, minSim = 0, nSims = 500
  real(dl) :: phi0 = 1._dl

! bare potential stationary solution
!  real(dl) :: a1 = 6._dl**0.5, a2 = 0.06905335_dl, a3 = -6._dl**0.5 !a2=0.13810679_dl for nLat=512

!~~~ TEMPERATURE 1
!  real(dl) :: a1 = 2.351890666415152_dl, a2 = 0.03999959950776118_dl, a3 = -2.3559422418640295_dl
!  real(dl) :: a1 = 2.3523170384796726_dl, a2 = 0.05667848829121438_dl, a3 = -2.352220561433891_dl
!  real(dl) :: a1 = 2.3522385932116916_dl, a2 = 0.05668420140629302_dl, a3 = -2.352029254592668_dl
!  real(dl) :: a1 = 2.352237791464387_dl, a2 = 0.05668445281875601_dl, a3 = -2.3520163774572587_dl
  real(dl) :: a1 = 2.3522377973632045_dl, a2 = 0.05668446180814009_dl, a3 = -2.352015076098187_dl
  real(dl) :: temp = 1._dl

!~~~ TEMPERATURE 06
!  real(dl) :: a1 = 2.419804913744879_dl, a2 = 0.063461234624571_dl, a3 = -2.4197823766081_dl
!  real(dl) :: a1 = 2.419776240423927_dl, a2 = 0.06346139223997574_dl, a3 = -2.419750486038284_dl
!  real(dl) :: a1 = 2.419774031256123_dl, a2 = 0.06346161298071637_dl, a3 = -2.41974779037758_dl
!  real(dl) :: temp = 0.6_dl
  
  integer, parameter :: inFile = 70, cpFile = 71
  real(dl), dimension(:,:), pointer :: fld

  fld(1:nLat,1:2) => yvec(1:nVar-1) ! store the field in yvec
  time => yvec(nVar) ! last position stores the time?
  call setup(nVar)

  call initialize_rand(93286123,12)  
  do sim = 0, nSims-1 ! run nSims simulations for the parameters, each with its output files
      call initialise_fields(fld, len, nyq, phi0, spec, temp, m2eff, a1, a2, a3)
      if (sim >= minSim) then
          call time_evolve(len, sim, temp, lambda, m2eff, gam)
          print*, "Simulation ", sim, " in ", nSims-1 , "with T =", temp, " done!"
      endif
  enddo
  print*, "All Done!"

contains

  subroutine initialise_fields(fld, lenLat, knyq, phi, kcut, tem, m2, a1, a2, a3)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: a1, a2, a3
    integer, intent(in) :: knyq, kcut
    real(dl), intent(in) :: lenLat, phi, tem, m2
    integer :: ii

    do ii = 1, nLat
        fld(ii,1) = a1*(TANH(a2 * (ii-nLat/4._dl)) - TANH(a2*(ii-3._dl*nLat/4._dl))) + a3
    enddo
    fld(:,2) = 0._dl
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here

    call initialize_thermal_fluctuations(fld, lenLat, m2, knyq, phi, kcut, tem)
  end subroutine initialise_fields

  subroutine initialize_thermal_fluctuations(fld, lenLat, m2, knyq, phi, kcut, tem)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: lenLat, m2, phi, tem
    integer, intent(in), optional :: knyq, kcut

    real(dl), dimension(1:size(fld(:,1)/2+1)) :: spec, w2eff
    real(dl), dimension(1:size(fld(:,1))) :: df

    integer :: i
    real(dl) :: norm

    norm = 1._dl / phi / sqrt(2._dl * lenLat)

    do i = 1, nyq
       w2eff(i) = m2 + dk**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:knyq) = norm / w2eff(2:knyq)**0.25 * sqrt( 2._dl/(exp(w2eff(2:knyq)**0.5/tem)-1._dl))
    call generate_1dGRF(df,spec(:kcut))
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(:kcut))
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_thermal_fluctuations



  subroutine time_evolve(lenLat, sim, tem, lambda, m2, gam)
    real(dl) :: dt, dtout, lenLat, tem, lambda, m2, gam
    integer :: i, j, sim, outsize
    
    dt = dx/alph
    if (dt > dx) print*, "Warning, violating Courant condition" !i.e. alph > 1
    outsize = alph !* nLat/nTime
    dtout = dt*outsize
    
    do i = 1, nTime ! this loops over time slices
       call output_fields(fld, dt, dtout, lenLat, sim, tem, m2)
       do j = 1, outsize ! how many integrations are needed to evolve by one time slice; depends on how many crosses
          call gl10(yvec, dt, lambda, m2, gam)
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


  subroutine output_fields(fld, dt, dtout, lenLat, sim, tem, m2)
    real(dl), dimension(1:nLat, 1:2) :: fld
    real(dl) :: dt, dtout, lenLat, tem, m2
    logical :: o
    integer :: m, sim
    integer, parameter :: oFile = 98
    inquire(file='/gpfs/dpirvu/thermal_bubbles/wall_fluctuations_x'//trim(str(nLat))//'_temp'//trim(real_str(tem))//'_sim'//trim(str(sim))//'_fields.dat', opened=o)
    if (.not.o) then
       open(unit=oFile,file='/gpfs/dpirvu/thermal_bubbles/wall_fluctuations_x'//trim(str(nLat))//'_temp'//trim(real_str(tem))//'_sim'//trim(str(sim))//'_fields.dat')
       write(oFile,*) "# Lattice Parameters dx = ", dx, 'len = ', lenLat
       write(oFile,*) "# Time Stepping parameters dt = ", dt, "dt_out = ", dtout
       write(oFile,*) "# Other Parameters m2eff = ", m2, "dk = ", dk
    endif

    do m = 1, nLat
       write(oFile,*) fld(m,:)
    enddo
  end subroutine output_fields



end program Gross_Pitaevskii_1d

