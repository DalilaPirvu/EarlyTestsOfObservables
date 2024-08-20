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
  integer :: sim, temploop, spec = nyq, nTime = nLat, minSim = 0, nSims = 50
  integer :: step = 1, alph = 8
  real(dl) :: temp, phi0 = 1._dl
  real(dl) :: m2eff, m2adj
 
  ! From m_eff at late time
  !real(dl), dimension(15) :: list_sigmasq = (/ 0.18300042184974438, 0.5672659480566075, 0.9940202797266966, 1.4374291642214636, 1.889074131084186, 2.3451148856978867, 2.803460811596109, 3.2628494674902173, 3.722467596210206, 4.181767980669381, 4.640370541351879, 5.0980047269897355, 5.554474128859033, 6.009633872211461, 6.463375715117995 /)
  
  ! From m_PS prediction
  ! real(dl), dimension(15) :: list_sigmasq = (/ 0.18083783870455825, 0.5627056797813478, 0.9886455321769851, 1.4330353130809275, 1.887558649241655, 2.3484286912167347, 2.8135811526525307, 3.2817646574261667, 3.7521683764508422, 4.224242563147683, 4.6976017051146375, 5.171968081570496, 5.647137056205483, 6.1229548736881085, 6.599304006792529 /)
  
  integer, parameter :: inFile = 70, cpFile = 71
  real(dl), dimension(:,:), pointer :: fld

  fld(1:nLat,1:2) => yvec(1:nVar-1) ! store the field in yvec
  time => yvec(nVar) ! last position stores the time?
  call setup(nVar)

  do temploop = 1, 15, 1
      call initialize_rand(93286123,12)
      temp = temploop * 1._dl
      do sim = 0, nSims-1 ! run nSims simulations for the parameters, each with its output files
          
          m2eff = 1._dl
          m2adj = m2eff!*(1._dl + 0.5_dl * lambda * list_sigmasq(temploop))

          call initialise_fields(fld, nyq, phi0, spec, temp, m2adj)
          
!          m2adj = m2eff*(1._dl + 0.5_dl * lambda * stddevsq(fld(:,1)))
          print*, temp, m2adj, m2eff!, list_sigmasq(temploop), stddevsq(fld(:,1))

          if (sim >= minSim) then
              call time_evolve(sim, step, temp, lambda, m2eff, gam)
              print*, "Simulation ", sim+1, " in ", nSims , "with", lambda, temp, " done!"
          endif
      enddo
  enddo
  print*, 'm2eff = ', m2eff, 'len = ', len, 'dt = ', dx/alph, 'dx = ', dx, 'dk = ', dk, 'phi0 = ', phi0, 'alph = ', alph, 'spec = ', spec

contains

  subroutine initialise_fields(fld,kmax,phi,klat,temp,m2)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: kmax, klat
    real(dl), intent(in) :: phi,temp,m2

    fld(:,1) = 0._dl ! free field
    fld(:,2) = 0._dl
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here

    call initialize_thermal_fluctuations(fld,len,m2,kmax,phi,klat,temp)
  end subroutine initialise_fields


  subroutine time_evolve(sim, step, temp, lambda, m2, gam)
    real(dl) :: dt, dtout, temp, lambda, gam, m2
    integer :: i, j, sim, step, outsize
    
    dt = dx/alph
    if (dt > dx) print*, "Warning, violating Courant condition" !i.e. alph > 1
    outsize = alph !* nLat/nTime
    dtout = dt*outsize

    do i = 1, nTime ! this loops over time slices
       call output_fields(fld, dt, dtout, sim, temp)
       do j = 1, outsize ! how many integrations are needed to evolve by one time slice; depends on how many crosses
          call gl10(yvec, dt, lambda, m2, gam)
       enddo
   enddo
  end subroutine time_evolve


  ! Calculates standard deviation for given set of values
  real(dl) function stddevsq(vals) result(sigma)
    real(dl), dimension(:), intent(in) :: vals
    real(dl) :: mean
    integer :: n

    n = size(vals)
    mean = sum(vals)/n
    sigma = sum((vals - mean)**2._dl)/n
  end function


  subroutine initialize_thermal_fluctuations(fld,len,m2,km,phi,kc,temp)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: len,m2,phi,temp
    integer, intent(in), optional :: km, kc

    real(dl), dimension(1:size(fld(:,1)/2+1)) :: spec, w2eff
    real(dl), dimension(1:size(fld(:,1))) :: df

    integer :: i
    real(dl) :: norm

    norm = 1._dl / phi / sqrt(2._dl * len)

    do i = 1, nyq
       w2eff(i) = m2 + dk**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:km) = norm / w2eff(2:km)**0.25 * sqrt( 2._dl/(exp(w2eff(2:km)**0.5/temp)-1._dl))
    call generate_1dGRF(df,spec(:kc))
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(:kc))
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_thermal_fluctuations


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


  subroutine output_fields(fld, dt, dtout, sim, temp)
    real(dl), dimension(1:nLat, 1:2) :: fld
    real(dl) :: dt, dtout, temp
    logical :: o; integer :: m, sim
    integer, parameter :: oFile = 98
    inquire(file='/gpfs/dpirvu/thermal_bubbles/phi4_m2_null_gam'//trim(real_str(gam))//'_t'//trim(str(nTime))//'_x'//trim(str(nLat))//'_temp'//trim(real_str(temp))//'_phi0'//trim(real_str(phi0))//'_lamb'//trim(real_str(lambda))//'_sim'//trim(str(sim))//'_fields.dat', opened=o)
    if (.not.o) then
       open(unit=oFile,file='/gpfs/dpirvu/thermal_bubbles/phi4_m2_null_gam'//trim(real_str(gam))//'_t'//trim(str(nTime))//'_x'//trim(str(nLat))//'_temp'//trim(real_str(temp))//'_phi0'//trim(real_str(phi0))//'_lamb'//trim(real_str(lambda))//'_sim'//trim(str(sim))//'_fields.dat')
       write(oFile,*) "# Lattice Parameters dx = ", dx, 'len = ', len
       write(oFile,*) "# Time Stepping parameters dt = ", dt, "dt_out = ", dtout
       write(oFile,*) "# Other Parameters m2eff = ", m2eff, "dk = ", dk
    endif

    do m = 1, nLat
       write(oFile,*) fld(m,:)
    enddo
  end subroutine output_fields

end program Gross_Pitaevskii_1d

