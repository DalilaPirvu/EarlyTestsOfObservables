!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS

#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
  use fftw3
  implicit none

  integer, parameter :: nFld = 1, nLat = 512, nyq = nLat/2+1, nVar = 2*nFld*nLat+1
  real(dl), dimension(1:nVar), target :: yvec
  real(dl), parameter :: len = 100._dl, dx = len/dble(nLat), dk = twopi/len
  type(transformPair1D) :: tPair

contains
  
  subroutine derivs(yc, yp, lamb, m2, gam)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), intent(in) :: lamb, m2, gam
    real(dl), dimension(:), intent(out) :: yp
 
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)

    !if (gam == 0) then ! positive effective mass; symmetry restored
    !    yp(DFLD) = - m2*yc(FLD) - lamb * yc(FLD)**3._dl / 6._dl
    !elseif (gam == 1) then ! negative effective mass; broken phase
    !    yp(DFLD) = + m2*yc(FLD) - lamb * yc(FLD)**3._dl / 6._dl
    !elseif (gam == 2) then ! quadratic potential; harmonic motion; free field
    !    yp(DFLD) = - m2*yc(FLD)
    !endif
    
    !!!!!!!!!!!! Other potentials

    !do x = 1, nLat
    !    if (yc(x) < 0) then
    !        yp(nLat+x) = + m2eff*yc(x) - lamb * yc(x)**3._dl / 6._dl + gam * yc(x)**2._dl / 2._dl
    !    else
    !        yp(nLat+x) = + m2eff*yc(x) - lamb * yc(x)**3._dl / 6._dl
    !    endif
    !enddo

    yp(DFLD) = + m2*yc(FLD) - lamb * yc(FLD)**3._dl / 6._dl + gam ! damped

    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
  end subroutine derivs

end module eom
