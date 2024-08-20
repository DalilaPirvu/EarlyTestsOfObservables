!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS

#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
  use fftw3
  implicit none

  integer, parameter :: nFld = 1, nLat = 128, nyq = nLat/2+1, nVar = 2*nFld*nLat+1
  real(dl), dimension(1:nVar), target :: yvec
  real(dl), parameter :: m2eff = 1._dl, lambda = 1._dl, gam = 0.25_dl
  real(dl), parameter :: len = 1._dl, dx = len/dble(nLat), dk = twopi/len
  type(transformPair1D) :: tPair

contains
  
  subroutine derivs(yc, yp, lambda, m2eff, gam)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), intent(in) :: lambda, m2eff, gam
    real(dl), dimension(:), intent(out) :: yp
    integer :: x
 
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)

    !do x = 1, nLat
    !    if (yc(x) < 0) then
    !        yp(nLat+x) = + m2eff*yc(x) - lambda * yc(x)**3._dl / 6._dl + gam * yc(x)**2._dl / 2._dl
    !    else
    !        yp(nLat+x) = + m2eff*yc(x) - lambda * yc(x)**3._dl / 6._dl
    !    endif
    !enddo
    
    yp(DFLD) = + m2eff*yc(FLD) - lambda * yc(FLD)**3._dl / 6._dl + gam ! damped gamma = 0.25
    
    ! correction gamma = 0.5

    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
  end subroutine derivs

end module eom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS
