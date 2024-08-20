!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS

#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
  use fftw3
  implicit none

  integer, parameter :: nFld = 1, nLat = 1024, nyq = nLat/2+1, nVar = 2*nFld*nLat+1
  real(dl), dimension(1:nVar), target :: yvec
  real(dl), parameter :: lambda = 1._dl, gam = 0._dl
  real(dl), parameter :: len = 100._dl, dx = len/dble(nLat), dk = twopi/len
  type(transformPair1D) :: tPair

contains
  
  subroutine derivs(yc, yp, lambda, m2, gam)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), intent(in) :: lambda, gam, m2
    real(dl), dimension(:), intent(out) :: yp
 
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)

    yp(DFLD) = - m2*yc(FLD) - lambda * yc(FLD)**3._dl / 6._dl

    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
  end subroutine derivs
  
end module eom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS
