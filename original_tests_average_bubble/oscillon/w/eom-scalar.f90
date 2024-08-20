!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS

#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
  use fftw3
  implicit none

  integer, parameter :: nFld = 1, nLat = 512, nyq = nLat/2+1, nVar = 2*nFld*nLat+1
  real(dl), dimension(1:nVar), target :: yvec
  integer, parameter :: typeP = 1, AA = 10, ppp = 3
  real(dl), parameter :: nu = 2.e-3
  real(dl), parameter :: lambda = 6._dl
  real(dl), parameter :: m2eff = 4._dl*nu*(-1._dl+lambda**2)
  real(dl), parameter :: phimax = 1.59858
  real(dl), parameter :: lenLat = 50._dl / (2.*nu)**0.5
  real(dl), parameter :: dx = lenLat/dble(nLat), dk = twopi/lenLat
!  real(dl), parameter :: dx = 1.5440808887540916, dk = 0.007947670612636881
  type(transformPair1D) :: tPair

contains
  
  subroutine derivs(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)

    if (typeP == 1) then ! lamb = 6
        yp(DFLD) = -4._dl*nu*( sin(yc(FLD)) + 0.5_dl*lambda**2*sin(2._dl*yc(FLD)) )
    
!    else if (typeP == 2) then
!        yp(DFLD) = -0.004*(1+36.*COS(yc(FLD)))*(2.+DTANH(-132.115+25.*yc(FLD)) - DTANH(-117.115+25.*yc(FLD)) + DTANH(-39.9644+25.*yc(FLD)) - DTANH(-24.9644+25.*yc(FLD)))*SIN(yc(FLD)) - 0.1*(COS(yc(FLD)) - 18.*SIN(yc(FLD))**2) * (COSH(24.9644-25.*yc(FLD))**-2 - COSH(39.9644-25.*yc(FLD))**-2 + COSH(117.115-25.*yc(FLD))**-2 - COSH(132.115-25.*yc(FLD))**-2) + 1.80139*COSH(25*(-5.28461+yc(FLD)))**-2 - 1.80139*COSH(25*(-4.68461+yc(FLD)))**-2 + 1.80139*COSH(25*(-1.59858+yc(FLD)))**-2 - 1.80139*COSH(25*(-0.998578+yc(FLD)))**-2
    
!    else if (typeP == 3) then
!        yp(DFLD) = -0.004*(1+36.*COS(yc(FLD)))*(2.+DTANH(-117.115+25.*yc(FLD)) - DTANH(-102.115+25.*yc(FLD)) + DTANH(-54.9644+25.*yc(FLD)) - DTANH(-39.9644+25.*yc(FLD)))*SIN(yc(FLD)) - 0.1*(COS(yc(FLD)) - 18.*SIN(yc(FLD))**2) * (COSH(39.9644-25.*yc(FLD))**-2 - COSH(54.9644-25.*yc(FLD))**-2 + COSH(102.115-25.*yc(FLD))**-2 - COSH(117.115-25.*yc(FLD))**-2) + 1.80139*COSH(25*(-4.68461+yc(FLD)))**-2 - 1.80139*COSH(25*(-4.08461+yc(FLD)))**-2 + 1.80139*COSH(25*(-2.19858+yc(FLD)))**-2 - 1.80139*COSH(25*(-1.59858+yc(FLD)))**-2
    
!    else if (typeP == 4) then    
!        yp(DFLD) = -0.004*(1+36.*COS(yc(FLD)))*(2.+DTANH(-132.115+25.*yc(FLD)) - DTANH(-102.115+25.*yc(FLD)) + DTANH(-54.9644+25.*yc(FLD)) - DTANH(-24.9644+25.*yc(FLD)))*SIN(yc(FLD)) - 0.1*(COS(yc(FLD)) - 18.*SIN(yc(FLD))**2) * (COSH(24.9644-25.*yc(FLD))**-2 - COSH(54.9644-25.*yc(FLD))**-2 + COSH(102.115-25.*yc(FLD))**-2 - COSH(132.115-25.*yc(FLD))**-2) + 1.80139*COSH(25*(-5.28461+yc(FLD)))**-2 - 1.80139*COSH(25*(-4.08461+yc(FLD)))**-2 + 1.80139*COSH(25*(-2.19858+yc(FLD)))**-2 - 1.80139*COSH(25*(-0.998578+yc(FLD)))**-2

    else if (typeP == 6) then
        where (yc(FLD) >= phimax)
            yp(DFLD) = -4._dl*nu*( sin(yc(FLD)) + 0.5_dl*lambda**2*sin(2._dl*yc(FLD)) )
        elsewhere
            yp(DFLD) = 4._dl*nu* AA * ABS(yc(FLD) - phimax)**(ppp-1.)
        end where
    end if

    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
  end subroutine derivs

end module eom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUBBLE CORRELATIONS
