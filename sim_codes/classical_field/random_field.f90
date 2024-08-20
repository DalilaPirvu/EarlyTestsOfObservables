module gaussianRandomField
  use, intrinsic :: iso_c_binding
  use constants
  use fftw3
  implicit none
contains

  subroutine generate_1dGRF(field, spectrum)
    real(C_DOUBLE), dimension(:), intent(inout) :: field
    real(dl), dimension(:), intent(in) :: spectrum
    integer :: nLat, nn, nnk
    complex(C_DOUBLE_COMPLEX), allocatable :: Fk(:)
    type(C_PTR) :: fft_plan

    nLat = size(field); nn = nLat/2 + 1; nnk = size(spectrum)
    if (nn > nnk) then
       print*,"Warning spectrum is smaller than the number of required Fourier modes in 1dGRF.  Additional high frequency modes will not be sampled."
    endif
    allocate(Fk(1:nn))
    fft_plan = fftw_plan_dft_c2r_1d(nLat, Fk, field, FFTW_ESTIMATE)

    Fk(1) = 0._dl
    Fk(1:nnk) = spectrum
    Fk(nnk+1:nn) = 0._dl
    call fftw_execute_dft_c2r(fft_plan, Fk, field)

  end subroutine generate_1dGRF
end module gaussianRandomField
