module gaussianRandomField
  use, intrinsic :: iso_c_binding
  use constants
  use fftw3
  implicit none
  logical :: seed_init = .false.
contains

  subroutine generate_1dGRF(field, spectrum, sym_in, initStride) !last two args are optional
    real(C_DOUBLE), dimension(:), intent(inout) :: field
    real(dl), dimension(:), intent(in) :: spectrum
    integer, optional :: sym_in, initStride
    integer :: symmetry, nLat, nn, stride, nnk !stride = 'a step or stage in progress toward an aim'
    logical :: cut
    real(dl), allocatable, dimension(:) :: amp, phase
    complex(dl), allocatable, dimension(:) :: deviate
    complex(C_DOUBLE_COMPLEX), allocatable :: Fk(:)
    type(C_PTR) :: fft_plan

    symmetry = 0
    if (present(sym_in)) symmetry = sym_in
    nLat = size(field); nn = nLat/2 + 1; nnk = size(spectrum); cut = .false.
    if (nn > nnk) then
       print*,"Warning spectrum is smaller than the number of required Fourier modes in 1dGRF.  Additional high frequency modes will not be sampled."
       cut = .true.
    endif
    if (.not.present(initStride)) then; stride = nnk; else; stride = initStride; endif
    if ( stride > nnk ) then; print*,"Stride is larger than number of Fourier modes.  Setting to nyquist"; stride = nnk; endif
    if (.not.seed_init) then
       print*,"Error, random number generator not initialized.  Call initialize_rand, using default seed values"
       call initialize_rand(75,13)
    endif
    allocate(Fk(1:nn))
    fft_plan = fftw_plan_dft_c2r_1d(nLat, Fk, field, FFTW_ESTIMATE)
    if (.not.seed_init) then
       print*,"Error, random number generator not initialized.  Calling initialize_rand, using default seed values"
       call initialize_rand(75,13)
    endif

    ! Generate Gaussian random deviates using Box-Muller
    ! Normalization chosen so that < Re(x)^2+Im(x)^2 > = 1
    ! generate random \rho and \phi for \psi = \sqrt(\rho) e^{i \phi}
    allocate( amp(1:nnk),phase(1:nnk),deviate(1:nnk) )
    call random_number(amp(1:stride)) 
    call random_number(phase(1:stride))
    if (stride < nnk) then
      call random_number( amp(stride+1:nnk) )
      call random_number( phase(stride+1:nnk) )
    endif
    select case (symmetry)
     case (1)
       deviate = dreal( sqrt(-2.*log(amp))*exp(iImag*twopi*phase) )
     case (2)
       deviate = iImag*dimag( sqrt(-2.*log(amp))*exp(iImag*twopi*phase) )
     case default
       deviate = sqrt(-log(amp))*exp(iImag*twopi*phase) 
    !construct gaussian random deviates with unit variance \alpha_n in paper
    end select
    
    ! Now produce the field with the given spectrum
    Fk(1) = 0._dl
    Fk(1:nnk) = deviate * spectrum
    Fk(nnk+1:nn) = 0._dl
    call fftw_execute_dft_c2r(fft_plan, Fk, field)

  end subroutine generate_1dGRF
  
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Warning, this subroutine isn't currently threadsafe
  subroutine initialize_rand(seed, seedfac)
    integer, intent(in) :: seed, seedfac
    integer :: nseed, i
    integer, allocatable, dimension(:) :: seeds

    seed_init=.true.
    call random_seed(size=nseed)
    allocate(seeds(1:nseed))
    seeds = seed + seedfac*(/ (i-1, i=1, nseed) /)
    call random_seed(put=seeds)
    deallocate(seeds)
  end subroutine initialize_rand
end module gaussianRandomField
