program cubet
  ! calculates cube of temperature only
  implicit none
  real, dimension(:,:,:), allocatable :: dd, xi, pp, te
  character(len=128) :: prefix
  real, parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  integer :: nx = 0, ny = 0, nz = 0
  integer :: itime
  character(len=4) :: ctime
  integer, parameter :: ilogging = 1

  print *, 'File prefix?'
  read '(a)', prefix

  print *, 'Save time?'
  read *, itime

  write(ctime, '(i4.4)') itime

  call read_var(dd, 'dd')
  call read_var(xi, 'xi')
  call read_var(pp, 'pp')

  ! Temperature
  allocate( te(nx, ny, nz) )
  te = (mp * mu / dd) * pp / (1.+xi) / boltzmann_k 
  call write_var(te, 'te')

  if (ilogging > 0) then
     print '(3(a,es10.2))', 'Temperature min/max/mean: ', &
          & minval(te), '/', maxval(te), '/', sum(te)/real(nx*ny*nz)
  end if
  
contains

  subroutine read_var(var, id)
    use wfitsutils, only: fitsread, fitscube
    real, intent(inout), dimension(:,:,:), allocatable :: var
    character(len=2), intent(in) :: id
    call fitsread(trim(prefix)//'-'//id//ctime//'.fits')
    if (nx==0) then
       nx = size(fitscube, 1)
       ny = size(fitscube, 2)
       nz = size(fitscube, 3)
    end if
    allocate(var(nx,ny,nz))
    var = fitscube
  end subroutine read_var

  subroutine write_var(var, id)
    use wfitsutils, only: fitswrite
    real, intent(in), dimension(:,:,:) :: var
    character(len=2), intent(in) :: id
    call fitswrite(var, trim(prefix)//'-'//id//ctime//'.fits')
  end subroutine write_var


end program cubet
