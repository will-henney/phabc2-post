program rotatetest
  use cuberotate, only: rotate, &
       & interpolation, interpolation_nearest, interpolation_linear
  implicit none
  character(len=128) :: prefix
  integer :: itime
  character(len=4) :: ctime
  character(len=8) :: varid
  character(len=11) :: rotid
  real :: theta, phi
  integer :: nx, ny, nz
  real, allocatable, dimension(:,:,:) :: original_cube, rotated_cube

  ! options are nearest/linear
!   interpolation = interpolation_nearest
  interpolation = interpolation_linear

  print *, 'File prefix?'
  read '(a)', prefix

  print *, 'Save time?'
  read *, itime
  write(ctime, '(i4.4)') itime

  print *, 'Variable ID?'
  read *, varid

  print *, 'Rotation theta, phi?'
  read *, theta, phi
  write(rotid, '(a,2i3.3,a)') '-rot', int(theta), int(phi), '-'

  call read_var(original_cube, trim(varid))
  call rotate(original_cube, theta, phi, rotated_cube)
  call write_var(rotated_cube, trim(varid)//rotid)

contains

  subroutine read_var(var, id)
    use wfitsutils, only: fitsread, fitscube
    real, intent(inout), dimension(:,:,:), allocatable :: var
    character(len=*), intent(in) :: id
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
    character(len=*), intent(in) :: id
    call fitswrite(var, trim(prefix)//'-'//id//ctime//'.fits')
  end subroutine write_var

end program rotatetest
