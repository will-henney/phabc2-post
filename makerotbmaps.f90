program makerotbmaps
  ! Make plane-of-sky maps of the integrated B-field through a simulation cube
  use wfitsutils, only: fitswrite
  use cuberotate, only: rotate, spin, &
       & interpolation, interpolation_nearest, interpolation_linear
  implicit none
  integer :: i,  nx=0, ny, nz, nxx, nyy, nzz
  real, dimension(:,:,:), allocatable :: bmapi, bmapn
  real, dimension(:,:,:), allocatable :: d, xi, bx, by, bz ! original datacubes
  real :: theta, phi
  real, dimension(:,:,:), allocatable :: dd, xxi, bbx, bby, bbz ! rotated datacubes
  real, dimension(:,:,:), allocatable :: bbxx, bbyy, bbzz ! rotated vectors
  character(len=128) :: prefix
  integer :: itime
  character(len=4) :: ctime
  character(len=13) :: rotid
  real ::  zsize, dzz
  character(len=1), dimension(3) :: axid = (/ 'x', 'y', 'z' /)

  ! options are nearest/linear
  interpolation = interpolation_linear

  print *, 'File prefix?'
  read '(a)', prefix

  print *, 'Save time?'
  read *, itime
  write(ctime, '(i4.4)') itime
  
  print *, 'Cube extent in z-direction (parsec)?'
  read *, zsize

  print *, 'Rotation theta, phi?'
  read *, theta, phi
  write(rotid, '(a,2(sp,i4.3),a)') '-rot', int(theta), int(phi), '-'

  call read_cube(d, 'dd')
  print '("Original grid size: ",i0,"x",i0,"x",i0)', shape(d)
  d = d/(1.3*1.67262158e-24)    ! convert to cm^-3
  dzz = zsize*3.085677582e18/real(nz) ! physical pixel scale

  ! WJH 17 Sep 2010 - port over the pre-allocation part from makerotmap.f90
  ! Pre-allocate output cube to prevent truncation along zz axis
  nxx = nx
  nyy = ny
  nzz = int(sqrt(real(nx**2 + ny**2 + nz**2))) ! longest diagonal possible
  allocate(dd(nxx,nyy,nzz))

  call rotate(d, theta, phi, dd)
  deallocate(d)
  print '("Rotated grid size: ",i0,"x",i0,"x",i0)', nxx, nyy, nzz

  ! Ionization fraction
  allocate(xxi(nxx,nyy,nzz))
  call read_cube(xi, 'xi')
  call rotate(xi, theta, phi, xxi)
  deallocate(xi)

  ! B field
  allocate(bbx(nxx,nyy,nzz))
  call read_cube(bx, 'bx')
  call rotate(bx, theta, phi, bbx)
  deallocate(bx)

  allocate(bby(nxx,nyy,nzz))
  call read_cube(by, 'by')
  call rotate(by, theta, phi, bby)
  deallocate(by)

  allocate(bbz(nxx,nyy,nzz))
  call read_cube(bz, 'bz')
  call rotate(bz, theta, phi, bbz)
  deallocate(bz)

  ! Rotate the vector of B to the new coordinate frame
  allocate(bbxx(nxx,nyy,nzz), bbyy(nxx,nyy,nzz), bbzz(nxx,nyy,nzz))
  call spin(bbx, bby, bbz, theta, phi, bbxx, bbyy, bbzz) 
  deallocate(bbx, bby, bbz)

 
  ! integration is along the new zz-axis
  allocate(bmapi(3, nxx, nyy), bmapn(3, nxx, nyy))

  bmapi(1, :, :) = sum(dd*xxi*bbxx, dim=3)
  bmapi(2, :, :) = sum(dd*xxi*bbyy, dim=3)
  bmapi(3, :, :) = sum(dd*xxi*bbzz, dim=3)

  bmapn(1, :, :) = sum(dd*(1.-xxi)*bbxx, dim=3)
  bmapn(2, :, :) = sum(dd*(1.-xxi)*bbyy, dim=3)
  bmapn(3, :, :) = sum(dd*(1.-xxi)*bbzz, dim=3)

  do i = 1, 3
     call fitswrite(bmapi(i, :, :), &
          & trim(prefix)//'-bmap-i-'//axid(i)//'-'//rotid//ctime//'.fits')
     call fitswrite(bmapn(i, :, :), &
          & trim(prefix)//'-bmap-n-'//axid(i)//'-'//rotid//ctime//'.fits')
  end do
  
contains

  subroutine read_cube(var, id)
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
    deallocate(fitscube)
  end subroutine read_cube

 
end program makerotbmaps
