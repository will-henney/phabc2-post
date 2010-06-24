program makerotbmaps
  ! Make plane-of-sky maps of the integrated B-field through a simulation cube
  use wfitsutils, only: fitswrite
  use cuberotate, only: rotate, spin, &
       & interpolation, interpolation_nearest, interpolation_linear
  implicit none
  integer :: i, j, k, nx=0, ny, nz, nxx, nyy, nzz
  real, dimension(:), allocatable :: tau
  real, dimension(:,:), allocatable :: map
  real, dimension(:,:,:), allocatable :: xi, bx, by, bz ! original datacubes
  real :: theta, phi
  real, dimension(:,:,:), allocatable :: xxi, bbx, bby, bbz ! rotated datacubes
  real, dimension(:,:,:), allocatable :: bbxx, bbyy, bbzz ! rotated vectors
  character(len=128) :: prefix
  integer :: itime
  character(len=4) :: ctime
  character(len=13) :: rotid
  real :: extinct_sigma, zsize, dzz

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
  call rotate(d, theta, phi, dd)
  deallocate(d)
  nxx = size(dd, 1)
  nyy = size(dd, 2)
  nzz = size(dd, 3)
  print '("Rotated grid size: ",i0,"x",i0,"x",i0)', nxx, nyy, nzz

  ! Ionization fraction
  call read_cube(xi, 'xi')
  call rotate(xi, theta, phi, xxi)
  deallocate(xi)

  ! B field
  call read_cube(bx, 'bx')
  call rotate(bx, theta, phi, bbx)
  deallocate(bx)
  call read_cube(by, 'by')
  call rotate(by, theta, phi, bby)
  deallocate(by)
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

  call fitswrite(map, trim(prefix)//'map'//rotid//emtype//ctime//'.fits')

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

  subroutine domap
    do i = 1, nx
       do j = 1, ny
          tau(1) = 0.0
          do k = 2, nz
             tau(k) = tau(k-1) + dd(i,j,k) 
          end do
          tau = tau*dzz*extinct_sigma
          map(i, j) = sum(ee(i,j,:)*exp(-tau))
       end do
    end do
  end subroutine domap

end program makerotbmaps
