program makerotbmaps
  ! Make plane-of-sky maps of the integrated B-field through a simulation cube
  use wfitsutils, only: fitswrite
  use cuberotate, only: rotate, spin, &
       & interpolation, interpolation_nearest, interpolation_linear
  use mod_molfrac, only: molfrac
  implicit none
  integer :: i,  nx=0, ny, nz, nxx, nyy, nzz
  real, dimension(:,:,:), allocatable :: bmapi, bmapn, bmapm
  real, dimension(:,:), allocatable :: coli, coln, colm ! column densities
  real, dimension(:,:,:), allocatable :: d, xi, bx, by, bz, AV ! original datacubes
  real :: theta, phi
  real, dimension(:,:,:), allocatable :: dd, xxi, bbx, bby, bbz, AAV, mm ! rotated datacubes
  real, dimension(:,:,:), allocatable :: bbxx, bbyy, bbzz ! rotated vectors
  character(len=128) :: prefix
  integer :: itime
  character(len=4) :: ctime
  character(len=13) :: rotid
  real ::  zsize, dzz
  character(len=1), dimension(5) :: axid = (/ 'x', 'y', 'z', 'U', 'V' /)

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

  ! Visual extinction
  allocate(AAV(nxx,nyy,nzz))
  call read_cube(AV, 'AV')
  call rotate(AV, theta, phi, AAV)
  deallocate(AV)
  ! Molecular fraction
  allocate(mm(nxx,nyy,nzz))
  mm = molfrac(AAV)
  deallocate(AAV)

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
  allocate(bmapi(3, nxx, nyy), bmapn(3, nxx, nyy), bmapm(3, nxx, nyy))
  allocate(coli(nxx, nyy), coln(nxx, nyy), colm(nxx, nyy))

  ! B field weighted by ionized column
  bmapi(1, :, :) = sum(dd*xxi*bbxx, dim=3)
  bmapi(2, :, :) = sum(dd*xxi*bbyy, dim=3)
  bmapi(3, :, :) = sum(dd*xxi*bbzz, dim=3)
  bmapi(4, :, :) = sum(dd*xxi*abs(bbxx), dim=3)
  bmapi(5, :, :) = sum(dd*xxi*abs(bbyy), dim=3)

  ! B field weighted by neutral column
  bmapn(1, :, :) = sum(dd*(1.-xxi)*(1.-mm)*bbxx, dim=3)
  bmapn(2, :, :) = sum(dd*(1.-xxi)*(1.-mm)*bbyy, dim=3)
  bmapn(3, :, :) = sum(dd*(1.-xxi)*(1.-mm)*bbzz, dim=3)
  bmapn(4, :, :) = sum(dd*(1.-xxi)*(1.-mm)*abs(bbxx), dim=3)
  bmapn(5, :, :) = sum(dd*(1.-xxi)*(1.-mm)*abs(bbyy), dim=3)

  ! B field weighted by molecular column
  bmapm(1, :, :) = sum(dd*(1.-xxi)*mm*bbxx, dim=3)
  bmapm(2, :, :) = sum(dd*(1.-xxi)*mm*bbyy, dim=3)
  bmapm(3, :, :) = sum(dd*(1.-xxi)*mm*bbzz, dim=3)
  bmapm(4, :, :) = sum(dd*(1.-xxi)*mm*abs(bbxx), dim=3)
  bmapm(5, :, :) = sum(dd*(1.-xxi)*mm*abs(bbyy), dim=3)

  ! ... and the columns themselves
  coli(:, :) = sum(dd*xxi, dim=3)
  coln(:, :) = sum(dd*(1.-xxi)*(1.-mm), dim=3)
  colm(:, :) = sum(dd*(1.-xxi)*mm, dim=3)

  do i = 1, 5
     call fitswrite(bmapi(i, :, :), &
          & trim(prefix)//'-bmap-i-'//axid(i)//rotid//ctime//'.fits')
     call fitswrite(bmapn(i, :, :), &
          & trim(prefix)//'-bmap-n-'//axid(i)//rotid//ctime//'.fits')
     call fitswrite(bmapm(i, :, :), &
          & trim(prefix)//'-bmap-m-'//axid(i)//rotid//ctime//'.fits')
  end do
  call fitswrite(coli, &
       & trim(prefix)//'-colmap-i'//rotid//ctime//'.fits')
  call fitswrite(coln, &
       & trim(prefix)//'-colmap-n'//rotid//ctime//'.fits')
  call fitswrite(colm, &
       & trim(prefix)//'-colmap-m'//rotid//ctime//'.fits')
  
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
