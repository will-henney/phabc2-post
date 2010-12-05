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
  real, dimension(:,:,:), allocatable :: bq, bu, alpha, bbxy           ! vectors in q-u space
  real, dimension(:,:), allocatable :: bmapiq, bmapnq, bmapmq
  real, dimension(:,:), allocatable :: bmapiu, bmapnu, bmapmu, qu, al
  character(len=128) :: prefix
  integer :: itime
  character(len=4) :: ctime
  character(len=13) :: rotid
  real ::  zsize, dzz
  integer, parameter :: naxes = 5
  character(len=1), dimension(naxes) :: axid = (/ 'x', 'y', 'z', 'U', 'V' /)

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
  allocate(bmapi(naxes, nxx, nyy), bmapn(naxes, nxx, nyy), bmapm(naxes, nxx, nyy))
  allocate(coli(nxx, nyy), coln(nxx, nyy), colm(nxx, nyy))

  ! B field weighted by ionized column
  bmapi(1, :, :) = sum(dd*xxi*bbxx, dim=3)
  bmapi(2, :, :) = sum(dd*xxi*bbyy, dim=3)
  bmapi(3, :, :) = sum(dd*xxi*bbzz, dim=3)

  ! B field weighted by neutral column
  bmapn(1, :, :) = sum(dd*(1.-xxi)*(1.-mm)*bbxx, dim=3)
  bmapn(2, :, :) = sum(dd*(1.-xxi)*(1.-mm)*bbyy, dim=3)
  bmapn(3, :, :) = sum(dd*(1.-xxi)*(1.-mm)*bbzz, dim=3)

  ! B field weighted by molecular column
  bmapm(1, :, :) = sum(dd*(1.-xxi)*mm*bbxx, dim=3)
  bmapm(2, :, :) = sum(dd*(1.-xxi)*mm*bbyy, dim=3)
  bmapm(3, :, :) = sum(dd*(1.-xxi)*mm*bbzz, dim=3)

  ! ... and the columns themselves
  coli(:, :) = sum(dd*xxi, dim=3)
  coln(:, :) = sum(dd*(1.-xxi)*(1.-mm), dim=3)
  colm(:, :) = sum(dd*(1.-xxi)*mm, dim=3)

  ! Now transform B_x, B_y to B_q, B_u
  ! B_q = B_x - B_y
  ! B_q^2 + B_u^2 = B_x^2 + B_y^2 = B_x^2 + B_y^2 - 2 B_x B_y + B_u^2
  ! B_u = +/- sqrt(2 B_x B_y) = sgn(By) sqrt(2 B_x B_y)
  allocate(bq(nxx,nyy,nzz), bu(nxx,nyy,nzz), bbxy(nxx,nyy,nzz), alpha(nxx,nyy,nzz))
  alpha = atan2(bbyy, bbxx)
  bbxy = sqrt(bbxx**2 + bbyy**2)
  bq = bbxy*cos(2.0*alpha)
  bu = bbxy*sin(2.0*alpha)
  deallocate(bbxy, alpha)

!!$  bu = sign(sqrt(2.0*bbxx*bbyy), bbyy)
  
  print '(a)', 'whole-cube statistics in microG [min/mean/max]:'
  print '(a,"[",f0.3,"/",f0.3,"/",f0.3,"]")', 'bbxx: ', &
       & minval(1.e6*bbxx), sum(1.e6*bbxx)/product(shape(bbxx)), maxval(1.e6*bbxx) 
  print '(a,"[",f0.3,"/",f0.3,"/",f0.3,"]")', 'bbyy: ', &
       & minval(1.e6*bbyy), sum(1.e6*bbyy)/product(shape(bbyy)), maxval(1.e6*bbyy) 
  print '(a,"[",f0.3,"/",f0.3,"/",f0.3,"]")', 'bbzz: ', &
       & minval(1.e6*bbzz), sum(1.e6*bbzz)/product(shape(bbzz)), maxval(1.e6*bbzz) 
  print '(a,"[",f0.3,"/",f0.3,"/",f0.3,"]")', 'bq(real): ', &
       & minval(1.e6*real(bq)), sum(1.e6*real(bq))/product(shape(bq)), maxval(1.e6*real(bq)) 
!!$  print '(a,"[",f0.3,"/",f0.3,"/",f0.3,"]")', 'bq(imag): ', &
!!$       & minval(1.e6*imag(bq)), sum(1.e6*imag(bq))/product(shape(bq)), maxval(1.e6*imag(bq)) 
  print '(a,"[",f0.3,"/",f0.3,"/",f0.3,"]")', 'bu(real): ', &
       & minval(1.e6*real(bu)), sum(1.e6*real(bu))/product(shape(bu)), maxval(1.e6*real(bu)) 
!!$  print '(a,"[",f0.3,"/",f0.3,"/",f0.3,"]")', 'bu(imag): ', &
!!$       & minval(1.e6*imag(bu)), sum(1.e6*imag(bu))/product(shape(bu)), maxval(1.e6*imag(bu)) 


  ! do the integration of the q, u components
  allocate(bmapiq(nxx, nyy), bmapnq(nxx, nyy), bmapmq(nxx, nyy))
  allocate(bmapiu(nxx, nyy), bmapnu(nxx, nyy), bmapmu(nxx, nyy))
  bmapiq = sum(dd*xxi*bq, dim=3)
  print *, 'q', 'i', maxval(real(bmapiq)), minval(real(bmapiq))
  bmapiu = sum(dd*xxi*bu, dim=3)
  print *, 'u', 'i', maxval(real(bmapiu)), minval(real(bmapiu))
  bmapnq = sum(dd*(1.-xxi)*(1.-mm)*bq, dim=3)
  print *, 'q', 'n', maxval(bmapnq), minval(bmapnq) 
  bmapnu = sum(dd*(1.-xxi)*(1.-mm)*bu, dim=3)
  print *, 'u', 'n', maxval(bmapnu), minval(bmapnu) 
  bmapmq = sum(dd*(1.-xxi)*mm*bq, dim=3)
  print *, 'q', 'm', maxval(bmapmq), minval(bmapmq) 
  bmapmu = sum(dd*(1.-xxi)*mm*bu, dim=3)
  print *, 'u', 'm', maxval(bmapmu), minval(bmapmu) 
  
  ! now switch back to x, y
  allocate(qu(nxx, nyy), al(nxx, nyy))
  qu = sqrt(bmapiq**2 + bmapiu**2)
  al = 0.5*atan2(bmapiu, bmapiq)
  bmapi(4,:,:) = qu*cos(al)
  bmapi(5,:,:) = qu*sin(al)
  qu = sqrt(bmapnq**2 + bmapnu**2)
  al = 0.5*atan2(bmapnu, bmapnq)
  bmapn(4,:,:) = qu*cos(al)
  bmapn(5,:,:) = qu*sin(al)
  qu = sqrt(bmapmq**2 + bmapmu**2)
  al = 0.5*atan2(bmapmu, bmapmq)
  bmapm(4,:,:) = qu*cos(al)
  bmapm(5,:,:) = qu*sin(al)
  
  do i = 1, naxes
     call fitswrite(bmapi(i, :, :), &
          & trim(prefix)//'-bmap-i-'//axid(i)//rotid//ctime//'.fits')
     call fitswrite(bmapn(i, :, :), &
          & trim(prefix)//'-bmap-n-'//axid(i)//rotid//ctime//'.fits')
     call fitswrite(bmapm(i, :, :), &
          & trim(prefix)//'-bmap-m-'//axid(i)//rotid//ctime//'.fits')
     print *, axid(i), 'i', maxval(bmapi(i, :, :)), minval(bmapi(i, :, :))
     print *, axid(i), 'n', maxval(bmapn(i, :, :)), minval(bmapn(i, :, :))
     print *, axid(i), 'm', maxval(bmapm(i, :, :)), minval(bmapm(i, :, :))
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
