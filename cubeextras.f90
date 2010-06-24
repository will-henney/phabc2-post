program cubeextras
  ! calculates cubes of various derived quantities
  ! cubes of density, ionization fraction, pressure, b field, velocity
  implicit none
  real, dimension(:,:,:), allocatable :: dd, xi, pp, te
  real, dimension(:,:,:), allocatable :: bx, by, bz, bb, be
  real, dimension(:,:,:), allocatable :: vx, vy, vz, mm, av, vv
!!$  real, dimension(:,:,:), allocatable :: x, y, z, r
  character(len=128) :: prefix
  real, parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  integer :: nx = 0, ny = 0, nz = 0
  integer :: i, j, k, itime
  character(len=4) :: ctime

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
  te = mp * mu * pp / dd / (1.+xi) / boltzmann_k 
  call write_var(te, 'te')

  deallocate(xi, te)

  ! B field
  call read_var(bx, 'bx')
  call read_var(by, 'by')
  call read_var(bz, 'bz')
  allocate( bb(nx, ny, nz), be(nx, ny, nz) )

  ! magnetic pressure
  bb = 0.5*(bx**2 + by**2 + bz**2)
  ! plasma beta
  be = pp/bb
  call write_var(bb, 'bb')
  call write_var(be, 'be')

  deallocate(bx, by, bz, be)

!!$  ! Position
!!$  allocate( x(nx, ny, nz), y(nx, ny, nz), z(nx, ny, nz), r(nx, ny, nz) )
!!$  forall(i=1:nx, j=1:ny, k=1:nz)
!!$     ! cartesian distances from the center
!!$     x(i,j,k) = real(i) - 0.5*real(nx+1)
!!$     y(i,j,k) = real(j) - 0.5*real(ny+1)
!!$     z(i,j,k) = real(k) - 0.5*real(nz+1)
!!$  end forall
!!$  r = sqrt(x**2 + y**2 + z**2)

  ! Velocities
  call read_var(vx, 'vx')
  call read_var(vy, 'vy')
  call read_var(vz, 'vz')

  allocate( vv(nx, ny, nz), mm(nx, ny, nz), av(nx, ny, nz) )
  ! KE density
  vv = 0.5*dd*(vx**2 + vy**2 + vz**2)
  ! Mach number squared
  mm = 2.0*vv/pp
  ! Alfven speed
  av = sqrt(2.0*bb/dd)

  call write_var(mm, 'mm')
  call write_var(vv, 'vv')
  ! Mac filesystem is case-insensitive, so we can't have AV and av
  ! WJH 25 Feb 2009
  call write_var(av, 'va')      

!!$  ! radial component of velocity
!!$  vr = (vx*x + vy*y + vz*z)/r
!!$  ! lateral component of velocity
!!$  vl = sqrt(vx**2 + vy**2 + vz**2 - vr**2)
!!$  call write_var(vr, 'vr')
!!$  call write_var(vl, 'vl')


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


end program cubeextras
