program globstats
  use modstats, only: find_zero
  use allstats_vars
  implicit none
  real, dimension(:,:,:), allocatable :: dd, xi, pp !, te
  real, dimension(:,:,:), allocatable :: bx, by, bz !, bb, be
  real, dimension(:,:,:), allocatable :: vx, vy, vz !, mm, av

  character(len=128) :: prefix
  real, parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  real, parameter :: pc = 3.085677582e18, msun = 1.989e33, km = 1.0e5
  real, parameter :: pi = 3.1415926535897932385e0
  real, parameter :: microgauss = 1.e-6/3.54490770182 ! sqrt(4.0*pi) 
  integer :: nx = 0, ny = 0, nz = 0
  integer :: i, j, k, itime
  character(len=4) :: ctime
  ! Size of simulation box in pc - hardcoded for simplicity
  real, parameter :: xmax = 2.0, ymax = 1.0, zmax = 1.0
  real :: dx, dy, dz            ! cell size along each axis
  ! length along each axis, in pc
  real, allocatable, dimension(:) :: x, y, z 
  ! Masses per unit length
  real, allocatable, dimension(:) :: lamx, lamy, lamz
  ! extra arrays that come in useful
  real, allocatable, dimension(:,:,:) :: dn
  real, allocatable, dimension(:) :: meanxi, meanvx
  integer :: j0, k0, j1, k1, i0, ishell
  real :: ay, az, a00, a10, a01, a11
  real :: cell_volume

  print *, 'File prefix?'
  read '(a)', prefix

  print *, 'Save time?'
  read *, itime
  write(ctime, '(i4.4)') itime

  call read_var(dd, 'dd')
  call read_var(xi, 'xi')
  call read_var(pp, 'pp')

  call read_var(bx, 'bx')
  call read_var(by, 'by')
  call read_var(bz, 'bz')
  ! put B field in microGauss
  bx = bx/microgauss
  by = by/microgauss
  bz = bz/microgauss


  call read_var(vx, 'vx')
  call read_var(vy, 'vy')
  call read_var(vz, 'vz')
  ! put velocities in km/s
  vx = vx/km
  vy = vy/km
  vz = vz/km

  dx = xmax/real(nx)
  dy = ymax/real(ny)
  dz = zmax/real(nz)

  ! Do the stats
  cell_volume = dx*dy*dz * (pc/msun)*pc*pc !avoid overflow

  ! We use the neutral mass a lot, so make a permanent array of it
  allocate(dn(nx,ny,nz))
  dn = dd*(1.0-xi)

  ! Total masses in Msun
  neutral_mass = sum(dn) * cell_volume
  ionized_mass = sum(dd*xi) * cell_volume

  ! Center of mass along each axis
  allocate(x(nx), lamx(nx), y(ny), lamy(ny), z(nz), lamz(nz))

!!$  x = (/ (dx*(real(i) - 0.5), i = 1, nx) /)
!!$  y = (/ (dy*(real(j) - 0.5), j = 1, ny) /)
!!$  z = (/ (dz*(real(k) - 0.5), k = 1, nz) /)
  ! should we take off 0.5 here or not? 
  ! We seem to get the globule COM being 0.5 more accurately if we don't
  x = (/ (dx*real(i), i = 1, nx) /)
  y = (/ (dy*real(j), j = 1, ny) /)
  z = (/ (dz*real(k), k = 1, nz) /)
  
  ! neutral mass per unit length along each axis
  lamx = sum(sum(dn, dim=3), dim=2)
  lamy = sum(sum(dn, dim=3), dim=1)
  lamz = sum(sum(dn, dim=1), dim=1)

  xcom = sum(x*lamx)/sum(lamx)
  ycom = sum(y*lamy)/sum(lamy)
  zcom = sum(z*lamz)/sum(lamz)

  ! mass-weighted mean velocity of neutral gas in km/s
  vnx = sum(vx*dn)/sum(dn)
  vny = sum(vy*dn)/sum(dn)
  vnz = sum(vz*dn)/sum(dn)

  ! Magnetic field 
  ! Total mean B
  bmeanx = sum(bx)/real(nx*ny*nz)
  bmeany = sum(by)/real(nx*ny*nz)
  bmeanz = sum(bz)/real(nx*ny*nz)
  ! Ionized mean B
  bmeanx_i = sum(bx*xi)/sum(xi)
  bmeany_i = sum(by*xi)/sum(xi)
  bmeanz_i = sum(bz*xi)/sum(xi)
  ! Neutral mean B
  bmeanx_n = sum(bx*(1.0-xi))/sum(1.0-xi)
  bmeany_n = sum(by*(1.0-xi))/sum(1.0-xi)
  bmeanz_n = sum(bz*(1.0-xi))/sum(1.0-xi)

  ! Total RMS B
  brms = sqrt(sum(bx**2 + by**2 + bz**2)/real(nx*ny*nz))

  ! mean betas
  beta = (2.0/microgauss**2)*sum(pp)/sum(bx**2 + by**2 + bz**2)
  beta_i = (2.0/microgauss**2)*sum(pp*xi)/sum(xi*(bx**2 + by**2 + bz**2))
  beta_n = (2.0/microgauss**2)*sum(pp*(1.-xi))/sum((1.-xi)*(bx**2 + by**2 + bz**2))

  ! Position of i-front along x-axis through globule COM
  allocate( meanxi(nx), meanvx(nx) )
  ! Linearly interpolate to get x-profile through COM
  j0 = maxloc(y, mask=y < ycom, dim=1)
  k0 = maxloc(z, mask=z < zcom, dim=1)
  j1 = j0 + 1
  k1 = k0 + 1
  ! linear interpolation factors in y and z directions
  ay = (ycom-y(j0))/dy
  az = (zcom-z(k0))/dz
  a00 = (1.0 - ay)*(1.0 - az)
  a10 = ay*(1.0 - az)
  a01 = (1.0 - ay)*az
  a11 = ay*az
  ! find mean ion fraction along x-axis
  meanxi = a00*xi(:,j0,k0) + a10*xi(:,j1,k0) + a01*xi(:,j0,k1) + a11*xi(:,j1,k1)
  ! reverse linear interpolation to find x value for front, assumed to be at xi = 0.5
  xifront = find_zero(0.5 - meanxi, x) ! function assumes that first arg is increasing f(x)

  ! Velocity behind i-front
  ! find mean vx along the same x-line as above
  meanvx = a00*vx(:,j0,k0) + a10*vx(:,j1,k0) + a01*vx(:,j0,k1) + a11*vx(:,j1,k1)
  ! look in a window 4 cells wide, starting from the i-front
  i0 = minloc(x, dim=1, mask=x>xifront)
  ishell = 4
  vxshell = maxval(meanvx(i0:i0+ishell))

  call write_stats()

  

contains
  subroutine write_stats()
    use allstats_vars, only: allstats
    open(11, &
         & file=trim(prefix)//'-'//'globstats'//ctime//'.dat', &
         & action='write')
    write(11, allstats)
    close(11)
  end subroutine write_stats

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

  

end program globstats
