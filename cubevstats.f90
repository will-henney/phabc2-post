program cubevstats
  ! calculates mean densities for the data cubes
  use wfitsutils, only: fitsread, fitscube
  implicit none
  real, parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  real, dimension(:,:,:), allocatable :: d, xi, vr, r, x, y, z, vx, vy, vz
  logical, dimension(:,:,:), allocatable :: m
  character(len=128) :: prefix, fitsfilename
  integer :: it1, it2, it, itstep
  integer :: nx, ny, nz, i, j, k
  real :: vr_vol_tot, vr_vol_n, vr_vol_i
  real :: vr_mass_tot, vr_mass_n, vr_mass_i
  real :: vr_em_tot, vr_em_n, vr_em_i, vr_em_if
  character(len=1), parameter :: TAB = achar(9)
  character(len=15) :: itstring
  real, parameter :: pi = 3.14159265358979, cubesize = 4.0*3.086e18
!   integer, parameter :: itsmall = 70
  real :: rmax

  print *, 'Run prefix (e.g., 30112005_c)?'
  read '(a)', prefix

  print *, 'First and last time index, and step?'
  read *, it1, it2, itstep

  write(itstring,'(3("-",i4.4))') it1, it2, itstep

  open(1, file=trim(prefix)//itstring//'.vstats', action='write')

  do it = it1, it2, itstep
     if (it==1) then 
        write(1, '("# ",11(a,"'//TAB//'"))') 'Time', &
             & 'Vr_vol_t', 'Vr_vol_n', 'Vr_vol_i', &
             & 'Vr_mass_t', 'Vr_mass_n', 'Vr_mass_i', &
             & 'Vr_em_t', 'Vr_em_n', 'Vr_em_i', 'Vr_em_if'

     end if

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-dd', it, '.fits'
     call fitsread(trim(fitsfilename))
     if (it ==it1) then
        ! first time setup
        nx = size(fitscube, 1)
        ny = size(fitscube, 2)
        nz = size(fitscube, 3)
        allocate( d(nx, ny, nz), xi(nx, ny, nz), vr(nx, ny, nz) )
        allocate( x(nx, ny, nz), y(nx, ny, nz), z(nx, ny, nz), r(nx, ny, nz) )
        allocate( vx(nx, ny, nz), vy(nx, ny, nz), vz(nx, ny, nz) )
        allocate( m(nx, ny, nz) )
     end if
     d = fitscube/mp/mu

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-xi', it, '.fits'
     call fitsread(trim(fitsfilename))
     xi = fitscube

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-vx', it, '.fits'
     call fitsread(trim(fitsfilename))
     vx = fitscube
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-vy', it, '.fits'
     call fitsread(trim(fitsfilename))
     vy = fitscube
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-vz', it, '.fits'
     call fitsread(trim(fitsfilename))
     vz = fitscube

     ! Position
     forall(i=1:nx, j=1:ny, k=1:nz)
        ! cartesian distances from the center
        x(i,j,k) = real(i) - 0.5*real(nx+1)
        y(i,j,k) = real(j) - 0.5*real(ny+1)
        z(i,j,k) = real(k) - 0.5*real(nz+1)
     end forall
     r = sqrt(x**2 + y**2 + z**2)
     ! radial component of velocity
     vr = (vx*x + vy*y + vz*z)/r

     ! For early times, we have some partially ionized material
     ! around the edges of the grid, which is skewing our statistics.
     ! We cut it out by only considering a sphere of radius 1 pc. 
     rmax = (1.0 + real(it)/50.0)*0.25*real(nx)
     m = r < rmax
     
     ! volume-weighted averages
     vr_vol_tot = sum(vr)/real(nx*ny*nz)
     vr_vol_n = sum(vr*(1.-xi))/sum(1.-xi)
     vr_vol_i = sum(vr*xi, mask=m)/sum(xi, mask=m)

     ! mass-weighted averages
     vr_mass_tot = sum(vr*d)/sum(d)
     vr_mass_n = sum(vr*(1.-xi)*d)/sum((1.-xi)*d)
     vr_mass_i = sum(vr*xi*d, mask=m)/sum(xi*d, mask=m)

     ! emission-weighted averages (assumed prop to n^2)
     vr_em_tot = sum(vr*d**2)/sum(d**2)
     vr_em_n = sum(vr*((1.-xi)*d)**2)/sum(((1.-xi)*d)**2)
     vr_em_i = sum(vr*(xi*d)**2, mask=m)/sum((xi*d)**2, mask=m)
     ! and finally, with weighting for partially ionized gas at i-front
     vr_em_if = sum(vr*xi*(1.-xi)*d**2, mask=m)/sum(xi*(1.-xi)*d**2, mask=m)

     write(1, '(i4.4,"'//TAB//'",10(es11.3,"'//TAB//'")))') it, &
             & vr_vol_tot, vr_vol_n, vr_vol_i, &
             & vr_mass_tot, vr_mass_n, vr_mass_i, &
             & vr_em_tot, vr_em_n, vr_em_i, vr_em_if

  end do
  
end program cubevstats
