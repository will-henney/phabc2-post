program cubedenstats
  ! calculates mean densities for the data cubes
  use wfitsutils, only: fitsread, fitscube
  use mod_molfrac, only: molfrac
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: boltzmann_k = 1.3806503e-16_dp, mp = 1.67262158e-24_dp, mu = 1.3_dp
  real(dp), dimension(:,:,:), allocatable :: d, x, r
  logical, dimension(:,:,:), allocatable :: ion_mask, m
  real(dp), allocatable, dimension(:,:,:) :: AV, xmol, wi, wn, wm
  character(len=128) :: prefix, fitsfilename
  integer :: it1, it2, it, itstep
  integer :: nx, ny, nz, i, j, k
  real(dp) :: dmean_tot, dmean_m, d2mean_m
  real(dp) :: dmean_n, dmean_i, d2mean_tot, d2mean_n, d2mean_i, d3mean_i
  real(dp) :: rmax
  character(len=1), parameter :: TAB = achar(9)
  character(len=15) :: itstring
  real(dp), parameter :: pi = 3.14159265358979_dp, cubesize = 4.0_dp*3.086e18_dp

  print *, 'Run prefix (e.g., 30112005_c)?'
  read '(a)', prefix

  print *, 'First and last time index, and step?'
  read *, it1, it2, itstep

  write(itstring,'(3("-",i4.4))') it1, it2, itstep

  open(1, file=trim(prefix)//itstring//'.dstats', action='write')
  write(1, '("# ",10(a,"'//TAB//'"))') 'Time', &
       & 'Dmean_tot', 'Dmean_m', 'Dmean_n', 'Dmean_i', &
       & 'D2mean_tot', 'D2mean_m', 'D2mean_n', 'D2mean_i', &
       & 'D3mean_i'

  do it = it1, it2, itstep

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-dd', it, '.fits'
     call fitsread(trim(fitsfilename))
     if (it == it1) then
        ! first time setup
        nx = size(fitscube, 1)
        ny = size(fitscube, 2)
        nz = size(fitscube, 3)
        allocate( d(nx, ny, nz), x(nx, ny, nz), ion_mask(nx, ny, nz) )
        allocate( r(nx, ny, nz), m(nx, ny, nz) )
        allocate( AV(nx, ny, nz), xmol(nx, ny, nz), wi(nx, ny, nz), wn(nx, ny, nz), wm(nx, ny, nz) ) 
     end if
     d = fitscube/mp/mu

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-xi', it, '.fits'
     call fitsread(trim(fitsfilename))
     x = fitscube
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-AV', it, '.fits'
     call fitsread(trim(fitsfilename))
     AV = fitscube

     forall(i=1:nx, j=1:ny, k=1:nz)
        ! radial distance from center in grid units
        r(i, j, k) = &
             & sqrt( &
             &        (real(i, dp) - 0.5_dp*real(nx+1, dp))**2 &
             &      + (real(j, dp) - 0.5_dp*real(ny+1, dp))**2 &
             &      + (real(k, dp) - 0.5_dp*real(nz+1, dp))**2 &
             &     )
     end forall

     ! For early times, we have some partially ionized material
     ! around the edges of the grid, which is skewing our statistics.
     ! We cut it out by only considering a sphere of radius 1 pc. 
     rmax = (1.0_dp + real(it, dp)/50.0_dp)*0.25_dp*real(nx, dp)
     m = r < rmax

     ! weights for ionized/neutral/molecular
     xmol = molfrac(AV)
     wi = max(x, 0.0_dp)
     wn = max((1.0_dp-x)*(1.0_dp-xmol), 0.0_dp)
     wm = max((1.0_dp-x)*xmol, 0.0_dp)

     dmean_tot = sum(d)/real(nx*ny*nz, dp)
     d2mean_tot = sum(d*d)/real(nx*ny*nz, dp)

     d = d / dmean_tot          ! rescale density to prevent overflow

     ! version 2: weight by the ion fraction
     dmean_m = sum(d*wm)/sum(wm)
     dmean_n = sum(d*wn)/sum(wn)
     dmean_i = sum(d*wi, mask=m)/sum(wi, mask=m)

     ! mean of density-squared
     d2mean_m = sum(d*d*wm)/sum(wm)
     d2mean_n = sum(d*d*wn)/sum(wn)
     d2mean_i = sum(d*d*wi, mask=m)/sum(wi, mask=m)

     ! mean of density-cubed - only for ionized gas
     d3mean_i = sum(d*d*d*wi, mask=m)/sum(wi, mask=m)

     ! put scaling back in
     dmean_m = dmean_m * dmean_tot
     dmean_n = dmean_n * dmean_tot
     dmean_i = dmean_i * dmean_tot
     d2mean_m = d2mean_m * dmean_tot * dmean_tot
     d2mean_n = d2mean_n * dmean_tot * dmean_tot
     d2mean_i = d2mean_i * dmean_tot * dmean_tot
     d3mean_i = d3mean_i * dmean_tot * dmean_tot * dmean_tot

     if (mod(it,10)==0) print *, 'Done timestep: ', it

     write(1, '(i4.4,"'//TAB//'",9(es11.3,"'//TAB//'")))') it, &
          & dmean_tot, dmean_m, dmean_n, dmean_i, &
          & d2mean_tot, d2mean_m, d2mean_n, d2mean_i, &
          & d3mean_i

  end do

end program cubedenstats
