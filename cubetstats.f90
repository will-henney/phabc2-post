program cubetstats
  ! calculates temperature fluctuations: tsquared
  ! 28 Jul 2008 - minimal modifications to work with the Fabio datasets
  use wfitsutils, only: fitsread, fitscube
  implicit none
  real, dimension(:,:,:), allocatable :: d, x, r, t, edeniden, p
  logical, dimension(:,:,:), allocatable :: hi_mask, m
  character(len=128) :: prefix, fitsfilename
  integer :: it1, it2, it, itstep
  integer :: nx, ny, nz, i, j, k
  real :: tmean, tsq, edeniden_sum
  real :: edeniden_sum_hi, edeniden_sum_lo, tmean_hi, tmean_lo, tsq_hi, tsq_lo
  real :: rmax
  real :: x0 = 0.99875 ! border between N+ and O++ zones
  character(len=1), parameter :: TAB = achar(9)
  character(len=15) :: itstring
  real, parameter :: pi = 3.14159265358979, cubesize = 4.0*3.086e18
  real, parameter :: mu = 1.3, mp = 1.67262158e-24, boltzmann_k = 1.3806503e-16
  real :: tmean_di
  
  print *, 'Run prefix (e.g., 30112005_c)?'
  read '(a)', prefix

  print *, 'First and last time index, and step?'
  read *, it1, it2, itstep

  write(itstring,'(3("-",i4.4))') it1, it2, itstep

  open(1, file=trim(prefix)//itstring//'.tstats', action='write')
  write(1, '("# ",9(2a))') 'Time', TAB, 'Tmean', TAB, 'tsq', TAB, &
       & 'Tmean(O++)', TAB, 'tsq(O++)', TAB, &
       & 'Tmean(N+)', TAB, 'tsq(N+)', TAB, 'Tdyna', TAB, 'VEM(O++)/VEM(H+)', TAB

  do it = it1, it2, itstep

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-dd', it, '.fits'
     call fitsread(trim(fitsfilename))
     if (it == it1) then
        ! first time setup
        nx = size(fitscube, 1)
        ny = size(fitscube, 2)
        nz = size(fitscube, 3)
        allocate( d(nx, ny, nz), x(nx, ny, nz), t(nx, ny, nz))
        allocate( edeniden(nx, ny, nz), hi_mask(nx, ny, nz) )
        allocate( r(nx, ny, nz), m(nx, ny, nz) )
     end if
     d = fitscube

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-xi', it, '.fits'
     call fitsread(trim(fitsfilename))
     x = fitscube
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-te', it, '.fits'
     call fitsread(trim(fitsfilename))
     t = fitscube

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-pp', it, '.fits'
     call fitsread(trim(fitsfilename))
     p = fitscube

     forall(i=1:nx, j=1:ny, k=1:nz)
        ! radial distance from center in grid units
        r(i, j, k) = &
             & sqrt( &
             &        (real(i) - 0.5*real(nx+1))**2 &
             &      + (real(j) - 0.5*real(ny+1))**2 &
             &      + (real(k) - 0.5*real(nz+1))**2 &
             &     )
     end forall

     ! For early times, we have some partially ionized material
     ! around the edges of the grid, which is skewing our statistics.
     ! We cut it out by only considering a sphere of radius 1 pc. 
     rmax = (1.0 + real(it)/50.0)*0.25*real(nx)
     m = r < rmax

     ! divide ionized gas into high and low ionization
     hi_mask = x > x0

     ! weight by the ion density x electron density
     edeniden = (d*x)**2
     edeniden_sum = sum(edeniden, mask=m)

     ! mean ionized temperature
     tmean = sum(t*edeniden, mask=m)/edeniden_sum
     ! ionized t-squared fluctuation
     tsq = sum((t-tmean)**2 * edeniden, mask=m)/(tmean**2 * edeniden_sum)

     ! repeat for highly ionized
     edeniden_sum_hi = sum(edeniden, mask=m.and.hi_mask)
     tmean_hi = sum(t*edeniden, mask=m.and.hi_mask)/edeniden_sum_hi
     tsq_hi = sum((t-tmean_hi)**2 * edeniden, mask=m.and.hi_mask)&
          & /(tmean_hi**2 * edeniden_sum_hi)
     ! repeat for lowly ionized
     edeniden_sum_lo = sum(edeniden, mask=m.and.(.not.hi_mask))
     tmean_lo = sum(t*edeniden, mask=m.and.(.not.hi_mask))/edeniden_sum_lo
     tsq_lo = sum((t-tmean_lo)**2 * edeniden, mask=m.and.(.not.hi_mask))&
          & /(tmean_lo**2 * edeniden_sum_lo)
     
     ! a different mean temperature, suitable for the dynamics
     tmean_di = (mp * mu / boltzmann_k )*sum(p*x/ (1.+x) )/sum(d*x)
     
     if (mod(it,10)==0) print *, 'Done timestep: ', it

     write(1, '(i4.4,a,7(es9.3,a))') it, TAB, tmean, TAB, tsq, TAB, &
          & tmean_hi, TAB, tsq_hi, TAB, tmean_lo, TAB, tsq_lo, TAB, &
          & tmean_di, TAB, &
          & edeniden_sum_hi/edeniden_sum, TAB

  end do

end program cubetstats
