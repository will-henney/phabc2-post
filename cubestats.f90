program cubestats
  ! calculates various statistics for the data cubes
  ! 17 Jul 2008 - minimal modifications to work with the Fabio datasets
  use wfitsutils, only: fitsread, fitscube
  implicit none
  real, parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  real, dimension(:,:,:), allocatable :: d, x, u, v, w, r
  logical, dimension(:,:,:), allocatable :: ion_mask, m
  character(len=128) :: prefix, fitsfilename
  integer :: it1, it2, it, itstep
  integer :: nx, ny, nz, i, j, k
  real :: rmax
  real :: rif_max, rif_min
  real :: vrms_vol_m, vrms_vol_n, vrms_vol_i
  real :: vrms_mass_m, vrms_mass_n, vrms_mass_i
  real :: rx1, rx2, rmean_vol_i, rmean_mass_i
  real :: frac_ion_vol1, frac_ion_vol2, frac_ion_mass
  character(len=1), parameter :: TAB = achar(9)
  character(len=15) :: itstring
  real, parameter :: pi = 3.14159265358979, cubesize = 4.0*3.086e18
  ! WJH 07 Jul 2010 - New distinction between neutral/molecular gas
  ! this is copied from the python code in mhd-pressures.py
  real, parameter :: mol_AV0 = 3.0                           ! position of molecular transition 
  real, parameter :: mol_sharpness = 4.0                     ! sharpness of molecular transition
  real, allocatable, dimension(:,:,:) :: AV, xmol, wi, wn, wm, xm_arg

  print *, 'Run prefix (e.g., 30112005_c)?'
  read '(a)', prefix


  print *, 'First and last time index, and step?'
  read *, it1, it2, itstep

  write(itstring,'(3("-",i4.4))') it1, it2, itstep

  open(1, file=trim(prefix)//itstring//'.stats', action='write')
  open(2, file=trim(prefix)//itstring//'.rstats', action='write')

  write(1, '("# ",10(a,"'//TAB//'"))') 'Time', &
       & 'Ifrac_v', 'Ifrac_v2', 'Ifrac_m', &
       & 'Vrms_vol_m', 'Vrms_vol_n', 'Vrms_vol_i', &
       & 'Vrms_mass_m', 'Vrms_mass_n', 'Vrms_mass_i'
  write(2, '("# ",7(a,"'//TAB//'"))') 'Time', &
       & 'rx1', 'rx2', 'rmean_vol_i', 'rmean_mass_i', 'rif_min', 'rif_max'

  do it = it1, it2, itstep
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-dd', it, '.fits'
     call fitsread(trim(fitsfilename))
     if (it ==it1) then
        ! first time setup
        nx = size(fitscube, 1)
        ny = size(fitscube, 2)
        nz = size(fitscube, 3)
        allocate( d(nx, ny, nz), x(nx, ny, nz), ion_mask(nx, ny, nz) )
        allocate( u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz) )
        allocate( r(nx, ny, nz), m(nx, ny, nz) )
        allocate( AV(nx, ny, nz) )
        allocate( xmol(nx, ny, nz), wi(nx, ny, nz), wn(nx, ny, nz), wm(nx, ny, nz), xm_arg(nx, ny, nz) )
     end if

     d = fitscube/mp/mu

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-xi', it, '.fits'
     call fitsread(trim(fitsfilename))
     x = fitscube

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-vx', it, '.fits'
     call fitsread(trim(fitsfilename))
     u = fitscube
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-vy', it, '.fits'
     call fitsread(trim(fitsfilename))
     v = fitscube
     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-vz', it, '.fits'
     call fitsread(trim(fitsfilename))
     w = fitscube

     write(fitsfilename, '(2a,i4.4,a)') trim(prefix), '-AV', it, '.fits'
     call fitsread(trim(fitsfilename))
     AV = fitscube

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
     ! We cut it out by only considering a sphere of radius 1 pc at start,
     ! growing linearly with time. 
     rmax = (1.0 + real(it)/50.0)*0.25*real(nx)
     m = r < rmax


     ! WJH 05 Jul 2010 - New distinction between neutral/molecular
     ! this is copied from the python code in mhd-pressures.py
     xm_arg = mol_sharpness*(AV-mol_AV0)
     where (xm_arg > 100.0) 
        xmol = 1.0
     elsewhere
        xmol = 1.0 - 1.0/(1.0 + exp(xm_arg))
     end where
     
     
     ! weights for ionized/neutral/molecular
     wi = x
     wn = (1.-x)*(1.-xmol)
     wm = (1.-x)*xmol


     ion_mask = x>0.5

     frac_ion_vol1 = count(ion_mask)/real(nx*ny*nz)
     ! Note that we do not apply the radius cutoff mask to the denominator
     ! since we are considering that the stuff at the edges shoud "really"
     ! be neutral
     frac_ion_vol2 = sum(x, mask=m)/real(nx*ny*nz)
     frac_ion_mass = sum(x*d, mask=m)/sum(d)

     ! direction-averaged radius of ionized volume
     rx1 = cubesize*(frac_ion_vol1*3.0/4.0/pi)**(1./3.)
     rx2 = cubesize*(frac_ion_vol2*3.0/4.0/pi)**(1./3.)

     ! min/max i-front radius
     rif_max = maxval(r, mask=ion_mask)
     rif_min = minval(r, mask=.not.ion_mask)

     ! mean radius of ionized gas
     rmean_vol_i = (cubesize/real(nx))*sum(r*x, mask=m)/sum(x, mask=m)
     rmean_mass_i = (cubesize/real(nx))*sum(r*d*x, mask=m)/sum(d*x, mask=m)

     ! the 1D RMS velocity - this is the one we will use
     vrms_vol_m = sqrt(sum((u*u + v*v + w*w)*wm)/(3*sum(wm)))
     vrms_mass_n = sqrt(sum((u*u + v*v + w*w)*d*wm)/(3*sum(d*wm)))

     vrms_vol_i = sqrt(sum((u*u + v*v + w*w)*wi, mask=m)/(3*sum(wi, mask=m)))
     vrms_mass_i = sqrt(sum((u*u + v*v + w*w)*d*wi, mask=m)/(3*sum(d*wi, mask=m)))

     vrms_vol_n = sqrt(sum((u*u + v*v + w*w)*wn)/(3*sum(wn)))
     vrms_mass_n = sqrt(sum((u*u + v*v + w*w)*d*wn)/(3*sum(d*wn)))

     if (mod(it,10)==0) print *, 'Done timestep: ', it
     write(1, '(i4.4,"'//TAB//'",9(es11.3,"'//TAB//'"))') it, &
          & frac_ion_vol1, frac_ion_vol2, frac_ion_mass, &
          & vrms_vol_m, vrms_vol_n, vrms_vol_i,  &
          & vrms_mass_m, vrms_mass_n, vrms_mass_i


     write(2, '(i4.4,"'//TAB//'",6(es11.3,"'//TAB//'"))') it, &
          & rx1, rx2, rmean_vol_i, rmean_mass_i, rif_min, rif_max
  end do

end program cubestats
