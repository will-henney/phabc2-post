! CALCULATES VARIOUS STATISTICS FOR THE DATA CUBES
!
! Usage: printf '%s\n' ID T1 T2 TSTEP| ./cubestats
!
! 17 Jul 2008 - minimal modifications to work with the Fabio datasets

program cubestats
  use wfitsutils, only: fitsread, fitscube
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: boltzmann_k = 1.3806503e-16_dp, mp = 1.67262158e-24_dp, mu = 1.3_dp
  real(dp), dimension(:,:,:), allocatable :: d, x, u, v, w, r
  logical, dimension(:,:,:), allocatable :: ion_mask, m
  character(len=128) :: prefix, fitsfilename
  integer :: it1, it2, it, itstep
  integer :: nx, ny, nz, i, j, k
  real(dp) :: rmax
  real(dp) :: rif_max, rif_min
  real(dp) :: vrms_vol_m, vrms_vol_n, vrms_vol_i
  real(dp) :: vrms_mass_m, vrms_mass_n, vrms_mass_i
  real(dp) :: rx1, rx2, rmean_vol_i, rmean_mass_i, rpdr
  real(dp) :: frac_ion_vol1, frac_ion_vol2, frac_ion_mass
  real(dp) :: frac_mol_vol, frac_mol_mass, frac_neut_vol, frac_neut_mass
  character(len=1), parameter :: TAB = achar(9)
  character(len=15) :: itstring
  real(dp), parameter :: pi = 3.14159265358979_dp, cubesize = 4.0_dp*3.086e18_dp
  ! WJH 07 Jul 2010 - New distinction between neutral/molecular gas
  ! this is copied from the python code in mhd-pressures.py
  real(dp), parameter :: mol_AV0 = 3.0_dp                           ! position of molecular transition 
  real(dp), parameter :: mol_sharpness = 4.0_dp                     ! sharpness of molecular transition
  real(dp), allocatable, dimension(:,:,:) :: AV, xmol, wi, wn, wm, xm_arg
  real(dp), allocatable, dimension(:,:,:) :: sumfrac_vol, sumfrac_mass ! sum of weights

  print *, 'Run prefix (e.g., 30112005_c)?'
  read '(a)', prefix


  print *, 'First and last time index, and step?'
  read *, it1, it2, itstep

  write(itstring,'(3("-",i4.4))') it1, it2, itstep

  open(1, file=trim(prefix)//itstring//'.stats', action='write')
  open(2, file=trim(prefix)//itstring//'.rstats', action='write')

  write(1, '("# ",13(a,"'//TAB//'"))') 'Time', &
       & 'Mfrac_v', 'Mfrac_m', &
       & 'Nfrac_v', 'Nfrac_m', &
       & 'Ifrac_v2', 'Ifrac_m', &
       & 'Vrms_vol_m', 'Vrms_vol_n', 'Vrms_vol_i', &
       & 'Vrms_mass_m', 'Vrms_mass_n', 'Vrms_mass_i'
  write(2, '("# ",8(a,"'//TAB//'"))') 'Time', &
       & 'rx1', 'rx2', 'rmean_vol_i', 'rmean_mass_i', 'rif_min', 'rif_max', 'rpdr'

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
        allocate( sumfrac_vol(nx, ny, nz), sumfrac_mass(nx, ny, nz) )
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
             &        (real(i, dp) - 0.5_dp*real(nx+1, dp))**2 &
             &      + (real(j, dp) - 0.5_dp*real(ny+1, dp))**2 &
             &      + (real(k, dp) - 0.5_dp*real(nz+1, dp))**2 &
             &     )
     end forall
     ! For early times, we have some partially ionized material
     ! around the edges of the grid, which is skewing our statistics.
     ! We cut it out by only considering a sphere of radius 1 pc at start,
     ! growing linearly with time. 
     rmax = (1.0_dp + real(it, dp)/50.0_dp)*0.25_dp*real(nx, dp)
     m = r < rmax


     ! WJH 05 Jul 2010 - New distinction between neutral/molecular
     ! this is copied from the python code in mhd-pressures.py
     xm_arg = mol_sharpness*(AV-mol_AV0)
     where (xm_arg > 10.0_dp) 
        xmol = 1.0_dp
     elsewhere
        xmol = 1.0_dp - 1.0_dp/(1.0_dp + exp(xm_arg))
     end where
     
     
     ! weights for ionized/neutral/molecular
     wi = max(x, 0.0_dp)
     wn = max((1._dp-x)*(1._dp-xmol), 0.0_dp)
     wm = max((1._dp-x)*xmol, 0.0_dp)


     ion_mask = x>0.5_dp

     frac_ion_vol1 = count(ion_mask)/real(nx*ny*nz, dp)
     ! Note that we do not apply the radius cutoff mask to the denominator
     ! since we are considering that the stuff at the edges shoud "really"
     ! be neutral
     frac_ion_vol2 = sum(x, mask=m)/real(nx*ny*nz, dp)
     frac_ion_mass = sum(x*d, mask=m)/sum(d)

     frac_mol_vol = sum(wm)/real(nx*ny*nz, dp)
     frac_mol_mass = sum(wm*d)/sum(d)

     frac_neut_vol = sum(wn)/real(nx*ny*nz, dp)
     frac_neut_mass = sum(wn*d)/sum(d)

     ! direction-averaged radius of ionized volume
     rx1 = cubesize*(frac_ion_vol1*3.0_dp/4.0_dp/pi)**(1._dp/3._dp)
     rx2 = cubesize*(frac_ion_vol2*3.0_dp/4.0_dp/pi)**(1._dp/3._dp)

     ! same for dissociation front
     rpdr = cubesize*((frac_ion_vol2 + frac_neut_vol)*3.0_dp/4.0_dp/pi)**(1._dp/3._dp)

     ! min/max i-front radius
     rif_max = maxval(r, mask=ion_mask)
     rif_min = minval(r, mask=.not.ion_mask)

     ! mean radius of ionized gas
     rmean_vol_i = (cubesize/real(nx, dp))*sum(r*x, mask=m)/sum(x, mask=m)
     rmean_mass_i = (cubesize/real(nx, dp))*sum(r*d*x, mask=m)/sum(d*x, mask=m)

     ! the 1D RMS velocity - this is the one we will use
     vrms_vol_m = sqrt(sum((u*u + v*v + w*w)*wm)/(3*sum(wm)))
     vrms_mass_m = sqrt(sum((u*u + v*v + w*w)*d*wm)/(3*sum(d*wm)))

     vrms_vol_i = sqrt(sum((u*u + v*v + w*w)*wi, mask=m)/(3*sum(wi, mask=m)))
     vrms_mass_i = sqrt(sum((u*u + v*v + w*w)*d*wi, mask=m)/(3*sum(d*wi, mask=m)))

     vrms_vol_n = sqrt(sum((u*u + v*v + w*w)*wn)/(3*sum(wn)))
     vrms_mass_n = sqrt(sum((u*u + v*v + w*w)*d*wn)/(3*sum(d*wn)))

     if (mod(it,10)==0) print *, 'Done timestep: ', it
     write(1, '(i4.4,"'//TAB//'",12(es11.3,"'//TAB//'"))') it, &
          & frac_mol_vol, frac_mol_mass, &
          & frac_neut_vol, frac_neut_mass, &
          & frac_ion_vol2, frac_ion_mass, &
          & vrms_vol_m, vrms_vol_n, vrms_vol_i,  &
          & vrms_mass_m, vrms_mass_n, vrms_mass_i  


     write(2, '(i4.4,"'//TAB//'",7(es11.3,"'//TAB//'"))') it, &
          & rx1, rx2, rmean_vol_i, rmean_mass_i, rif_min, rif_max, rpdr

     
     ! checking that all the fractions add up to unity globally
     print *, 'Sum of global volume fractions = ', frac_ion_vol2 + frac_neut_vol + frac_mol_vol 
     print *, 'Sum of global mass fractions = ', frac_ion_mass + frac_neut_mass + frac_mol_mass 

     ! check the same locally
     sumfrac_vol = wi + wn + wm 
     print *, 'min/max summed volume fracs = ', minval(sumfrac_vol), maxval(sumfrac_vol)

  end do

end program cubestats
