program enro2fits
  ! Converts Enrique data files to FITS format for use as initial conditions
  use wfitsutils, only: fitswrite
  implicit none
  integer, parameter :: dp = kind(1.0d0), sp = kind(1.0)
  integer :: enronumber, n
  integer, parameter :: ne = 7
  character(len=2) :: enrovarid(ne) = (/'ro', 'vx', 'vy', 'vz', 'bx', 'by', 'bz'/)
  real(sp), allocatable :: enrovars(:,:,:,:) 
  character(len=10) :: enrofile
  real(sp), allocatable, dimension(:,:,:) :: pp, vx, vy, vz, bx, by, bz, dd, xi, xn
  ! Parameters of turbulence simulation
  real(kind=dp),parameter :: Jeans=4.0_dp ! Jeans number
  real(kind=dp),parameter :: csound=0.2e5_dp ! sound speed
  real(kind=dp),parameter :: mu_enro=2.4_dp ! mean weight per particle
  real(kind=dp),parameter :: mp=1.67262158e-24_dp ! proton mass
  real(kind=dp),parameter :: boxsize=4.0_dp ! parsecs
  real(kind=dp),parameter :: gamma=5.0_dp/3.0_dp
  real(kind=dp),parameter :: beta = 0.1_dp
  real(kind=sp),parameter :: pi = 3.1415926535897932385e0
  real(dp) :: dscale, bscale, cellsize
  character(len=32) :: outid
  integer, dimension(3) :: center_here
  integer :: idim, shft, i, j, k
  real(sp) :: cutrad
  logical, allocatable, dimension(:,:,:) :: cutmask
  real(sp) :: ek_mean, eb_mean, p_mean, bx_mean, by_mean, bz_mean, b_rms
  real(sp) :: beta_ordered, beta_total

  print *, 'Run number?'
  read *, enronumber
  print *, 'Gridsize?'
  read *, n
  cellsize = boxsize*3.085677582e18/n

  ! read in all the enrique vars into one big array
  allocate( enrovars(n,n,n,7) )
  allocate( pp(n,n,n), vx(n,n,n), vy(n,n,n), vz(n,n,n), bx(n,n,n), &
       &    by(n,n,n), bz(n,n,n), dd(n,n,n), xi(n,n,n), xn(n,n,n), &
       &    cutmask(n,n,n))
  do i = 1, ne
     write(enrofile,'(a,i4.4,a)') enrovarid(i), enronumber, '.dat'
     open(1, file=enrofile, form='unformatted', action='read')
     read(1) enrovars(:,:,:,i)
  end do
  
  ! Take care of centering
  center_here=maxloc(enrovars(:,:,:,1))
  do idim=1,3
     shft=center_here(idim)-n/2
     enrovars=cshift(enrovars,shft,idim)
  enddo
    
  write(*,*) 'Maximum density: ',maxval(enrovars(:,:,:,1)), &
       ' at ',maxloc(enrovars(:,:,:,1))
    
  ! Subtract velocity of source point
  enrovars(:, :, :, 2) = enrovars(:, :, :, 2) - enrovars(n/2, n/2, n/2, 2)
  enrovars(:, :, :, 3) = enrovars(:, :, :, 3) - enrovars(n/2, n/2, n/2, 3)
  enrovars(:, :, :, 4) = enrovars(:, :, :, 4) - enrovars(n/2, n/2, n/2, 4)

  dscale = 500.0_dp*((Jeans/boxsize)**2)*mu_enro*mp
  dd = dscale*enrovars(:,:,:,1)
  vx = csound*enrovars(:,:,:,2)
  vy = csound*enrovars(:,:,:,3)
  vz = csound*enrovars(:,:,:,4)

  ! Scaling of the magnetic field
  ! This should put it in CGS units [Gauss]
  bscale = sqrt(4.0d0*pi*dscale)*csound
  print *, 'bscale is ', bscale
  bx = bscale*enrovars(:,:,:,5)
  print *, 'Mean Bx in code units: ', sum(real(enrovars(:,:,:,5), dp))/real(n*n*n, dp)
  by = bscale*enrovars(:,:,:,6)
  bz = bscale*enrovars(:,:,:,7)

  ! csound is isothermal sound speed
  pp = dd*(csound**2) ! /gamma
  ! start with neutral conditions
  xi = 1.e-6
  xn = 1.0 - xi

  ! Work out some mean properties...

  ! Kinetic energy
  ek_mean = 0.5*sum(dd*(vx**2 + vy**2 + vz**2))/real(n**3)
  ! Pressure
  p_mean = sum(pp)/real(n**3)
  ! Vector mean of field gives ordered component
  bx_mean = sum(bx)/real(n**3)
  by_mean = sum(by)/real(n**3)
  bz_mean = sum(bz)/real(n**3)
  beta_ordered = 8.0*pi*p_mean/(bx_mean**2 + by_mean**2 + bz_mean**2)
  ! RMS mean of field gives total of ordered + disordered components
  b_rms = sqrt(sum(bx**2 + by**2 + bz**2)/real(n**3))
  eb_mean = b_rms**2/(8.0*pi)
  beta_total = p_mean/eb_mean

  print '(a,"[",f0.3,",",f0.3,",",f0.3,"]")', &
       & 'Mean B field (micro Gauss): ', &
       & 1.e6*(/bx_mean, by_mean, bz_mean/)
  print '(a,f0.3)', 'RMS B field (micro Gauss): ', 1.e6*b_rms
  print '(a,es10.4,"/",es10.4,"/",es10.4)', &
       & 'Mean energies (Thermal/Kinetic/Magnetic): ', &
       & p_mean, ek_mean, eb_mean
  print '(a,f0.5,"/",f0.5)', 'Mean beta (ordered/total): ', beta_ordered, beta_total

  ! Remove mass from the central cells
  print *, 'Radius of region to remove mass from (pixels)?' 
  read *, cutrad
  forall(i=1:n, j=1:n, k=1:n)
     ! use a sharp-bounded spherical volume
     cutmask(i, j, k) = sum(((/i, j, k/) - n/2)**2) <= cutrad**2
  end forall
  print *, 'Mass removed is ',&
       &  0.9*sum(dd, mask=cutmask)*cellsize**3 / 1.989e33, ' solar masses'
  where(cutmask)
     dd = 0.1*dd
     pp = 0.1*pp
  end where

  print *, 'Run id?'
  read '(a)', outid

  call fitswrite(dd, trim(outid)//'-dd0000.fits')
  call fitswrite(vx, trim(outid)//'-vx0000.fits')
  call fitswrite(vy, trim(outid)//'-vy0000.fits')
  call fitswrite(vz, trim(outid)//'-vz0000.fits')
  ! put B field in Fabio code units
  call fitswrite(bx/sqrt(4.0*pi), trim(outid)//'-bx0000.fits')
  call fitswrite(by/sqrt(4.0*pi), trim(outid)//'-by0000.fits')
  call fitswrite(bz/sqrt(4.0*pi), trim(outid)//'-bz0000.fits')
  call fitswrite(pp, trim(outid)//'-pp0000.fits')
  call fitswrite(xi, trim(outid)//'-xi0000.fits')
  call fitswrite(xn, trim(outid)//'-xn0000.fits')

end program enro2fits
