program cubeextras
  ! calculates cubes of various derived quantities
  ! cubes of density, ionization fraction, pressure, b field, velocity
  use emissmod, only: emtype, emissivity, dp, pah_emissivity
  implicit none
  real, dimension(:,:,:), allocatable :: dd, xi, te, av
  real(dp), dimension(:,:,:), allocatable :: em
  character(len=128) :: prefix
  real, parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  integer :: nx = 0, ny = 0, nz = 0
  integer :: itime
  character(len=4) :: ctime
  integer, parameter :: idebug = 1

  print *, 'File prefix?'
  read '(a)', prefix

  print *, 'Save time?'
  read *, itime

  write(ctime, '(i4.4)') itime

  call read_var(dd, 'dd')
  dd = dd / (mu*mp)
  call read_var(xi, 'xi')
  if (idebug >= 2) then
     print *, 'Density max, min: ', maxval(dd), minval(dd)
     print *, 'Ion frac max, min: ', maxval(xi), minval(xi)
  end if
  
  ! emissivity
  print *, 'Emissivity type?'
  read '(a)', emtype

  allocate( em(nx, ny, nz) )
  if (emtype(1:3)=='PAH') then
     ! PAH emission depends on UV flux, and hence on A_V array
     call read_var(AV, 'AV')
     em = real(pah_emissivity(dd, xi, av), dp)
  else
     ! others depend on temperature
     call read_var(te, 'te')
     if (idebug >= 2) then
        print *, 'Temperature max, min: ', maxval(te), minval(te), minloc(te)
     endif
     em = emissivity(dd, xi, te)
  end if
  
  if (idebug >= 1) then
     print '(3a)', 'Emissivity stats for ', emtype, ':'
     print '(a,es12.3,a,3(i0,tr1))', 'Max value: ', maxval(em), ' at ', maxloc(em)
     print '(a,es12.3,a,3(i0,tr1))', 'Min value: ', minval(em), ' at ', minloc(em)
     print '(a,es12.3)', 'Mean value: ', sum(em)/product(shape(em))
  endif


  call write_var(real(em), 'e-'//emtype)

contains

  subroutine read_var(var, id)
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
  end subroutine read_var

  subroutine write_var(var, id)
    use wfitsutils, only: fitswrite
    real, intent(in), dimension(:,:,:) :: var
    character(len=*), intent(in) :: id
    call fitswrite(var, trim(prefix)//'-'//id//ctime//'.fits')
  end subroutine write_var


end program cubeextras
