program omp_test
  ! Testbed for OpenMP stuff
  !
  ! This is based on cubeextras
  implicit none
  real, dimension(:,:,:), allocatable :: dd, xi, pp
  real, dimension(:,:,:), allocatable, save :: te
  real :: tesum
  character(len=128) :: prefix
  real, parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  integer :: nx = 0, ny = 0, nz = 0
  integer :: itime, i
  character(len=4) :: ctime
  real :: time_end, time_start
  integer:: icount_end, icount_start, icount_rate
  integer :: ntrips

  !$OMP THREADPRIVATE(te)

  print *, 'Number of trips?'
  read *, ntrips

  ! Fix selection of simulation
  prefix = 'Bstar-ep'
  itime = 134
  write(ctime, '(i4.4)') itime

  call read_var(dd, 'dd')
  call read_var(xi, 'xi')
  call read_var(pp, 'pp')

  ! Temperature
  !$OMP PARALLEL
  allocate( te(nx, ny, nz) )
  !$OMP END PARALLEL

  call system_clock(icount_start, icount_rate)
  call cpu_time(time_start)


  tesum = 0.0
  !$OMP PARALLEL
  !$OMP DO
  do i = 1, ntrips
     te = real(i) * mp * mu * pp / dd / (1.+xi) / boltzmann_k 
     !$OMP ATOMIC
     tesum = tesum + sum(te)
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  call cpu_time(time_end)
  call system_clock(icount_end, icount_rate)

  print *, 'Grand total sum: ', tesum

  print *, 'Number of trips: ', ntrips
  print *, 'Time to calculate each cube (seconds): '
  print *, 'CPU time = ', (time_end - time_start)/real(ntrips)
  print *, 'Wall time = ', real(icount_end - icount_start)/icount_rate/real(ntrips)

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


end program omp_test
