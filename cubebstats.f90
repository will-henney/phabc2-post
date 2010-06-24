program cubevstats
  ! calculates mean densities for the data cubes
  use wfitsutils, only: fitsread, fitscube
  implicit none
  real, parameter :: boltzmann_k = 1.3806503e-16, mp = 1.67262158e-24, mu = 1.3
  real, dimension(:,:,:), allocatable :: d, xi, bx, by, bz, bb
  logical, dimension(:,:,:), allocatable :: m
  print *, 'Run prefix (e.g., 30112005_c)?'
  read '(a)', prefix

  print *, 'First and last time index, and step?'
  read *, it1, it2, itstep

  write(itstring,'(3("-",i4.4))') it1, it2, itstep

  open(1, file=trim(prefix)//itstring//'.vstats', action='write')

  do it = it1, it2, itstep
     if (it==1) then 
