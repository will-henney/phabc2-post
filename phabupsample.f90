program phabupsample
  use wfitsutils, only: fitswrite, fitsread, fitscube, nx_fits
  implicit none
  integer, parameter :: neq = 10
  integer :: n, nn
  integer :: i, j, k, ieq, it = 0, idebug = 1
  character(len=2), dimension(neq) :: varid = (/ 'pp', 'vx', 'vy', 'vz', 'bx', 'by', 'bz', 'dd', 'xn', 'xi' /)
  character(len=256) :: runid, newrunid, outfile
  
  real, allocatable, dimension(:,:,:,:) :: newprim, newu, prim, u
  real, parameter :: gamma = 1.66666666667

  ! Get input parameters
  print *, 'This program rebins a sim (in current dir) to double resolution'
  print *, 'Original Run ID?'
  read '(a)', runid
  print *, 'Save time?'
  read *, it
  print *, 'New Run ID?'
  read '(a)', newrunid

  ! Read in prims from hi-res FITS file
  do ieq = 1, neq
     write(outfile,'(3a,i4.4,a)')  &
          & trim(runid), '-', varid(ieq), it, '.fits' 
     call fitsread(trim(outfile))

     ! On first trip, must allocate arrays before anything else
     if (ieq == 1) then
        n = nx_fits
        nn = 2*n
        allocate( newprim(neq,nn,nn,nn), newu(neq,nn,nn,nn), prim(neq,n,n,n), u(neq,n,n,n) )
        if (idebug > 0) print *, 'Allocated arrays: converting ', n, ' -> ', nn
     end if
     
     prim(ieq,:,:,:) = fitscube
     if (idebug > 0) print *, varid(ieq), ' read in'
  end do
  
  ! Convert to conserved variables
  forall ( i = 1:n, j = 1:n, k = 1:n )
     u(:,i,j,k) = primu(prim(:,i,j,k))
  end forall
  
  if (idebug > 0) print *, 'primu done'
  
  ! rebin
  forall( ieq = 1:neq, i = 1:nn, j = 1:nn, k = 1:nn )
     newu(ieq, i, j, k) = sum( u(ieq, 2*i-1:2*i, 2*j-1:2*j, 2*k-1:2*k) ) / 8.0
  end forall
  if (idebug > 0) print *, 'rebin done'

  ! Convert back to primitive variables
  forall ( i = 1:nn, j = 1:nn, k = 1:nn )
     newprim(:,i,j,k) = uprim(newu(:,i,j,k))
  end forall
  if (idebug > 0) print *, 'uprim done'

  ! Do the pressure specially by averaging the primitive
  forall( i = 1:nn, j = 1:nn, k = 1:nn )
     newprim(1, i, j, k) = sum( prim(1, 2*i-1:2*i, 2*j-1:2*j, 2*k-1:2*k) ) / 8.0
  end forall
  

  ! Write out low-res FITS files
  do ieq = 1, neq
     write(outfile,'(3a,i4.4,a)')  &
          & trim(newrunid), '-', varid(ieq), it, '.fits' 
     call fitswrite(newprim(ieq,:,:,:), trim(outfile))
     if (idebug > 0) print *, varid(ieq), ' written out'
  end do

  contains
    pure function primu(p) result(uu)
      real, intent(in) :: p(:)
      real :: uu(size(p))
      real :: btot, vv
      uu(5:8)=p(5:8)
      uu(2:4)=p(2:4)*p(8)
      btot=p(5)**2+p(6)**2+p(7)**2
      vv=p(2)**2+p(3)**2+p(4)**2
      uu(1)=0.5*p(8)*vv+p(1)/(gamma-1.0)+btot*0.5
      uu(9:10)=p(9:10)*p(8)       ! scalars
    end function primu
    
    
    pure function uprim(u) result(p)
      real, intent(in) :: u(:)
      real :: p(size(u))
      real :: btot, vv
      p(9:10)=u(9:10)/u(8)        ! scalars
      p(5:8)=u(5:8)
      p(2:4)=u(2:4)/p(8)
      btot=p(5)**2+p(6)**2+p(7)**2
      vv=p(2)**2+p(3)**2+p(4)**2
      p(1)=(u(1)-0.5*p(8)*vv-0.5*btot)*(gamma-1.0)
    end function uprim
    

end program phabupsample
