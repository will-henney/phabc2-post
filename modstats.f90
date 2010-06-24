module modstats
  ! WJH 11 Sep 2004
  ! Various helper routines to find statistics
  implicit none
  integer, private, parameter :: DP = kind(1.0d0)
  type pair
    real(DP) :: a, b
  end type pair
  type pair_s
    real :: a, b
  end type pair_s

  interface find_peak
    module procedure find_peak_d, find_peak_s
  end interface
  interface find_halfmax
    module procedure find_halfmax_d, find_halfmax_s
  end interface
  interface find_zero
    module procedure find_zero_d, find_zero_s
  end interface
  
contains

  function find_zero_d( y, x ) result (x0)
    real(DP) :: x0
    real(DP), intent(in), dimension(:) :: y, x
    !! function that semi-accurately finds the root y(x) = 0 for
    !! discrete arrays x_i, y_i by linear interpolation
    !! NOTE: assumes y(x) increasing with x and only one zero-crossing !!!!
    integer :: i1
    i1 = maxloc(x,dim=1,mask=y<0.0_dp)
    if (i1>=size(x)) then
       x0 = x(size(x))
    else
       x0 = x(i1) - y(i1)*(x(i1+1)-x(i1))/(y(i1+1)-y(i1))
    end if
  end function find_zero_d
  
  function find_peak_d( y, x ) result (xpeak)
    real(DP) :: xpeak
    real(DP), intent(in), dimension(:) :: y, x
    !! function that accurately finds the peak of an array
    !! by fitting a parabola through the largest element plus
    !! its two neighbours
    integer :: i1, i2, i3
    real(DP) :: x1, x2, x3, y1, y2, y3, alpha, beta, gamma

    ! peak and points to either side of peak
    i2 = maxloc(y,dim=1)
    if (i2==size(y)) then
       xpeak = x(size(x))
       return
    end if
    if (i2==1) then
       xpeak = x(1)
       return
    end if
    i1 = i2-1
    i3 = i2+1
    ! corresponding velocities
    x1 = x(i1)
    x2 = x(i2)
    x3 = x(i3)
    ! corresponding intensities
    y1 = y(i1)
    y2 = y(i2)
    y3 = y(i3)
    if (x1/=x2 .and. x2/=x3) then
       alpha = y1/(x1-x2)/(x1-x3)
       beta =  y2/(x2-x3)/(x2-x1)
       gamma = y3/(x3-x2)/(x3-x1)
       xpeak = 0.5 * ( x1*(beta+gamma) + x2*(alpha+gamma) + x3*(alpha+beta) ) / (alpha+beta+gamma)
    else
       xpeak = x2
    end if
    
  end function find_peak_d


  function find_halfmax_d ( y, x, fraction ) result (halfmax)
    type(pair) :: halfmax
    real(DP), intent(in), dimension(:) :: y, x
    real(DP), optional, intent(in) :: fraction
    integer :: ipeak, ileft, iright
    real(DP) :: ypeak, xleft, xright, frac

    if (present(fraction)) then
       frac = fraction
    else
       ! unless otherwise indicated, find width at HALF maximum
       frac = 0.5
    end if

    ipeak = maxloc(y,dim=1)
    ypeak = y(ipeak)
    ! find index of leftmost and rightmost pixel with flux > 1/2 max
    ileft = minloc(x, dim=1, mask=y>=frac*ypeak)
    iright = maxloc(x, dim=1, mask=y>=frac*ypeak)
    if (ileft>1 .and. ileft< size(y)) then
       ! linear interpolation to find x of lefthand half power point
       xleft = x(ileft) -&
            & (y(ileft)-frac*ypeak)/(y(ileft)-y(ileft-1))&
            & *(x(ileft)-x(ileft-1))
    else
       xleft = x(1)
    end if
    ! same for righthand side
    if (iright>1 .and. iright < size(y)) then
       xright =  x(iright) +&
            & (y(iright)-frac*ypeak)/(y(iright)-y(iright+1))&
            & *(x(iright+1)-x(iright))
    else
       xright =  x(max(1,iright)) 
    end if

!!$    fwhm = xright - xleft
    halfmax%a = xleft
    halfmax%b = xright
  end function find_halfmax_d


!!!
!!! Single precision versions
!!!  

  function find_zero_s( y, x ) result (x0)
    real :: x0
    real, intent(in), dimension(:) :: y, x
    x0 = real(find_zero_d( real(y,DP), real(x,DP) ) )
  end function find_zero_s
  

  function find_peak_s( y, x ) result (xpeak)
    real :: xpeak
    real, intent(in), dimension(:) :: y, x
    xpeak = real(find_peak_d( real(y,DP), real(x,DP) ) )
  end function find_peak_s
  
  function find_halfmax_s ( y, x, fraction ) result (halfmax)
    type(pair_s) :: halfmax
    real, intent(in), dimension(:) :: y, x
    real, optional, intent(in) :: fraction
    type(pair) :: halfmax_d
    if (present(fraction)) then
       halfmax_d = find_halfmax_d( real(y,DP), real(x,DP), real(fraction,DP) )
    else
       halfmax_d = find_halfmax_d( real(y,DP), real(x,DP) )
    end if
    halfmax%a = real(halfmax_d%a)
    halfmax%b = real(halfmax_d%b)
  end function find_halfmax_s
  

end module modstats
