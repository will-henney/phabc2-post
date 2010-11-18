! 
! Fake a molecular hydrogen fraction

! In the absence of proper treatment of molecules in Phab-C2 we are
! reduced to this simple function of column density, through its proxy:
! V-band extinction, AV. The function has been vaguely calibrated
! against nameless Cloudy models, probably the ones used in the appendix
! of the 2009 magnetic globule paper

! History: 

! 18 Sep 2010 - Ported back to F90 as a mini-module

! 24 Aug 2010 - Package the function into its own module, since I am
! using it all over the place

! """
module mod_molfrac
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real, private :: mol_AV0 = 3.0         ! position of molecular transition 
  real, private :: mol_sharpness = 4.0   ! sharpness of molecular transition
  interface molfrac
     module procedure molfrac_sp, molfrac_dp
  end  interface molfrac
contains

  elemental function molfrac_sp(AV) result (molfrac)
    real :: molfrac
    real, intent(in) :: AV
    molfrac = 1.0/(1.0 + exp(-mol_sharpness*(AV-mol_AV0)))
  end function molfrac_sp

  elemental function molfrac_dp(AV) result (molfrac)
    real(dp) :: molfrac
    real(dp), intent(in) :: AV
    molfrac = 1.0_dp/(1.0_dp + exp(-real(mol_sharpness,dp)*(AV-real(mol_AV0,dp))))
  end function molfrac_dp

end module mod_molfrac
