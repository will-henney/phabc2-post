! """
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
  real, private :: mol_AV0 = 3.0         ! position of molecular transition 
  real, private :: mol_sharpness = 4.0   ! sharpness of molecular transition

contains

  real elemental function molfrac(AV)
    real, intent(in) :: AV
    molfrac = 1.0 - 1.0/(1.0 + exp(mol_sharpness*(AV-mol_AV0)))
  end function molfrac

end module mod_molfrac
