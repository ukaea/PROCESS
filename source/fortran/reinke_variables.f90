module reinke_variables
  !! author: S. Muldrew (UKAEA)
  !!
  !! This module contains global variables relating to the minimum impurity fraction
  !! for detached divertor conditions Reinke criterion. It furthermore uses
  !! several parameters from Kallenbach model like netau and empurity_enrichment.
  !!
  !!### References
  !!
  !! - M.L. Reinke 2017 Nucl. Fusion 57 034004

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  integer :: impvardiv
  !! Index of impurity to be iterated for Reinke divertor detachment criterion

  real(dp) :: lhat
  !! Connection length factor L|| = lhat qstar R for Reinke criterion, default value from
  !! Post et al. 1995 J. Nucl. Mat.  220-2 1014

  real(dp) :: fzmin
  !! Minimum impurity fraction necessary for detachment. This is the impurity at the SOL/Div.

  real(dp) :: fzactual
  !! Actual impurity fraction of divertor impurity (impvardiv) in the SoL (taking
  !! impurity_enrichment into account) (`iteration variable 148`)

  integer :: reinke_mode
  !! Switch for Reinke criterion H/I mode:
  !!
  !! - =0 H-mode
  !! - =1 I-mode
end module reinke_variables
