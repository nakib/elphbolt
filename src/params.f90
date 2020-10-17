module params
  !! Module containing various parameters and constants.
  
  implicit none

  integer, parameter :: dp = selected_real_kind(14,200)
  !! Custom real precision.
  integer, parameter :: k4 = selected_real_kind(8)
  !! Custom integer precision.

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  !! Value of pi
  real(dp), parameter :: twopi = 2.0_dp*pi
  !! Value of 2pi
end module params
