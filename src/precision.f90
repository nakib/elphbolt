module precision
   use iso_fortran_env, only: real64, real128, int64

   implicit none
    
   integer, parameter :: r64 = real64
   integer, parameter :: r128 = real128
   integer, parameter :: i64 = int64

end module precision