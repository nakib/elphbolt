module precision
   !! Module containing short hands for precisions of various data type defined in 
   !! `iso_fortran_env`

   use iso_fortran_env, only: real64, real128, int64

   implicit none
    
   ! Floating point precisions
   integer, parameter :: r64 = real64
   !! 64-bit float
   integer, parameter :: r128 = real128
   !! 128-bit float

   ! Integer precisions
   integer, parameter :: i64 = int64
   !! 64-bit integer (single)

end module precision