program V3_permutations
  !! TODO add description
  
  use precision, only: i64, r64
  use misc, only: print_message, subtitle, timer, exit_with_message, mux_vector, demux_state, &
       mux_state, twonorm, permutations
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use phonon_module, only: phonon
      
  implicit none
  
  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(phonon) :: ph
  !type(timer) :: t_event

  if(this_image() == 1) then
     write(*, '(A)')  'V3 permutations playground'
     write(*, '(A, I5)') 'Number of coarray images = ', num_images()
  end if
     
  !Set up crystal
  call crys%initialize

  !Set up numerics data
  call num%initialize(crys)

  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)

  !Calculate phonons
  call ph%initialize(crys, sym, num)

  !Generate the permutations of interacting triplets
  !TODO

contains

  !your functions/subroutine
  
end program V3_permutations
