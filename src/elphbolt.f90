program elphBolt
  !! Author: Nakib H. Protik
  !! Summary: Main driver program.
  !!
  !! elphBolt is a program for solving the coupled electron-phonon Boltzmann
  !! transport equations (e-ph BTEs) as formulated in Phys. Rev. B 101, 075202 (2020)
  !! and arXiv:2008.08722 (2020) with both the electron-phonon and phonon-phonon
  !! interactions computed ab initio.

  use config, only: initialize_system
  !use mesh, only: calculate_electrons, calculate_phonons
  !use wannier, only: calculate_g_mixed, calculate_g_bloch, gmixed_epw_gamma
  use wannier, only: epw_wannier
  
  implicit none

  type(epw_wannier) :: wann
  
  if(this_image() == 1) then
     print*, 'Number of images = ', num_images()
  end if


  call wann%read
  
  !Read inputs, find symmetry operations, and set up calculation environment.
  call initialize_system

  !Test electron bands, phonon dispersions, and g along path.
  if(this_image() == 1) then
     call wann%test_wannier
  end if
  
!!$  !Calculate electrons
!!$  call calculate_electrons
!!$
!!$  !Calculate phonons
!!$  call calculate_phonons
!!$
!!$  !Calculate mixed Bloch-Wannier space e-ph vertex
!!$  call calculate_g_mixed
!!$
!!$  !Calculate Bloch space e-ph vertex
!!$  call calculate_g_bloch

  !Calculate ph-ph vertex
  !TODO call calculate_V

  !Calculate transition probabilities

  !Iterate BTEs
end program elphBolt
