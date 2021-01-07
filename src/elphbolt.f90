program elphBolt
  !! Author: Nakib H. Protik
  !! Summary: Main driver program.
  !!
  !! elphBolt is a program for solving the coupled electron-phonon Boltzmann
  !! transport equations (e-ph BTEs) as formulated in Phys. Rev. B 101, 075202 (2020)
  !! and Phys. Rev. B 102, 245202 (2020) with both the electron-phonon and phonon-phonon
  !! interactions computed ab initio.

  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use electron_module, only: electron
  use phonon_module, only: phonon
  use wannier_module, only: epw_wannier
  use interactions, only:  calculate_g_mixed, calculate_g2_bloch
  
  implicit none
  
  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(epw_wannier) :: wann
  type(electron) :: el
  type(phonon) :: ph
  
  if(this_image() == 1) then
     print*, 'Number of images = ', num_images()
  end if

  !Set up numerics data
  call num%initialize
  
  !Set up crystal
  call crys%initialize

  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)

  !Read EPW Wannier data
  call wann%read

!!$  !Test electron bands, phonon dispersions, and g along path.
!!$  if(this_image() == 1) then
!!$     call wann%test_wannier(crys, num)
!!$  end if
  
  !Calculate electrons
  call el%initialize(wann, crys, sym, num)

  !Calculate phonons
  call ph%initialize(wann, crys, sym, num)

  !Calculate mixed Bloch-Wannier space e-ph vertex
  call calculate_g_mixed(wann, el, num)

  !Calculate Bloch space e-ph vertex
  call calculate_g2_bloch(wann, crys, el, ph, num)

  !Calculate ph-ph vertex
  !TODO call calculate_V_3ph

  !Calculate transition probabilities

  !Iterate BTEs
end program elphBolt
