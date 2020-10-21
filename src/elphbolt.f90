program elphBolt
  !! Author: Nakib H. Protik
  !! Summary: Main driver program.
  !!
  !! elphBolt is a program for solving the coupled electron-phonon Boltzmann
  !! transport equations (e-ph BTEs) as formulated in Phys. Rev. B 101, 075202 (2020)
  !! with both the electron-phonon and phonon-phonon interactions computed ab initio.

  use derived_types, only: crystal_data, reciprocal_lattice_data, symmetry_data, &
       electron_data, phonon_data, control_data, EPW_data
  use config, only: initialize_system
  use wannier, only: read_EPW_Wannier
  
  implicit none

  type(crystal_data) :: crys
  !! Contains the crystal information.
  type(reciprocal_lattice_data) :: reclat
  !! Contains the reciprocal lattice information.
  type(symmetry_data) :: sym
  !! Contains the symmetry information.
  type(electron_data) :: el
  !! Contains the electronic information.
  type(phonon_data) :: ph
  !! Contains the phononic information.
  type(control_data) :: con
  !! Contains the system i/o and run-specific information.
  type(EPW_data) :: epw
  !! Contains the epw information.
  
  if(this_image() == 1) then
     print*, 'Number of images = ', num_images()
  end if

  !Read inputs, find symmetry operations, and set up calculation environment.
  call initialize_system(crys,reclat,sym,el,ph,con)
  
  !Read EPW Wannier data.
  call read_EPW_Wannier(epw)
  
  !Calculate electrons
  !TODO call calculate_electrons(el)

  !Calculate phonons
  !TODO call calculate_phonons(ph)

  !Calculate real space e-ph vertex
  !TODO call calculate_g_real(wann,eph)

  !Calculate Bloch space e-ph vertex
  !TODO call calculate_g_bloch(wann,eph)

  !Calculate ph-ph vertex
  !TODO call calculate_V

  !Calculate transition probabilities

  !Iterate BTEs
end program elphBolt
