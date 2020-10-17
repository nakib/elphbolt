program elphbolt
  !! Author: Nakib H. Protik
  !! Summary: Main driver program.
  !!
  !! elphBolt is a program for solving the coupled electron-phonon Boltzmann
  !! transport equations (e-ph BTEs) as formulated in Phys. Rev. B 101, 075202 (2020).

  use derived_types, only: crystal_data, reciprocal_lattice_data, symmetry_data, electron_data
  use config, only: initialize_system
  
  implicit none

  type(crystal_data) :: crys
  !! Contains the crystal information.
  type(reciprocal_lattice_data) :: reclat
  !! Contains the reciprocal lattice information.
  type(symmetry_data) :: sym
  !! Contains the symmetry information.
  type(electron_data) :: el
  !! Contains the electronic information.
  
  call initialize_system(crys,reclat,sym,el)

  if(this_image() == 1) then
     print*, 'Number of images = ', num_images()
!!$  print*, 'Number of elements = ', crystal%numelements
!!$  print*, 'Number of atoms = ', crystal%numatoms
!!$  print*, 'Elements = ', crystal%elements
!!$  print*, 'Atom types = ', crystal%atom_types
!!$  print*, 'Lattice vectors', crystal%lattvecs
!!$  print*, 'Basis', crystal%basis
!!$  print*, 'Basis (Cartesian)', crystal%basis_cart
     print*, 'Volume = ', crys%volume, 'nm^3'
     print*, 'BZ Volume = ', reclat%volume_bz, 'nm^-3'
     print*, 'Reciprocal lattice vectors ', reclat%reclattvecs
  end if
  
  !call get_crystal(crystal)
  
  !Read input
  
  !Read read electron and phonon real space data

  !Calculate electrons

  !Calculate phonons

  !Calculate e-ph vertex

  !Calculate ph-ph vertex

  !Calculate transition probabilities

  !Iterate BTEs
end program elphbolt
