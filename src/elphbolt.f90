program elphBolt
  !! Author: Nakib H. Protik
  !! Summary: Main driver program.
  !!
  !! elphBolt is a program for solving the coupled electron-phonon Boltzmann
  !! transport equations (e-ph BTEs) as formulated in Phys. Rev. B 101, 075202 (2020)
  !! and Phys. Rev. B 102, 245202 (2020) with both the electron-phonon and phonon-phonon
  !! interactions computed ab initio.

  use params, only: k4, dp
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use electron_module, only: electron
  use phonon_module, only: phonon
  use wannier_module, only: epw_wannier
  use bte_module, only: bte
  use bz_sums, only: calculate_dos, calculate_chempot
  use interactions, only: calculate_gReq, calculate_gkRp, calculate_3ph_interaction, &
       calculate_eph_interaction_ibzq, calculate_eph_interaction_ibzk
  
  implicit none
  
  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(epw_wannier) :: wann
  type(electron) :: el
  type(phonon) :: ph
  type(bte) :: bt

  if(this_image() == 1) then
     print*, 'Number of images = ', num_images()
  end if

  !Set up crystal
  call crys%initialize
  
  !Set up numerics data
  call num%initialize
  
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

  !Calculate electron density of states
  call calculate_dos(el, num%tetrahedra)

  !Calculate chemical potential
  call calculate_chempot(el, crys%T, crys%volume)
  
  !Calculate phonons
  call ph%initialize(wann, crys, sym, num)

  !Calculate phonon density of states
  call calculate_dos(ph, num%tetrahedra)

  if(num%phe) then
     if(.not. num%read_gq2) then
        !Calculate mixed Bloch-Wannier space e-ph vertex g(Re,q)
        call calculate_gReq(wann, ph, num)

        !Calculate Bloch space e-ph vertex g(k,q) for IBZ q
        call calculate_eph_interaction_ibzq(wann, crys, el, ph, num, 'g')
     end if
     
     !Calculate ph-e transition probabilities
     call calculate_eph_interaction_ibzq(wann, crys, el, ph, num, 'Y')
  end if

  if(.not. num%read_gk2) then
     !Calculate mixed Bloch-Wannier space e-ph vertex g(k,Rp)
     call calculate_gkRp(wann, el, num)

     !Calculate Bloch space e-ph vertex g(k,q) for IBZ k
     call calculate_eph_interaction_ibzk(wann, crys, el, ph, num, 'g')
  end if

  !Calculate ph-e transition probabilities
  call calculate_eph_interaction_ibzk(wann, crys, el, ph, num, 'X')
     
  if(.not. num%read_V) then
     !Calculate ph-ph vertex
     call calculate_3ph_interaction(ph, crys, num, 'V')
  end if
  
  !Calculate ph-ph transition probabilities
  call calculate_3ph_interaction(ph, crys, num, 'W')
  
  !RTA solution
  if(num%phbte) then
     if(num%phe) then
        call bt%solve_bte(num, crys, sym, ph, el)
     else
        call bt%solve_bte(num, crys, sym, ph)
     end if
  end if

  if(num%ebte) then
     call bt%solve_bte(num, crys, sym, ph, el)
  end if
end program elphBolt
