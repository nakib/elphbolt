! Copyright 2022 elphbolt contributors.
! This file is part of elphbolt <https://github.com/nakib/elphbolt>.
!
! elphbolt is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! elphbolt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with elphbolt. If not, see <http://www.gnu.org/licenses/>.

program superconda
  !! Author: Nakib Haider Protik
  !! Summary: Superconductivity solver app.
  !! Version: 
  !!
  !! superconda is the superconductivity suite of elphbolt.
  !! This is a solver for the Migdal-Eliashberg equations as described in
  !! E. R. Margine and F. Giustino, Phys. Rev. B 87, 024505 (2013).

  use misc, only: print_message, subtitle, timer, exit_with_message
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use electron_module, only: electron
  use phonon_module, only: phonon
  use wannier_module, only: wannier
  use MigEl_sc_module, only: migel_sc
  use bz_sums, only: calculate_el_dos_fermi, calculate_el_dos_fermi_Gaussian, &
       calculate_el_Ws, calculate_el_Ws_Gaussian
  use interactions, only: calculate_gkRp, calculate_eph_interaction_ibzk 
  use eliashberg, only: calculate_a2F
  
  implicit none

  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(wannier) :: wann
  type(electron) :: el
  type(phonon) :: ph
  type(MigEl_sc) :: migel
  type(timer) :: t_all, t_event

  !Print banner and other information
  call welcome

  call t_all%start_timer('superconda')

  call t_event%start_timer('Initialization')

  !Set up crystal
  call crys%initialize

  !Set up numerics data
  call num%initialize(crys)

  !For now...
  if(num%runlevel /= 3) call exit_with_message('Set runlevel = 3 for superconductivity mode. Exiting.')
  sync all

  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)

  !Read EPW Wannier data
  call wann%read(num)

  !Calculate electrons
  call el%initialize(wann, crys, sym, num)
  
  call t_event%end_timer('Initialization')

  if(num%plot_along_path) then
     call t_event%start_timer('Plots along path')

     call subtitle("Plotting along high-symmetry path...")

     !Plot electron bands, phonon dispersions, and g along path.
     call wann%plot_along_path(crys, num, el%scissor)

     call t_event%end_timer('Plots along path')
  end if
  
  call t_event%start_timer('Phonons')

  !Calculate phonons
  if(num%use_Wannier_ifc2s) then
     call ph%initialize(crys, sym, num, wann)
  else
     call ph%initialize(crys, sym, num)
  end if

  call t_event%end_timer('Phonons')

  call t_event%start_timer('Migdal-Eliashberg setup')

  !Initialize Migdal-Eliashberg environment
  call migel%initialize(maxval(ph%ens(:,:)))

  call t_event%end_timer('Migdal-Eliashberg setup')

  call t_event%start_timer('Density of states')

  call subtitle("Calculating density of states...")

  !Calculate electron density of states at the Fermi level
  !call calculate_el_dos_Fermi(el, num%tetrahedra)
  call calculate_el_dos_Fermi_Gaussian(el, crys%reclattvecs)

  !Calculate the scaled electron delta functions
  !call calculate_el_Ws(el, num%tetrahedra)
  call calculate_el_Ws_Gaussian(el, crys%reclattvecs)

  call t_event%end_timer('Density of states')

  if(.not. num%read_gk2) then
     call t_event%start_timer('IBZ k e-ph interactions')

     call subtitle("Calculating e-ph interactions...")

     !Calculate mixed Bloch-Wannier space e-ph vertex g(k,Rp)
     call calculate_gkRp(wann, el, num)

     !Calculate Bloch space e-ph vertex g(k,q) for IBZ k
     call calculate_eph_interaction_ibzk(wann, crys, el, ph, num, 'g')

     call t_event%end_timer('IBZ k e-ph interactions')
  end if

  !Deallocate Wannier quantities
  call wann%deallocate_wannier(num)
  
  !After this point the electron eigenvectors are not needed
  call el%deallocate_eigenvecs

  call t_event%start_timer('IBZ a2F')

  call subtitle("Calculating a2F for all IBZ states...")

  !Calculate anisotropic a2F for all IBZ states
  call calculate_a2F(wann, el, ph, num, migel%omegas, migel%iso_lambda0, migel%omegalog, &
       migel%use_external_eps)

  call t_event%end_timer('IBZ a2F')

  if(migel%iso_lambda0 <= migel%mustar) then
     call print_message("There is no superconductivity since e-ph coupling <= Coulomb pseudopotential.")
  else
     call t_event%start_timer('Superconductivity')

     call subtitle("Solving Migdal-Eliashberg equations...")

     !Calculate McMillan-Allen-Dynes theory
     call migel%calculate_MAD_theory

     !Calculate Migdal-Eliashberg theory
     call migel%calculate_MigEl_theory(el, wann, num, maxval(ph%ens(:,:)))

     call t_event%end_timer('Superconductivity')
  end if

  call t_all%end_timer('superconda')

  call print_message('______________________Thanks for using elphbolt->superconda. Bye!______________________')

contains

  subroutine welcome
    !! Subroutine to print a pretty banner.

    if(this_image() == 1) then
       write(*,'(A75)') "+-------------------------------------------------------------------------+"
       write(*,'(A75)') "| \                                                                       |"
       write(*,'(A75)') "|  \                                                                      |"
       write(*,'(A75)') "|   \   \                                                                 |"
       write(*,'(A75)') "|    \   \                                                                |"
       write(*,'(A75)') "|   __\   \                                                     _         |"                               
       write(*,'(A75)') "|   \      \     ___ _   _ _ __   ___  _ ___ ___ ___  _ __   __|.| __ _   |"
       write(*,'(A75)') "|    \    __\   / __| | | |.'_ \ / _ \|.///// __/ _ \| '_ \ / _`.|/ _` |  |"
       write(*,'(A75)') "|     \  \      \__ \ |_| |.|_) :  __/|.|  | (_| (_) : | | : (_|.| (_| |  |"
       write(*,'(A75)') "|      \ \      |___/\__,_|. __/ \___||_|   \___\___/|_| |_|\__,_|\__,_|  |"
       write(*,'(A75)') "|       \ \               |_|                                             |"                          
       write(*,'(A75)') "|        \\                                                               |"
       write(*,'(A75)') "|         \\                                                              |"
       write(*,'(A75)') "|          \                                                              |"
       write(*,'(A75)') "|           \                                                             |"
       write(*,'(A75)') "| A solver for the Eliashberg superconductivity equations.                |"
       write(*,'(A75)') "| Copyright 2022 elphbolt contributors.                                   |"
       write(*,'(A75)') "|                                                                         |"
       write(*,'(A75)') "| This is 'free as in freedom'[*] software, distributed under the GPLv3.  |"
       write(*,'(A75)') "| [*] https://www.gnu.org/philosophy/free-sw.en.html                      |"
       write(*,'(A75)') "+-------------------------------------------------------------------------+" 
       print*, ' '

       write(*, '(A, I5)') 'Number of coarray images = ', num_images()
    end if
  end subroutine welcome
end program superconda
