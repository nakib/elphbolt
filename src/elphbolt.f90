! Copyright (C) 2020- Nakib Haider Protik <nakib.haider.protik@gmail.com>
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

program elphbolt
  !! Author: Nakib Haider Protik
  !! Summary: Main driver program.
  !! Version: 1.0.1
  !!
  !! elphbolt is a program for solving the coupled electron-phonon Boltzmann transport equations
  !! (e-ph BTEs) as formulated in https://arxiv.org/abs/2109.08547 (2021) with both the
  !! electron-phonon and phonon-phonon interactions computed ab initio.
  
  use misc, only: welcome, print_message, subtitle, timer, exit_with_message
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use electron_module, only: electron
  use phonon_module, only: phonon
  use wannier_module, only: epw_wannier
  use MigEl_sc_module, only: migel_sc
  use bte_module, only: bte
  use bz_sums, only: calculate_dos, calculate_qTF, calculate_el_dos_fermi, calculate_el_Ws
  use interactions, only: calculate_gReq, calculate_gkRp, calculate_3ph_interaction, &
       calculate_eph_interaction_ibzq, calculate_eph_interaction_ibzk, &
       calculate_echimp_interaction_ibzk, calculate_bound_scatt_rates, &
       calculate_defect_scatt_rates
  use eliashberg, only: calculate_a2F
  use phonon_defect_module, only: phonon_defect
  use Green_function, only: calculate_retarded_phonon_D0
  
  implicit none
  
  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(epw_wannier) :: wann
  type(electron) :: el
  type(phonon) :: ph
  type(bte) :: bt
  type(MigEl_sc) :: migel
  type(phonon_defect) :: ph_def
  type(timer) :: t_all, t_event
  
  !Print banner and other information
  call welcome

  call t_all%start_timer('elphbolt')

  call t_event%start_timer('Initialization')
  
  !Set up crystal
  call crys%initialize
  
  !Set up numerics data
  call num%initialize(crys)
  
  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)

  if(num%onlyebte .or. num%drag .or. num%phe &
       .or. num%plot_along_path .or. num%runlevel == 3) then
     !Read EPW Wannier data
     call wann%read(num)

     !Calculate electrons
     call el%initialize(wann, crys, sym, num)
  end if
  
  call t_event%end_timer('Initialization')
  
  call t_event%start_timer('Phonons')
  
  !Calculate phonons
  call ph%initialize(crys, sym, num)

  call t_event%end_timer('Phonons')

  !Create phonon defect
  if(num%phdef_Tmat) call ph_def%initialize(ph, crys)
  
  select case(num%runlevel)
  case(1) !BTE workflow
     call t_event%start_timer('Density of states and one-particle scattering rates')

     call subtitle("Calculating density of states...")
     if(num%onlyebte .or. num%drag) then
        !Calculate electron density of states
        call calculate_dos(el, num%tetrahedra)

        !Calculate Thomas-Fermi screening
        call calculate_qTF(crys, el)

        !Calculate boundary scattering rates.
        call calculate_bound_scatt_rates(el%prefix, num%elbound, crys%bound_length, &
             el%vels, el%indexlist_irred, bt%el_rta_rates_bound_ibz)
     end if

     !Calculate phonon density of states and, if needed, phonon-isotope
     !and/or phonon-substitution scattering rates.
     call calculate_dos(ph, num%tetrahedra, crys%gfactors, crys%subs_gfactors, &
          crys%atomtypes, bt%ph_rta_rates_iso_ibz, bt%ph_rta_rates_subs_ibz, &
          num%phiso, num%phsubs)

     !Calculate boundary scattering rates.
     call calculate_bound_scatt_rates(ph%prefix, num%phbound, crys%bound_length, &
          ph%vels, ph%indexlist_irred, bt%ph_rta_rates_bound_ibz)

     call t_event%end_timer('Density of states and one-particle scattering rates')
     
     if(num%plot_along_path) then
        call t_event%start_timer('Plots along path')
        
        call subtitle("Plotting along high-symmetry path...")

        !Plot electron bands, phonon dispersions, and g along path.
        call wann%plot_along_path(crys, num)

        call t_event%end_timer('Plots along path')
     end if
     
     call subtitle("Calculating interactions...")

     !Set chemical potential dependent directory
     call num%create_chempot_dirs(el%chempot)

     !
     if(num%phdef_Tmat) then
        !Calculate phonon-defect interactions
        call t_event%start_timer("Phonon-defect transition rates")

        !Calculate phonon Green's function on defective space
        call calculate_retarded_phonon_D0(ph, crys, ph_def%cell_pos_intvec, ph_def%pcell_atom_label, ph_def%D0, &
             ph_def%dimp_cell_pos_intvec, ph_def%pcell_atom_dof)

        call ph_def%calculate_phonon_Tmatrix(ph, crys)
        
        call t_event%end_timer("Phonon-defect transition rates")
        
        call exit
     end if
     !!
     
     if(num%onlyphbte .and. num%phe .or. num%drag) then
        if(.not. num%read_gq2) then
           call t_event%start_timer('IBZ q e-ph interactions')
           
           !Calculate mixed Bloch-Wannier space e-ph vertex g(Re,q)
           call calculate_gReq(wann, ph, num)
           
           !Calculate Bloch space e-ph vertex g(k,q) for IBZ q
           call calculate_eph_interaction_ibzq(wann, crys, el, ph, num, 'g')

           call t_event%end_timer('IBZ q e-ph interactions')
        end if
        
        call t_event%start_timer('IBZ ph-e transition probilities')
        
        !Calculate ph-e transition probabilities
        call calculate_eph_interaction_ibzq(wann, crys, el, ph, num, 'Y')
                
        call t_event%end_timer('IBZ ph-e transition probilities')
     end if

     if(num%onlyebte .or. num%drag) then
        if(.not. num%read_gk2) then
           call t_event%start_timer('IBZ k e-ph interactions')
           
           !Calculate mixed Bloch-Wannier space e-ph vertex g(k,Rp)
           call calculate_gkRp(wann, el, num)

           !Calculate Bloch space e-ph vertex g(k,q) for IBZ k
           call calculate_eph_interaction_ibzk(wann, crys, el, ph, num, 'g')

           call t_event%end_timer('IBZ k e-ph interactions')
        end if
        
        call t_event%start_timer('IBZ e-ph transition probabilities')

        !Calculate e-ph transition probabilities
        call calculate_eph_interaction_ibzk(wann, crys, el, ph, num, 'X')
        
        call t_event%end_timer('IBZ e-ph transition probabilities')
     end if

     if(num%onlyebte .or. num%drag .or. num%phe .or. num%drag &
          .or. num%plot_along_path) then
        !Deallocate Wannier quantities
        call wann%deallocate_wannier(num)
     end if

     if(num%onlyebte .or. num%drag .or. num%phe) then
        !After this point the electron eigenvectors are not needed
        call el%deallocate_eigenvecs
     end if

     if(num%onlyebte .or. num%drag) then
        if(num%elchimp) then
           call t_event%start_timer('e-ch. imp. interactions')
           
           !Calculate e-ch. imp. transition probabilities
           call calculate_echimp_interaction_ibzk(crys, el, num)

           call t_event%end_timer('e-ch. imp. interactions')
        end if
     end if

     if(num%onlyphbte .or. num%drag) then
        if(.not. num%read_V) then
           call t_event%start_timer('IBZ q ph-ph interactions')
           
           !Calculate ph-ph vertex
           call calculate_3ph_interaction(ph, crys, num, 'V')

           call t_event%end_timer('IBZ q ph-ph interactions')
        end if

        if(.not. num%read_W) then
           call t_event%start_timer('IBZ ph-ph scattering rates')
           
           !Calculate ph-ph transition probabilities
           call calculate_3ph_interaction(ph, crys, num, 'W')
           
           call t_event%end_timer('IBZ ph-ph scattering rates')
        end if
     end if

     if(num%onlyphbte .or. num%drag) then
        !After this point the phonon eigenvectors and other quantities are not needed
        call ph%deallocate_phonon_quantities
     end if
     
     !Solve BTEs
     if(num%onlyphbte .and. .not. num%phe) then
        call bt%solve_bte(num, crys, sym, ph)
     else
        call bt%solve_bte(num, crys, sym, ph, el)
     end if
  case(2) !BTE Post-processing case
     call subtitle("Post-processing...")
        
     !Read RTA response functions from finished calculation
     if(num%onlyphbte .and. .not. num%phe) then
        call bt%post_process(num, crys, sym, ph)
     else
        call bt%post_process(num, crys, sym, ph, el)
     end if
  case(3) !Superconductivity mode worklflow
     call t_event%start_timer('Migdal-Eliashberg setup')

     !Initialize Migdal-Eliashberg environment
     call migel%initialize(maxval(ph%ens(:,:)))

     call t_event%end_timer('Migdal-Eliashberg setup')
     
     call t_event%start_timer('Density of states')
     
     call subtitle("Calculating density of states...")
     
     !Calculate electron density of states at the Fermi level
     call calculate_el_dos_Fermi(el, num%tetrahedra)
     
     !Calculate the scaled electron delta functions
     call calculate_el_Ws(el, num%tetrahedra)
     
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

     !After this point we can release some of the phonon quantities
     call ph%deallocate_phonon_quantities

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
  case default
     call exit_with_message('Unknown runlevel. Exiting.')
  end select

  call t_all%end_timer('elphbolt')
  
  call print_message('______________________Thanks for using elphbolt. Bye!______________________')
end program elphbolt
