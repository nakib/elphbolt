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

module electron_module
  !! Module containing types and procedures related to the electronic properties.

  use params, only: dp, k8
  use particle_module, only: particle
  use misc, only: exit_with_message, print_message, demux_state, sort, &
       binsearch, subtitle, Fermi, write2file_rank2_real, write2file_rank3_real 
  use numerics_module, only: numerics
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal, calculate_wavevectors_full
  use symmetry_module, only: symmetry, find_irred_wedge, create_fbz2ibz_map
  use delta, only: form_tetrahedra_3d, fill_tetrahedra_3d, form_triangles, &
       fill_triangles
  
  implicit none

  private
  public electron

  type, extends (particle) :: electron
     !! Data and procedures related to the electronic properties.

     character(len = 2) :: prefix = 'el'
     !! Prefix idenitfying particle type.
     integer(k8) :: spindeg
     !! Spin degeneracy.
     integer(k8) :: numtransbands
     !! Total number of transport active bands.
     integer(k8) :: indlowband
     !! Lowest transport band index.
     integer(k8) :: indhighband
     !! Highest transport band index.
     integer(k8) :: indlowconduction
     !! Lowest conduction band index.
     integer(k8) :: indhighvalence
     !! Highest valence band index.
     integer(k8), allocatable :: bandlist(:)
     !! List of transport active band indices.
     integer(k8) :: mesh_ref
     !! Electron mesh refinement factor compared to the phonon mesh.
     integer(k8) :: mesh_ref_array(3)
     !! The same as above, but in array form. This is useful for 3d vs 2d cases.
     integer(k8) :: nstates_inwindow
     !! Number of electron wave vectors within transport window.
     integer(k8) :: nstates_irred_inwindow
     !! Number of IBZ wedge electron wave vectors within transport window.
     integer(k8), allocatable :: IBZ_inwindow_states(:,:)
     !! List of irreducible wedge states within transport window.
     real(dp) :: enref
     !! Electron reference energy (eV).
     !! This is the center of the transport energy window.
     real(dp) :: fsthick
     !! Fermi surface thickness (eV).
     real(dp) :: chempot
     !! Chemical potential in (eV).
     real(dp), allocatable :: conc(:)
     !! Band resolved carrier concentration.
     real(dp) :: conc_el
     !! Total electron carrier concentration.
     real(dp) :: conc_hole
     !! Total hole carrier concentration.
     real(dp) :: chimp_conc_n
     !! Concentration of donor impurities.
     real(dp) :: chimp_conc_p
     !! Concentration of acceptor impurities.
     real(dp) :: Zn
     !! Ionization number of donor dopant.
     real(dp) :: Zp
     !! Ionization number of acceptor dopant.
     logical :: metallic
     !! Is the system metallic?
     character(len = 1) :: dopingtype
     !! Type of doping. This is needed for runlevel 0 only.
     integer(k8) :: numconc
     !! Number of concentration points. This is needed for runlevel 0 only.
     real(dp), allocatable :: conclist(:)
     !! List of concentrations. This is needed for runlevel 0 only.
     integer(k8) :: numT
     !! Number of temperature points. This is needed for runlevel 0 only.
     real(dp), allocatable :: Tlist(:)
     !! List of temperatures. This is needed for runlevel 0 only.
     real(dp) :: spinnormed_dos_fermi
     !! Spin-normalized density of states at the Fermi level
     real(dp), allocatable :: Ws_irred(:, :), Ws(:, :)
     !! Electron delta functions normalized by spinnormed_dos_fermi
     
   contains

     procedure, public :: initialize=>read_input_and_setup, deallocate_eigenvecs
     procedure, private :: calculate_electrons, calculate_carrier_conc, &
          calculate_chempot
  end type electron

contains

  subroutine read_input_and_setup(self, wann, crys, sym, num)
    !! Read input file and setup groundstate electronic system.

    class(electron), intent(out) :: self
    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(numerics), intent(in) :: num

    !Local variables
    real(dp) :: enref, Zn, Zp, chempot
    real(dp), allocatable :: Tlist(:), conclist(:)
    integer(k8) :: ib, spindeg, numbands, indlowband, indhighband, &
         indlowconduction, indhighvalence, numT, numconc
    logical :: metallic
    character(len = 6) :: concunits
    character(len = 1) :: dopingtype

    namelist /electrons/ enref, spindeg, numbands, &
         indlowband, indhighband, metallic, chempot, Zn, Zp, &
         indlowconduction, indhighvalence, dopingtype, numT, numconc, &
         Tlist, conclist
         
    call subtitle("Setting up electrons...")
    
    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Read electrons information
    spindeg = 2 !Default calculation is non-spin polarized
    numbands = 0
    indlowband = 0
    indhighband = 0
    indlowconduction = 0
    indhighvalence = 0
    metallic = .false.
    Zn = 0.0_dp
    Zp = 0.0_dp
    chempot = -999999.99999_dp !Something crazy
    enref = -999999.99999_dp !Something crazy
    numT = 100 !Something crazy big
    numconc = 100 !Something crazy big
    dopingtype = 'x'
    allocate(Tlist(nuMT), conclist(numconc))
    Tlist = -1.0_dp !Something crazy
    conclist = 0.0_dp
    read(1, nml = electrons)
    if(spindeg < 1 .or. spindeg > 2) then
       call exit_with_message('spindeg can be 1 or 2.')
    end if
    if(numbands < 1) then
       call exit_with_message('numbands should be > 0.')
    end if
    if(indlowband < 1) then
       call exit_with_message('indlowband should be > 0.')
    end if
    if(indhighband < 1) then
       call exit_with_message('indhighband should be > 0.')
    end if
    if(.not. metallic) then
       if(indlowconduction < 1 .and. indhighvalence < 1) then
          call exit_with_message(&
               'For non-metals, must provide lowest conduction or highest valence band.')
       end if
    end if
    if(num%runlevel == 0) then
       if(numT <= 0 .or. numconc <= 0) then
          call exit_with_message('numT or numconc should be > 0.')
       end if
       if(numT > 100 .or. numconc > 100) then
          call exit_with_message('numT or numconc > 1000 is not supported.')
       end if
       if(any(Tlist(1:numT) <= 0.0_dp)) then
          call exit_with_message('Unphysical Tlist provided.')
       end if
       if(dopingtype /= 'n' .and. dopingtype /= 'p') then
          print*, dopingtype, len(dopingtype)
          call exit_with_message("dopingtype must be 'n' or 'p'.")
       end if
    end if
    
    self%spindeg = spindeg
    self%numbands = numbands
    self%indlowband = indlowband
    self%indhighband = indhighband
    self%numtransbands = self%indhighband - self%indlowband + 1
    allocate(self%bandlist(self%numtransbands))
    do ib = 1, self%numtransbands
       self%bandlist(ib) = indlowband + ib - 1
    end do
    self%metallic = metallic
    self%indlowconduction = indlowconduction
    self%indhighvalence = indhighvalence
    self%enref = enref
    self%chempot = chempot
    self%Zn = Zn
    self%Zp = Zp
    if(self%metallic) then
       self%Zn = 0
       self%Zp = 0
    end if
    if(num%runlevel == 0) then
       self%numT = numT
       self%numconc = numconc
       allocate(self%Tlist(self%numT), self%conclist(self%numconc))
       self%Tlist(:) = Tlist(1:numT)
       self%conclist(:) = conclist(1:numconc)
       self%dopingtype = dopingtype
    end if
    
    !Close input file
    close(1)

    !Set some electronic properties from the numerics object
    self%mesh_ref = num%mesh_ref
    self%mesh_ref_array = (/num%mesh_ref, num%mesh_ref, num%mesh_ref/)
    if(crys%twod) then
       self%wvmesh(3) = 1_k8
       self%mesh_ref_array(3) = 1_k8
    end if
    self%wvmesh = self%mesh_ref_array*num%qmesh
    self%fsthick = num%fsthick
    
    !Print out information.
    if(this_image() == 1) then
       write(*, "(A, I1)") "Spin degeneracy = ", self%spindeg
       write(*, "(A, I5)") "Number of Wannier electronic bands = ", self%numbands
       write(*, "(A, I5)") "Number of transport active electronic bands = ", self%numtransbands
       write(*, "(A, I5, I5)") "Lowest and highest transport active electronic bands = ", &
            self%bandlist(1), self%bandlist(self%numtransbands)
       write(*, "(A, 1E16.8, A)") "Reference electron energy = ", self%enref, ' eV'
       write(*, "(A, L)") "System is metallic: ", self%metallic
       if(indlowconduction > 0) then
          write(*, "(A, I5)") "Lowest conduction band index = ", self%indlowconduction
       end if
       if(indhighvalence > 0) then
          write(*, "(A, I5)") "Highest valence band index = ", self%indhighvalence
       end if
    end if
    
    !Calculate electrons
    call calculate_electrons(self, wann, crys, sym, num)
    
    !Set total number of charged impurities
    if(.not. self%metallic) then
       self%chimp_conc_n = 0.0_dp
       self%chimp_conc_p = 0.0_dp
       if(self%Zn > 0) self%chimp_conc_n = self%chimp_conc_n + self%conc_el/self%Zn
       if(self%Zp > 0) self%chimp_conc_p = self%chimp_conc_p + self%conc_hole/self%Zp
    end if

    !Print out information.
    call print_message("Electron calculations summary:")
    call print_message("------------------------------")
    if(this_image() == 1) then
       if(crys%twod) then
          concunits = ' cm^-2'
       else
          concunits = ' cm^-3'
       end if
       write(*, "(A, 1E16.8, A)") "Chemical potential = ", self%chempot, ' eV'
       if(.not. self%metallic) then
          write(*, "(A, 1E16.8)") 'Band resolved carrier concentration (+/- = hole/electron):'
          do ib = self%indlowband, self%indhighband
             write(*, "(A, I5, A, 1E16.8, A)") ' Band: ', ib, ', concentration: ', &
                  self%conc(ib), concunits
          end do
          write(*, "(A, 1E16.8, A)") "Absolute total electron concentration = ", self%conc_el, &
               concunits
          write(*, "(A, 1E16.8, A)") "Absolute total hole concentration = ", self%conc_hole, &
               concunits
          write(*, "(A, 1E16.8)") "Ionization of donor impurity = ", self%Zn
          write(*, "(A, 1E16.8)") "Ionization of acceptor impurity = ", self%Zp
          write(*, "(A, 1E16.8, A)") "Donor impurity concentration = ", self%chimp_conc_n, &
               concunits
          write(*, "(A, 1E16.8, A)") "Acceptor impurity concentration = ", self%chimp_conc_p, &
               concunits
       end if
    end if
  end subroutine read_input_and_setup
  
  subroutine calculate_electrons(self, wann, crys, sym, num)
    !! Calculate electron energy window restricted wave vector meshes
    !! and the electronic properties on them

    class(electron), intent(inout) :: self
    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(numerics), intent(in) :: num
    
    !Some utitlity variables
    integer(k8) :: i, l, s, il, ii, jj, kk, ib, count, istate, aux
    real(dp), allocatable :: el_ens_tmp(:, :), el_vels_tmp(:, :, :)

    !Switch for mesh utilites with or without energy restriction
    logical :: blocks

    !I/O related
    character(len = 1024) :: filename, numcols

    call print_message("Energy unrestricted calculation:")
    call print_message("--------------------------------")
    
    !Set initial FBZ total number of wave vectors
    self%nwv = product(self%wvmesh)
    
    !The electronic mesh setup proceeds in multiple steps:
    ! 1. Calculate full electron wave vector mesh
    call print_message("Calculating FBZ...")
    blocks = .false.
    call calculate_wavevectors_full(self%wvmesh, self%wavevecs, blocks)
    
    ! 2. Calculate the IBZ
    call print_message("Calculating IBZ and IBZ -> FBZ mappings...")
    call find_irred_wedge(self%wvmesh, self%nwv_irred, self%wavevecs_irred, &
         self%indexlist_irred, self%nequiv, sym%nsymm_rot, sym%qrotations, &
         self%ibz2fbz_map, self%equiv_map, blocks)
    
    ! 3. Calculate IBZ quantities
    call print_message("Calculating IBZ energies...")
    allocate(self%ens_irred(self%nwv_irred, wann%numwannbands), &
         self%vels_irred(self%nwv_irred, wann%numwannbands, 3), &
         self%evecs_irred(self%nwv_irred, wann%numwannbands, wann%numwannbands))
    call wann%el_wann_epw(crys, self%nwv_irred, self%wavevecs_irred, self%ens_irred, &
         self%vels_irred, self%evecs_irred)
    
    ! 4. Map out FBZ quantities from IBZ ones
    call print_message("Mapping out FBZ energies...")
    allocate(self%indexlist(self%nwv), self%ens(self%nwv, wann%numwannbands), &
         self%vels(self%nwv, wann%numwannbands, 3))
    
    do i = 1,self%nwv_irred !an irreducible point
       do l = 1,self%nequiv(i) !number of equivalent points of i
          il = self%ibz2fbz_map(l, i, 2) ! (i, l) -> il
          s = self%ibz2fbz_map(l, i, 1) ! mapping rotation

          !index list
          self%indexlist(il) = il

          !energy
          self%ens(il,:) = self%ens_irred(i,:)
          
          !velocity
          do ib = 1, self%numtransbands !wann%numwannbands
             !here use real space (Cartesian) rotations
             self%vels(il, ib, :) = matmul(sym%crotations(:, :, s), self%vels_irred(i, ib, :))
          end do
       end do
    end do

    if(.not. self%metallic) then
       if(num%runlevel == 0) then !Calculate chemical potentials for
          !the given temperatures and concentrations
          if(crys%twod) then
             call calculate_chempot(self, crys%volume, self%dopingtype, self%Tlist, self%conclist, &
                  crys%thickness)
          else
             call calculate_chempot(self, crys%volume, self%dopingtype, self%Tlist, self%conclist)
          end if
          call exit_with_message("Chemical potentials calculated. Runlevel 0 finished. Exiting.")
       else !Calculate carrier concentration for non-metals
          call print_message("Calculating carrier concentrations...")
          if(crys%twod) then
             call calculate_carrier_conc(self, crys%T, crys%volume, crys%thickness)
          else
             call calculate_carrier_conc(self, crys%T, crys%volume)
          end if
       end if
    end if
    
    call print_message("Transport energy window restricted calculation:")
    call print_message("-----------------------------------------------")
    
    ! 5. Find energy window restricted FBZ blocks.
    !    After this step, self%nwv, self%indexlist will refer
    !    to the energy restricted mesh.
    call print_message("Calculating Fermi window restricted FBZ blocks...")
    call apply_energy_window(self%nwv, self%indexlist, self%ens, self%enref, self%fsthick)
    
    ! 6. Sort index list and related quanties of FBZ blocks
    call print_message("Sorting FBZ blocks index list...")
    call sort(self%indexlist)

    ! 7. Get FBZ blocks wave vectors, energies, velocities and eigenvectors.
    !    After this step, self%wavevecs, self%ens, self%vels, and self%evecs
    !    will refer to the energy restricted mesh.
    call print_message("Calcutating FBZ blocks quantities...")
    
    !wave vectors
    deallocate(self%wavevecs)
    
    blocks = .true.
    call calculate_wavevectors_full(self%wvmesh, self%wavevecs, blocks, self%indexlist) !wave vectors

    !Print electron FBZ mesh
    call write2file_rank2_real("el.wavevecs_fbz", self%wavevecs)
    
    !energies and velocities
    call fbz_blocks_quantities(self%indexlist, self%ens, self%vels)

    !Get FBZ blocks eigenvectors from direct calculations since we are
    !not getting these from IBZ quantities via symmetry rotations
    allocate(self%evecs(self%nwv, wann%numwannbands, wann%numwannbands))
    allocate(el_ens_tmp(self%nwv, wann%numwannbands), el_vels_tmp(self%nwv, wann%numwannbands, 3))
    call wann%el_wann_epw(crys, self%nwv, self%wavevecs, el_ens_tmp, el_vels_tmp, self%evecs)
    deallocate(el_ens_tmp, el_vels_tmp) !free up memory
    
    ! 8. Find IBZ of energy window restricted blocks
    !    After this step, self%nwv_irred, self%indexlist_irred, 
    !    self%wavevecs_irred, self%nequiv, and self%ibz2fbz_map
    !    will refer to the energy restricted mesh 
    call print_message("Calculating IBZ blocks...")
    deallocate(self%wavevecs_irred, self%indexlist_irred, self%nequiv, &
         self%ibz2fbz_map, self%equiv_map)
    blocks = .true.
    call find_irred_wedge(self%wvmesh, self%nwv_irred, self%wavevecs_irred, &
         self%indexlist_irred, self%nequiv, sym%nsymm_rot, sym%qrotations, &
         self%ibz2fbz_map, self%equiv_map, blocks, self%indexlist)

    !Print electron IBZ mesh
    call write2file_rank2_real("el.wavevecs_ibz", self%wavevecs_irred)
    
    !Create symmetrizers of wave vector dependent vectors ShengBTE style
    allocate(self%symmetrizers(3, 3, self%nwv))
    self%symmetrizers = 0.0_dp
    do i = 1, self%nwv
       ii = self%indexlist(i)
       kk = 0
       do jj = 1, sym%nsymm
          if(self%equiv_map(jj, i) == ii) then
             self%symmetrizers(:, :, i) = self%symmetrizers(:, :, i) + &
                  sym%crotations_orig(:, :, jj)
             kk = kk + 1
          end if
       end do
       if(kk > 1) then
          self%symmetrizers(:, :, i) = self%symmetrizers(:, :, i)/kk
       end if
    end do
    
    ! 9. Get IBZ blocks energies, velocities, and eigen vectors.
    call print_message("Calcutating IBZ blocks quantities...")
    deallocate(self%ens_irred, self%vels_irred, self%evecs_irred)
    allocate(self%ens_irred(self%nwv_irred, wann%numwannbands), &
         self%vels_irred(self%nwv_irred, wann%numwannbands, 3), &
         self%evecs_irred(self%nwv_irred, wann%numwannbands, wann%numwannbands))
    call wann%el_wann_epw(crys, self%nwv_irred, self%wavevecs_irred, self%ens_irred, &
         self%vels_irred, self%evecs_irred)
    
    ! 10. Calculate the number of FBZ blocks electronic states
    !     available for scattering
    self%nstates_inwindow = 0
    do i = 1,self%nwv !over FBZ blocks
       do ib = 1,wann%numwannbands !bands
          if(abs(self%ens(i, ib) - self%enref) <= self%fsthick) &
               self%nstates_inwindow = self%nstates_inwindow + 1
       end do
    end do
    if(this_image() == 1) write(*, "(A, I10)") &
         " Number of energy restricted FBZ blocks states = ", self%nstates_inwindow

    ! 11. Create FBZ blocks to IBZ blocks map
    call print_message("Calculating FBZ -> IBZ mappings...")
    call create_fbz2ibz_map(self%fbz2ibz_map,self%nwv,self%nwv_irred, &
         self%indexlist,self%nequiv,self%ibz2fbz_map)
    
    do i = 1, self%nwv_irred !IBZ
       do l = 1, self%nequiv(i) !number of equivalent points of i
          il = self%ibz2fbz_map(l, i, 2) ! (i, l) -> il
          s = self%ibz2fbz_map(l, i, 1) ! symmetry
          call binsearch(self%indexlist, il, aux)

          !energy
          self%ens(aux,:) = self%ens_irred(i,:)
          
          !velocity
          do ib = 1,wann%numwannbands
             !here use real space (Cartesian) rotations
             self%vels(aux, ib, :) = matmul(sym%crotations(:, :, s), self%vels_irred(i, ib, :))
          end do
          self%vels(aux,:,:) = transpose(&
               matmul(self%symmetrizers(:,:,aux),transpose(self%vels(aux,:,:))))
       end do
    end do
        
    ! 12. Calculate the number of IBZ electronic states available for scattering
    self%nstates_irred_inwindow = 0
    do istate = 1,self%nwv_irred*wann%numwannbands
       !Demux state index into band (ib) and wave vector (i) indices
       call demux_state(istate, wann%numwannbands, ib, i)
       if(abs(self%ens_irred(i, ib) - self%enref) <= self%fsthick) then 
          self%nstates_irred_inwindow = self%nstates_irred_inwindow + 1
       end if
    end do
    if(this_image() == 1) write(*, "(A, I10)") " Number of energy restricted IBZ blocks states = ", &
         self%nstates_irred_inwindow
    
    !Calculate list of IBZ in-window states = (wave vector index, band index)
    allocate(self%IBZ_inwindow_states(self%nstates_irred_inwindow,2))
    count = 0
    do istate = 1, self%nwv_irred*wann%numwannbands
       !Demux state index into band (ib) and wave vector (i) indices
       call demux_state(istate, wann%numwannbands, ib, i)
       if(abs(self%ens_irred(i, ib) - self%enref) <= self%fsthick) then
          count = count + 1
          self%IBZ_inwindow_states(count,:) = [i, ib]
       end if
    end do

    !Write IBZ in-window states as text data to file
    if(this_image() == 1) then
       call chdir(num%cwd)
       filename = 'el.inwindow_states_ibz'
       write(numcols,"(I0)") 2
       open(1,file=trim(filename),status='replace')
       write(1,*) "#k-vec index     band index"
       do i = 1, self%nstates_irred_inwindow
          write(1,"("//trim(adjustl(numcols))//"I10)") self%IBZ_inwindow_states(i,:)
       end do
       close(1)
    end if
    !Deallocating this here since this is not used later in the program
    deallocate(self%IBZ_inwindow_states)

    !Print out irreducible electron energies and velocities
    call write2file_rank2_real("el.ens_ibz", self%ens_irred)
    call write2file_rank3_real("el.vels_ibz", self%vels_irred)
    
    !Calculate electron tetrahedra
    if(num%tetrahedra) then
       call print_message("Calculating electron mesh tetrahedra...")
       call form_tetrahedra_3d(self%nwv, self%wvmesh, self%tetra, self%tetracount, &
            self%tetramap, .true., self%indexlist)
       call fill_tetrahedra_3d(self%tetra, self%ens, self%tetra_evals)
    else
       call print_message("Calculating electron mesh triangles...")
       call form_triangles(self%nwv, self%wvmesh, self%triang, self%triangcount, &
            self%triangmap, .true., self%indexlist)
       call fill_triangles(self%triang, self%ens, self%triang_evals)
    end if
  end subroutine calculate_electrons

  subroutine apply_energy_window(nk, indexlist, energies, enref, fsthick)
    !! Subroutine to find the Fermi window restricted blocks of BZ.
    !! This could be used for FBZ and IBZ.
    !!
    !! nk is the number of mesh points - will be updated to
    !! the number of mesh points in energy restricted blocks
    !!
    !! indexlist is the list of wave vector indices - will be updated
    !! to the list of indices in the blocks

    integer(k8), intent(inout) :: nk
    integer(k8), allocatable, intent(inout) :: indexlist(:)
    real(dp), intent(in) :: energies(:,:), enref, fsthick

    integer(k8) :: ik, count, numbands, inwindow(nk)
    real(dp), allocatable :: aux(:)
    
    numbands = size(energies(1,:))    
    allocate(aux(numbands))
    
    count = 0
    do ik = 1, nk
       aux = energies(ik, :)
       !Check if any band energy is within the Fermi window
       if(any(abs(aux(:) - enref) <= fsthick)) then
          count = count + 1
          inwindow(count) = ik !save index of in-window points
       end if
    end do

    if(count == 0) call exit_with_message("No states found within Fermi window.")
    
    !Update index list
    deallocate(indexlist)
    allocate(indexlist(count))
    indexlist(1:count) = inwindow(1:count)

    !Update number of irreducible points
    nk = count
  end subroutine apply_energy_window

  subroutine fbz_blocks_quantities(indexlist, energies, velocities)
    !! Subroutine to find FBZ quanties the lie within the Fermi window.

    integer(k8), intent(in) :: indexlist(:)
    real(dp), allocatable, intent(inout) :: energies(:,:), velocities(:,:,:)
    integer(k8) :: i, nk, numbands
    real(dp), allocatable :: energies_tmp(:,:), velocities_tmp(:,:,:)

    nk = size(indexlist)
    numbands = size(energies(1,:))

    allocate(energies_tmp(nk, numbands), velocities_tmp(nk, numbands, 3))
    
    do i = 1, nk
       energies_tmp(i,:) = energies(indexlist(i),:)
       velocities_tmp(i,:,:) = velocities(indexlist(i),:,:)
    end do

    deallocate(energies, velocities)
    allocate(energies(nk, numbands), velocities(nk, numbands, 3))

    energies(1:nk, :) = energies_tmp(1:nk, :)
    velocities(1:nk, :, :) = velocities_tmp(1:nk, :, :)
  end subroutine fbz_blocks_quantities

  subroutine calculate_carrier_conc(self, T, vol, h)
    !! Subroutine to calculate the band resolved carrier concentration
    !! for a given chemical potential and temperature.

    class(electron), intent(inout) :: self
    real(dp), intent(in) :: T, vol
    real(dp), intent(in), optional :: h

    !Local variables
    real(dp) :: const
    integer(k8) :: ib, ik

    !Allocate conc
    allocate(self%conc(self%numbands))
    self%conc = 0.0_dp

    self%conc_el = 0.0_dp
    self%conc_hole = 0.0_dp
    
    !Normalization and units factor
    const = self%spindeg/dble(product(self%wvmesh))/vol/(1.0e-21_dp)

    do ik = 1, self%nwv
       !Electron concentration
       !By convention, the electron carrier concentration will have a negative sign.
       if(self%indlowconduction > 0) then !Calculation includes conduction bands
          do ib = self%indlowconduction, self%indhighband !Conduction bands manifold
             self%conc(ib) = self%conc(ib) - Fermi(self%ens(ik, ib), self%chempot, T)
          end do
          !Total electron concentration
          self%conc_el = abs(sum(self%conc(self%indlowconduction:self%indhighband)))
       end if

       !Hole concentration
       !By convention, the hole carrier concentration will have a positive sign.
       if(self%indhighvalence > 0) then !Calculation includes valence bands
          do ib = self%indlowband, self%indhighvalence !Valence bands manifold
             self%conc(ib) = self%conc(ib) + (1.0_dp - Fermi(self%ens(ik, ib), self%chempot, T))
          end do
          !Total hole concentration
          self%conc_hole = sum(self%conc(self%indlowband:self%indhighvalence))
       end if
    end do
    self%conc = self%conc*const !cm^-3
    self%conc_el = self%conc_el*const !cm^-3
    self%conc_hole = self%conc_hole*const !cm^-3

    !If h is present that means the system is 2d
    if(present(h)) then
       self%conc = self%conc*h*1.0e-7_dp !cm^-2
       self%conc_el = self%conc_el*h*1.0e-7_dp !cm^-2
       self%conc_hole = self%conc_hole*h*1.0e-7_dp !cm^-2
    end if    
  end subroutine calculate_carrier_conc

  subroutine calculate_chempot(self, vol, dopingtype, Tlist, conclist, h)
    !! Subroutine to calculate the chemical potential for a
    !! given carrier concentration.

    class(electron), intent(in) :: self
    real(dp), intent(in) :: vol, Tlist(:), conclist(:)
    character(len = 1), intent(in) :: dopingtype
    real(dp), intent(in), optional :: h

    !Local variables
    integer(k8) :: ib, ik, it, ngrid, maxiter, itemp, iconc, &
         high, low, numtemp, numconc
    real(dp) :: a, b, aux, const, mu, thresh
    real(dp), allocatable :: chempot(:,:)

    call print_message("Calculating chemical potential...")

    !Number of temperature points
    numtemp = size(Tlist)

    !Number of cocentration points
    numconc = size(conclist)

    !Allocate chemical potential array
    allocate(chempot(numtemp, numconc))
    chempot = -99.99_dp
    
    !Total number of points in full mesh
    ngrid = product(self%wvmesh)

    !Normalization and units factor
    const = self%spindeg/dble(ngrid)/vol/(1.0e-21_dp)
    if(present(h)) then !2d system
       const = const*h*1.0e-7_dp
    end if

    !Maximum number of iterations
    maxiter = 5000

    !Convergence threshold
    thresh = 1.0e-12_dp

    !Check doping type
    if(dopingtype == 'n') then
       low = self%indlowconduction
       high = self%indhighband
    else
       low =  self%indlowband
       high = self%indhighvalence
    end if

    !Loop over temperatures
    do itemp = 1, numtemp
       if(this_image() == 1) then
          write(*,"(A, F7.2, A)") 'Crystal temperature = ', Tlist(itemp), ' K:'
          if(present(h)) then !2d system
             write(*, "(A)") 'Carrier conc. [cm^-2]    Chemical potential [eV]'
          else
             write(*, "(A)") 'Carrier conc. [cm^-3]    Chemical potential [eV]'
          end if
       end if
       !Loop over concentrations
       do iconc = 1, numconc
          if(dopingtype == 'n') then
             a = self%enref - 12.0_dp !guess lower bound
             b = self%enref + 12.0_dp !guess upper bound
          else
             a = self%enref + 12.0_dp !guess lower bound
             b = self%enref - 12.0_dp !guess upper bound
          end if
          do it = 1, maxiter
             mu = 0.5_dp*(a + b)
             aux = 0.0_dp
             do ib = low, high
                do ik = 1, self%nwv
                   if(dopingtype == 'n') then
                      aux = aux + Fermi(self%ens(ik, ib), mu, Tlist(itemp))
                   else
                      aux = aux + 1.0_dp - Fermi(self%ens(ik, ib), mu, Tlist(itemp))
                   end if
                end do
             end do
             aux = aux*const !cm^-3 for 3d, cm^-2 for 2d
             if(abs(aux - conclist(iconc))/conclist(iconc) < thresh) then
                exit
             else if(aux < conclist(iconc)) then
                a = mu
             else
                b = mu
             end if
          end do
          chempot(itemp,iconc) = mu
          
          if(abs(aux - conclist(iconc))/conclist(iconc) > thresh) then
             call exit_with_message(&
                  "Could not converge to correct chemical potential. Exiting.")
          end if

          if(this_image() == 1) then
             write(*, "(1E16.8, A, 1E16.8)") conclist(iconc), '         ', chempot(itemp,iconc)
          end if
       end do
    end do
  end subroutine calculate_chempot

  subroutine deallocate_eigenvecs(self)
    !! Deallocate the electron eigenvectors

    class(electron), intent(inout) :: self

    deallocate(self%evecs, self%evecs_irred)
  end subroutine deallocate_eigenvecs
end module electron_module
