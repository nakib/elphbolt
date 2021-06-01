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
  use misc, only: exit_with_message, print_message, demux_state, sort, &
       binsearch, subtitle
  use numerics_module, only: numerics
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal, calculate_wavevectors_full
  use symmetry_module, only: symmetry, find_irred_wedge, create_fbz2ibz_map
  use delta, only: form_tetrahedra_3d, fill_tetrahedra_3d
  
  implicit none

  private
  public electron

  type electron
     !! Data and procedures related to the electronic properties.

     character(len = 2) :: prefix = 'el'
     !! Prefix idenitfying particle type.
     integer(k8) :: spindeg
     !! Spin degeneracy.
     integer(k8) :: numbands
     !! Total number of electronic bands (from DFT/Wannier).
     integer(k8) :: numtransbands
     !! Total number of transport active bands.
     integer(k8) :: indlowband
     !! Lowest transport band index.
     integer(k8) :: indhighband
     !! Highest transport band index.
     integer(k8) :: indhighvalence
     !! Highest valence band index.
     integer(k8) :: indlowconduction
     !! Lowest conduction band index.
     integer(k8) :: mesh_ref
     !! Electron mesh refinement factor compared to the phonon mesh.
     integer(k8) :: kmesh(3)
     !! Electron wave vector mesh.
     integer(k8) :: nk
     !! Number of fine electron wave vectors in the full Brillouin zone (FBZ).
     integer(k8) :: nk_irred
     !! Number of fine electron wave vectors in the irreducible wedge of Brillouin zone (IBZ).
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
     !! Fermi surface thickness in (eV).
     real(dp) :: chempot
     !! Chemical potential
     real(dp) :: conc
     !! Electron(-)/hole(+) carrier concentration.
     real(dp) :: chimp_conc
     !! Concentration of charged impurities.
     real(dp) :: Z
     !! Ionization number of dopant.
     real(dp), allocatable :: wavevecs(:,:)
     !! List of all electron wave vectors (crystal coordinates).
     real(dp), allocatable :: wavevecs_irred(:,:)
     !! List of irreducible electron wave vectors (crystal coordinates).
     integer(k8), allocatable :: indexlist(:)
     !! List of muxed indices of the FBZ wave vectors.
     integer(k8), allocatable :: indexlist_irred(:)
     !! List of muxed indices of the IBZ wedge.
     integer(k8), allocatable :: nequiv(:)
     !! List of the number of equivalent points for each IBZ point.
     integer(k8), allocatable :: ibz2fbz_map(:,:,:)
     !! Map from an IBZ electron point to its images.
     !! The third axis contains the pair (symmetry index, image).
     integer(k8), allocatable :: fbz2ibz_map(:)
     !! Map from an FBZ electron point to its IBZ wedge image.
     integer(k8), allocatable :: equiv_map(:,:)
     !! Map of equivalent points under rotations.
     !! Axis 1 runs over rotations.
     !! Axis 2 runs over wave vectors.
     real(dp), allocatable :: symmetrizers(:,:,:)
     !! Symmetrizers of wave vector dependent vectors.
     integer(k8), allocatable :: tetra(:,:)
     !! List of all the wave vector mesh tetrahedra vertices.
     !! First axis list tetraheda and the second axis list the vertices.
     integer(k8), allocatable :: tetracount(:)
     !! The number of tetrahedra in which a wave vector belongs.
     integer(k8), allocatable :: tetramap(:,:,:)
     !! Mapping from a wave vector to the (tetrahedron, vertex) where it belongs.
     real(dp), allocatable :: tetra_evals(:,:,:)
     !! Tetrahedra vertices filled with eigenvalues.
     real(dp), allocatable :: ens(:,:)
     !! List of electron energies on FBZ.
     real(dp), allocatable :: ens_irred(:,:)
     !! List of electron energies on IBZ.
     real(dp), allocatable :: vels(:,:,:)
     !! List of electron velocities on FBZ.
     real(dp), allocatable :: vels_irred(:,:,:)
     !! List of electron velocites on IBZ.
     complex(dp), allocatable :: evecs(:,:,:)
     !! List of all electron eigenvectors.
     complex(dp), allocatable :: evecs_irred(:,:,:)
     !! List of IBZ wedge electron eigenvectors.
     logical :: metallic
     !! Is the system metallic?
     !! Semimetals and semiconductors require providing conduction and/or valence bands.
     real(dp), allocatable :: dos(:,:)
     !! Band resolved density of states.
     
   contains

     procedure :: initialize=>read_input_and_setup, deallocate_eigenvecs
  end type electron

contains

  subroutine read_input_and_setup(el, wann, crys, sym, num)
    !! Read input file and setup groundstate electronic system.

    class(electron), intent(out) :: el
    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(numerics), intent(in) :: num

    !Local variables
    real(dp) :: enref, conc, chimp_conc
    integer(k8) :: spindeg, numbands, numtransbands, indlowband, indhighband, &
         indhighvalence, indlowconduction
    logical :: metallic

    namelist /electrons/ enref, spindeg, numbands, numtransbands, &
         indlowband, indhighband, metallic, indhighvalence, indlowconduction, &
         conc, chimp_conc

    call subtitle("Setting up electrons...")
    
    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Read electrons information
    spindeg = 2 !Default calculation is non-spin polarized
    numbands = -1
    numtransbands = -1
    indlowband = -1
    indhighband = -1
    indhighvalence = -1
    indlowconduction = -1
    metallic = .false.
    conc = 0.0_dp
    chimp_conc = 0.0_dp
    read(1, nml = electrons)
    if(spindeg < 1 .or. spindeg > 2) then
       call exit_with_message('spindeg can be 1 or 2.')
    end if
    if(numbands < 1) then
       call exit_with_message('numbands should be > 0.')
    end if
    if(numtransbands < 1) then
       call exit_with_message('numtransbands should be > 0.')
    end if
    if(indlowband < 1) then
       call exit_with_message('indlowband should be > 0.')
    end if
    if(indhighband < 1) then
       call exit_with_message('indhighband should be > 0.')
    end if
    if(.not. metallic) then
       if(indhighvalence < 1 .and. indlowconduction < 1) then
          call exit_with_message('Either indhighvalence or indlowconduction should be >= 1.')
       end if
       if(indhighvalence >= 1) then
          if(indhighvalence < indlowband) then
             call exit_with_message('indhighvalence should be >= indlowband.')
          end if
       end if
       if(indlowconduction >= 1) then
          if(indlowconduction > indhighband) then
             call exit_with_message('indlowconduction should be <= indhighband.')
          end if
       end if
    end if
    el%spindeg = spindeg
    el%numbands = numbands
    el%numtransbands = numtransbands
    el%indlowband = indlowband
    el%indhighband = indhighband
    el%indhighvalence = indhighvalence
    el%indlowconduction = indlowconduction
    el%metallic = metallic
    el%enref = enref
    el%conc = conc
    el%chimp_conc = chimp_conc
    if(el%chimp_conc /= 0.0_dp) then
       el%Z = abs(el%conc/el%chimp_conc)
    else
       el%Z = 0.0_dp
    end if
    
    !Close input file
    close(1)

    !Set some electronic properties from the numerics object
    el%mesh_ref = num%mesh_ref
    el%kmesh = el%mesh_ref*num%qmesh
    el%fsthick = num%fsthick

    !Print out information.
    if(this_image() == 1) then
       write(*, "(A, I1)") "Spin degeneracy = ", el%spindeg
       write(*, "(A, I5)") "Number of electronic bands = ", el%numbands
       write(*, "(A, 1E16.8)") "Reference electron energy = ", el%enref
       write(*, "(A, 1E16.8)") "Specified carrier concentration = ", el%conc
       write(*, "(A, 1E16.8)") "Charged impurity concentration = ", el%chimp_conc
       write(*, "(A, 1E16.8)") "Ionization of impurity = ", el%Z
    end if
    
    !Calculate electrons
    call calculate_electrons(el, wann, crys, sym, num)
  end subroutine read_input_and_setup
  
  subroutine calculate_electrons(el, wann, crys, sym, num)
    !! Calculate electron energy window restricted wave vector meshes
    !! and the electronic properties on them

    class(electron), intent(inout) :: el
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
    el%nk = product(el%kmesh)
    
    !The electronic mesh setup proceeds in multiple steps:
    ! 1. Calculate full electron wave vector mesh
    call print_message("Calculating FBZ...")
    blocks = .false.
    call calculate_wavevectors_full(el%kmesh, el%wavevecs, blocks)
    
    ! 2. Calculate the IBZ
    call print_message("Calculating IBZ and IBZ -> FBZ mappings...")
    call find_irred_wedge(el%kmesh, el%nk_irred, el%wavevecs_irred, &
         el%indexlist_irred, el%nequiv, sym%nsymm_rot, sym%qrotations, &
         el%ibz2fbz_map, el%equiv_map, blocks)
    
    ! 3. Calculate IBZ quantities
    call print_message("Calculating IBZ energies...")
    allocate(el%ens_irred(el%nk_irred, wann%numwannbands), &
         el%vels_irred(el%nk_irred, wann%numwannbands, 3), &
         el%evecs_irred(el%nk_irred, wann%numwannbands, wann%numwannbands))
    call wann%el_wann_epw(crys, el%nk_irred, el%wavevecs_irred, el%ens_irred, &
         el%vels_irred, el%evecs_irred)
    
    ! 4. Map out FBZ quantities from IBZ ones
    call print_message("Mapping out FBZ energies...")
    allocate(el%indexlist(el%nk), el%ens(el%nk, wann%numwannbands), el%vels(el%nk, wann%numwannbands, 3))
    
    do i = 1,el%nk_irred !an irreducible point
       do l = 1,el%nequiv(i) !number of equivalent points of i
          il = el%ibz2fbz_map(l, i, 2) ! (i, l) -> il
          s = el%ibz2fbz_map(l, i, 1) ! mapping rotation

          !index list
          el%indexlist(il) = il

          !energy
          el%ens(il,:) = el%ens_irred(i,:)
          
          !velocity
          do ib = 1,wann%numwannbands
             !here use real space (Cartesian) rotations
             el%vels(il, ib, :) = matmul(sym%crotations(:, :, s), el%vels_irred(i, ib, :))
          end do
       end do
    end do
    
    call print_message("Transport energy window restricted calculation:")
    call print_message("-----------------------------------------------")
    
    ! 5. Find energy window restricted FBZ blocks.
    !    After this step, el%nk, el%indexlist will refer
    !    to the energy restricted mesh.
    call print_message("Calculating Fermi window restricted FBZ blocks...")
    call apply_energy_window(el%nk, el%indexlist, el%ens, el%enref, el%fsthick)
    
    ! 6. Sort index list and related quanties of FBZ blocks
    call print_message("Sorting FBZ blocks index list...")
    call sort(el%indexlist)

    ! 7. Get FBZ blocks wave vectors, energies, velocities and eigenvectors.
    !    After this step, el%wavevecs, el%ens, el%vels, and el%evecs
    !    will refer to the energy restricted mesh.
    call print_message("Calcutating FBZ blocks quantities...")
    
    !wave vectors
    deallocate(el%wavevecs)
    
    blocks = .true.
    call calculate_wavevectors_full(el%kmesh, el%wavevecs, blocks, el%indexlist) !wave vectors
    
    !energies and velocities
    call fbz_blocks_quantities(el%indexlist, el%ens, el%vels)

    !Get FBZ blocks eigenvectors from direct calculations since we are
    !not getting these from IBZ quantities via symmetry rotations
    allocate(el%evecs(el%nk, wann%numwannbands, wann%numwannbands))
    allocate(el_ens_tmp(el%nk, wann%numwannbands), el_vels_tmp(el%nk, wann%numwannbands, 3))
    call wann%el_wann_epw(crys, el%nk, el%wavevecs, el_ens_tmp, el_vels_tmp, el%evecs)
    deallocate(el_ens_tmp, el_vels_tmp) !free up memory

    ! 8. Find IBZ of energy window restricted blocks
    !    After this step, el%nk_irred, el%indexlist_irred, 
    !    el%wavevecs_irred, el%nequiv, and el%ibz2fbz_map
    !    will refer to the energy restricted mesh 
    call print_message("Calculating IBZ blocks...")
    deallocate(el%wavevecs_irred, el%indexlist_irred, el%nequiv, &
         el%ibz2fbz_map, el%equiv_map)
    blocks = .true.
    call find_irred_wedge(el%kmesh, el%nk_irred, el%wavevecs_irred, &
         el%indexlist_irred, el%nequiv, sym%nsymm_rot, sym%qrotations, &
         el%ibz2fbz_map, el%equiv_map, blocks, el%indexlist)

    !Create symmetrizers of wave vector dependent vectors ShengBTE style
    allocate(el%symmetrizers(3, 3, el%nk))
    el%symmetrizers = 0.0_dp
    do i = 1, el%nk
       ii = el%indexlist(i)
       kk = 0
       do jj = 1, sym%nsymm
          if(el%equiv_map(jj, i) == ii) then
             el%symmetrizers(:, :, i) = el%symmetrizers(:, :, i) + &
                  sym%crotations_orig(:, :, jj)
             kk = kk + 1
          end if
       end do
       if(kk > 1) then
          el%symmetrizers(:, :, i) = el%symmetrizers(:, :, i)/kk
       end if
    end do
    
    ! 9. Get IBZ blocks energies, velocities, and eigen vectors.
    call print_message("Calcutating IBZ blocks quantities...")
    deallocate(el%ens_irred, el%vels_irred, el%evecs_irred)
    allocate(el%ens_irred(el%nk_irred, wann%numwannbands), &
         el%vels_irred(el%nk_irred, wann%numwannbands, 3), &
         el%evecs_irred(el%nk_irred, wann%numwannbands, wann%numwannbands))
    call wann%el_wann_epw(crys, el%nk_irred, el%wavevecs_irred, el%ens_irred, &
         el%vels_irred, el%evecs_irred)
    
    ! 10. Calculate the number of FBZ blocks electronic states
    !     available for scattering
    el%nstates_inwindow = 0
    do i = 1,el%nk !over FBZ blocks
       do ib = 1,wann%numwannbands !bands
          if(abs(el%ens(i, ib) - el%enref) <= el%fsthick) &
               el%nstates_inwindow = el%nstates_inwindow + 1
       end do
    end do
    if(this_image() == 1) write(*, "(A, I10)") &
         " Number of energy restricted FBZ blocks states = ", el%nstates_inwindow

    ! 11. Create FBZ blocks to IBZ blocks map
    call print_message("Calculating FBZ -> IBZ mappings...")
    !call create_fbz2ibz_map
    call create_fbz2ibz_map(el%fbz2ibz_map,el%nk,el%nk_irred,el%indexlist,el%nequiv,el%ibz2fbz_map)
    
    do i = 1, el%nk_irred !IBZ
       do l = 1, el%nequiv(i) !number of equivalent points of i
          il = el%ibz2fbz_map(l, i, 2) ! (i, l) -> il
          s = el%ibz2fbz_map(l, i, 1) ! symmetry
          call binsearch(el%indexlist, il, aux)

          !energy
          el%ens(aux,:) = el%ens_irred(i,:)
          
          !velocity
          do ib = 1,wann%numwannbands
             !here use real space (Cartesian) rotations
             el%vels(aux, ib, :) = matmul(sym%crotations(:, :, s), el%vels_irred(i, ib, :))
          end do
          el%vels(aux,:,:) = transpose(&
               matmul(el%symmetrizers(:,:,aux),transpose(el%vels(aux,:,:))))
       end do
    end do
        
    ! 12. Calculate the number of IBZ electronic states available for scattering
    el%nstates_irred_inwindow = 0
    do istate = 1,el%nk_irred*wann%numwannbands
       !Demux state index into band (ib) and wave vector (i) indices
       call demux_state(istate, wann%numwannbands, ib, i)
       if(abs(el%ens_irred(i, ib) - el%enref) <= el%fsthick) then 
          el%nstates_irred_inwindow = el%nstates_irred_inwindow + 1
       end if
    end do
    if(this_image() == 1) write(*, "(A, I10)") " Number of energy restricted IBZ blocks states = ", &
         el%nstates_irred_inwindow
    
    !Calculate list of IBZ in-window states = (wave vector index, band index)
    allocate(el%IBZ_inwindow_states(el%nstates_irred_inwindow,2))
    count = 0
    do istate = 1, el%nk_irred*wann%numwannbands
       !Demux state index into band (ib) and wave vector (i) indices
       call demux_state(istate, wann%numwannbands, ib, i)
       if(abs(el%ens_irred(i, ib) - el%enref) <= el%fsthick) then
          count = count + 1
          el%IBZ_inwindow_states(count,:) = (/i, ib/)
       end if
    end do

    !Write IBZ in-window states as text data to file
    if(this_image() == 1) then
       call chdir(num%cwd)
       filename = 'IBZ_inwindow_states'
       write(numcols,"(I0)") 2
       open(1,file=trim(filename),status='replace')
       write(1,*) "#k-vec index     band index"
       do i = 1, el%nstates_irred_inwindow
          write(1,"("//trim(adjustl(numcols))//"I10)") el%IBZ_inwindow_states(i,:)
       end do
       close(1)
    end if

    !Print out irreducible electron energies
    write(numcols, "(I0)") el%numbands
    if(this_image() == 1) then
       open(1, file = "el.ens", status = "replace")
       do i = 1, el%nk_irred
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               el%ens_irred(i, :)
       end do
       close(1)
    end if

    !Calculate electron tetrahedra
    if(num%tetrahedra) then
       call print_message("Calculating electron mesh tetrahedra...")
       call form_tetrahedra_3d(el%nk, el%kmesh, el%tetra, el%tetracount, &
            el%tetramap, .true., el%indexlist)
       call fill_tetrahedra_3d(el%tetra, el%ens, el%tetra_evals)
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

  subroutine deallocate_eigenvecs(el)
    !! Deallocate the electron eigenvectors

    class(electron), intent(inout) :: el

    deallocate(el%evecs, el%evecs_irred)
  end subroutine deallocate_eigenvecs
end module electron_module
