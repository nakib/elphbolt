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

module phonon_module
  !! Module containing type and procedures related to the phononic properties.

  use params, only: dp, k8, bohr2nm, pi, twopi, Ryd2eV, oneI
  use misc, only: print_message, subtitle, expi, distribute_points, &
       write2file_rank2_real, exit_with_message
  use numerics_module, only: numerics
  use wannier_module, only: epw_wannier
  use crystal_module, only: crystal, calculate_wavevectors_full
  use symmetry_module, only: symmetry, find_irred_wedge, create_fbz2ibz_map
  use delta, only: form_tetrahedra_3d, fill_tetrahedra_3d, form_triangles, &
       fill_triangles
  
  implicit none

  private
  public phonon, phonon_espresso
  
  type phonon
     !! Data and procedures related to phonons.
     
     character(len = 2) :: prefix = 'ph'
     !! Prefix idenitfying particle type.
     integer(k8) :: numbranches
     !! Total number of phonon branches.
     integer(k8) :: nq
     !! Number of phonon wave vectors in the full Brillouin zone (FBZ).
     integer(k8) :: nq_irred
     !! Number of phonon wave vectors in the irreducible wedge of Brillouin zone (IBZ).
     integer(k8) :: qmesh(3) 
     !! Phonon wave vector mesh.
     real(dp), allocatable :: wavevecs(:,:)
     !! List of all phonon wave vectors (crystal coordinates).
     real(dp), allocatable :: wavevecs_irred(:,:)
     !! List of irreducible phonon wave vectors (crystal coordinates).
     integer(k8), allocatable :: indexlist_irred(:)
     !! List of muxed indices of the IBZ wedge.
     integer(k8), allocatable :: nequiv(:)
     !! List of the number of equivalent points for each IBZ point.
     integer(k8), allocatable :: ibz2fbz_map(:,:,:)
     !! Map from an IBZ phonon point to its images.
     !! The third axis contains the pair (symmetry index, image).
     integer(k8), allocatable :: equiv_map(:,:)
     !! Map of equivalent points under rotations.
     !! Axis 1 runs over rotations.
     !! Axis 2 runs over wave vectors (full Brillouin zone).
     real(dp), allocatable :: symmetrizers(:,:,:)
     !! Symmetrizers of wave vector dependent vectors.
     integer(k8), allocatable :: tetra(:,:)
     !! List of all the wave vector mesh tetrahedra vertices.
     !! First axis lists tetraheda and the second axis lists the vertices.
     integer(k8), allocatable :: tetracount(:)
     !! The number of tetrahedra in which a wave vector belongs.
     integer(k8), allocatable :: tetramap(:,:,:)
     !! Mapping from a wave vector to the (tetrahedron, vertex) where it belongs.
     real(dp), allocatable :: tetra_evals(:,:,:)
     !! Tetrahedra vertices filled with eigenvalues.
     integer(k8), allocatable :: triang(:,:)
     !! List of all the wave vector mesh triangles vertices.
     !! First axis lists triangles and the second axis lists the vertices.
     integer(k8), allocatable :: triangcount(:)
     !! The number of triangles in which a wave vector belongs.
     integer(k8), allocatable :: triangmap(:,:,:)
     !! Mapping from a wave vector to the (triangle, vertex) where it belongs.
     real(dp), allocatable :: triang_evals(:,:,:)
     !! Triangles vertices filled with eigenvalues.
     real(dp), allocatable :: ens(:,:)
     !! List of phonon energies on FBZ.
     real(dp), allocatable :: vels(:,:,:)
     !! List of phonon velocities on IBZ.
     complex(dp), allocatable :: evecs(:,:,:)
     !! List of all phonon eigenvectors on IBZ.
     integer(k8) :: scell(3)
     !! q-mesh used in DFPT or, equivalently, supercell used in finite displencement
     !! method for calculating the 2nd order force constants.
     real(dp), allocatable :: ifc2(:,:,:,:,:,:,:)
     !! Second order force constants (ifc2) tensor.
     real(dp), allocatable :: ifc3(:,:,:,:)
     !! Third order force constants (ifc3) tensor.
     integer(k8) :: numtriplets
     !! Number of triplets in the ifc3 file.
     real(dp), allocatable :: R_j(:,:), R_k(:,:)
     !! Position of the 2nd and 3rd atoms in supercell for an ifc3 triplet.
     integer(k8), allocatable :: Index_i(:), Index_j(:), Index_k(:)
     !! Label of primitive cell atoms in the ifc3 triplet.
     real(dp), allocatable :: dos(:,:)
     !! Branch resolved density of states.

     !Data read from ifc2 file. This will be used in the phonon calculation.
     real(dp) :: rws(124, 0:3), cell_r(1:3, 0:3), cell_g(1:3, 0:3)
     real(dp), allocatable :: mm(:,:), rr(:,:,:)
      
   contains

     procedure :: initialize, deallocate_phonon_quantities
     
  end type phonon

contains

  subroutine initialize(ph, wann, crys, sym, num)
    !! Initialize the phonon data type, calculate ground state phonon properties,
    !! and read 3rd order force constants data. 

    class(phonon), intent(out) :: ph
    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(numerics), intent(in) :: num

    call subtitle("Setting up phonons...")

    !Set phonon branches
    ph%numbranches = crys%numatoms*3
    !Set wave vector mesh
    ph%qmesh = num%qmesh
    !Set number of phonon wave vectors
    ph%nq = product(ph%qmesh(:))

    !Read ifc2 and related quantities
    call read_ifc2(ph, crys)
    
    !Calculate harmonic properties
    call calculate_phonons(ph, wann, crys, sym, num)

    if(.not. num%onlyebte) then
       !Read ifc3s and related quantities
       call read_ifc3(ph, crys)
    end if
  end subroutine initialize

  subroutine deallocate_phonon_quantities(ph)
    !! Deallocate the electron eigenvectors

    class(phonon), intent(inout) :: ph

    deallocate(ph%evecs, ph%ifc2, ph%ifc3, ph%Index_i, ph%Index_j, ph%Index_k, &
         ph%mm, ph%rr)
  end subroutine deallocate_phonon_quantities
  
  subroutine calculate_phonons(ph, wann, crys, sym, num)
    !! Calculate phonon quantities on the FBZ and IBZ meshes.

    class(phonon), intent(inout) :: ph
    type(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(numerics), intent(in) :: num
    
    !Local variables
    integer(k8) :: i, iq, ii, jj, kk, l, il, s, ib, im, chunk, &
         num_active_images
    integer(k8), allocatable :: start[:], end[:]
    real(dp), allocatable :: ens_chunk(:,:)[:], vels_chunk(:,:,:)[:], &
         symmetrizers_chunk(:,:,:)[:]
    complex(dp), allocatable :: evecs_chunk(:,:,:)[:]
    !Switch for mesh utilites with or without energy restriction
    logical :: blocks
    character(len = 1024) :: numcols

    blocks = .false.

    call print_message("Calculating phonon FBZ quantities...")

    !Allocate start and end coarrays
    allocate(start[*], end[*])

    !Divide wave vectors among images
    call distribute_points(ph%nq, chunk, start, end, num_active_images)

    !Allocate small work variable chunk for each image
    allocate(ens_chunk(chunk, ph%numbranches)[*])
    allocate(vels_chunk(chunk, ph%numbranches, 3)[*])
    allocate(evecs_chunk(chunk, ph%numbranches, ph%numbranches)[*])

    !Calculate FBZ mesh
    call calculate_wavevectors_full(ph%qmesh, ph%wavevecs, blocks)

    !Print phonon FBZ mesh
    call write2file_rank2_real("ph.wavevecs_fbz", ph%wavevecs)
    
    !Calculate FBZ phonon quantities
    call phonon_espresso(ph, crys, chunk, ph%wavevecs(start:end, :), &
         ens_chunk, evecs_chunk, vels_chunk)

    !Gather the chunks from the images
    allocate(ph%ens(ph%nq, ph%numbranches))
    allocate(ph%vels(ph%nq, ph%numbranches, 3))
    allocate(ph%evecs(ph%nq, ph%numbranches, ph%numbranches))
    sync all
    do im = 1, num_active_images
       ph%ens(start[im]:end[im], :) = ens_chunk(:,:)[im]
       ph%vels(start[im]:end[im], :, :) = vels_chunk(:,:,:)[im]
       ph%evecs(start[im]:end[im], :, :) = evecs_chunk(:,:,:)[im]
    end do
    sync all
    deallocate(ens_chunk, vels_chunk, evecs_chunk)
    
    !Calculate IBZ mesh
    call print_message("Calculating IBZ and IBZ -> FBZ mappings...")
    call find_irred_wedge(ph%qmesh, ph%nq_irred, ph%wavevecs_irred, &
         ph%indexlist_irred, ph%nequiv, sym%nsymm_rot, sym%qrotations, &
         ph%ibz2fbz_map, ph%equiv_map, blocks)

    !Print phonon IBZ mesh
    call write2file_rank2_real("ph.wavevecs_ibz", ph%wavevecs_irred)
    
    !Create symmetrizers of wave vector dependent vectors ShengBTE style
    allocate(symmetrizers_chunk(3, 3, chunk)[*])
    symmetrizers_chunk = 0.0_dp
    do iq = start, end
       kk = 0
       do jj = 1, sym%nsymm
          if(ph%equiv_map(jj, iq) == iq) then
             symmetrizers_chunk(:, :, iq - start + 1) = &
                  symmetrizers_chunk(:, :, iq - start + 1) + &
                  sym%crotations_orig(:, :, jj)
             kk = kk + 1
          end if
       end do
       if(kk > 1) then
          symmetrizers_chunk(:, :, iq - start + 1) = &
               symmetrizers_chunk(:, :, iq - start + 1)/kk
       end if
    end do
    
    allocate(ph%symmetrizers(3, 3, ph%nq))
    sync all
    do im = 1, num_active_images
       ph%symmetrizers(:, :, start[im]:end[im]) = symmetrizers_chunk(:,:,:)[im]
    end do
    sync all
    deallocate(symmetrizers_chunk)
    
    !Symmetrize phonon energies and velocities.
    do i = 1, ph%nq_irred !an irreducible point
       ii = ph%indexlist_irred(i)
       ph%vels(ii,:,:)=transpose(&
            matmul(ph%symmetrizers(:,:,ii),transpose(ph%vels(ii,:,:))))
       do l = 1, ph%nequiv(i) !number of equivalent points of i
          il = ph%ibz2fbz_map(l, i, 2) ! (i, l) -> il
          s = ph%ibz2fbz_map(l, i, 1) ! mapping rotation

          !energy
          ph%ens(il,:) = ph%ens(ii,:)

          !velocity
          do ib = 1, ph%numbranches
             !here use real space (Cartesian) rotations
             ph%vels(il, ib, :) = matmul(sym%crotations(:, :, s), ph%vels(ii, ib, :))
          end do
       end do
    end do
    
    !Print out irreducible phonon energies and velocities
    if(this_image() == 1) then
       write(numcols, "(I0)") ph%numbranches
       open(1, file = "ph.ens_ibz", status = "replace")
       do iq = 1, ph%nq_irred
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               ph%ens(ph%indexlist_irred(iq), :)
       end do
       close(1)

       write(numcols, "(I0)") 3*ph%numbranches
       open(1, file = "ph.vels_ibz", status = "replace")
       do iq = 1, ph%nq_irred
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               ph%vels(ph%indexlist_irred(iq), :, :)
       end do
       close(1)
    end if
    
    !Calculate phonon tetrahedra
    if(num%tetrahedra) then
       call print_message("Calculating phonon mesh tetrahedra...")
       call form_tetrahedra_3d(ph%nq, ph%qmesh, ph%tetra, ph%tetracount, &
            ph%tetramap, .false.)
       call fill_tetrahedra_3d(ph%tetra, ph%ens, ph%tetra_evals)
    else
       call print_message("Calculating phonon mesh triangles...")
       call form_triangles(ph%nq, ph%qmesh, ph%triang, ph%triangcount, &
            ph%triangmap, .false.)
       call fill_triangles(ph%triang, ph%ens, ph%triang_evals)
    end if
  end subroutine calculate_phonons
  
  subroutine read_ifc2(ph, crys)
    !! Read the 2nd order force constants from the Quantum Espresso format.
    !! This is adapted from ShengBTE.
    
    class(phonon), intent(inout) :: ph
    type(crystal), intent(in) :: crys

    !Local variables
    integer(k8) :: qscell(3), tipo(crys%numatoms), t1, t2, t3, i, j, &
         iat, jat, ibrav, ipol, jpol, m1, m2, m3, ntype, nat, nfc2 
    real(dp) :: r(crys%numatoms, 3), wscell(3,0:3), celldm(6), at(3,3), &
         mass(crys%numelements), zeff(crys%numatoms, 3, 3), eps(3, 3), &
         dnrm2
    character(len = 1) :: polar_key
    character(len = 6) :: label(crys%numelements)
    real(dp), parameter :: massfactor=1.8218779_dp*6.022e-4_dp

    allocate(ph%mm(crys%numatoms, crys%numatoms))
    allocate(ph%rr(crys%numatoms, crys%numatoms, 3))
    
    open(1,file="espresso.ifc2",status="old")
    !Read some stuff that will not be used in the code.
    read(1,*) ntype, nat, ibrav, celldm(1:6)
    if (ibrav==0) then
       read(1,*) ((at(i,j),i=1,3),j=1,3)
    end if

    do i = 1, ntype
       read(1, *) j, label(i), mass(i)
    end do
    mass = crys%masses/massfactor
    
    do i = 1, nat
       read(1, *) j, tipo(i), r(i, 1:3)
    end do
    r = transpose(matmul(crys%lattvecs, crys%basis))/bohr2nm
    
    read(1, *) polar_key
    if(polar_key == "T") then
       do i = 1, 3
          read(1, *) eps(i, 1:3)
       end do
       do i = 1, nat
          read(1, *)
          do j = 1, 3
             read(1, *) zeff(i, j, 1:3)
          end do
       end do
    end if
    read(1,*) qscell(1:3)

    ph%scell = qscell
    
    !Read the force constants.
    allocate(ph%ifc2(3, 3, nat, nat, ph%scell(1), ph%scell(2), ph%scell(3)))
    nfc2 = 3*3*nat*nat
    do i = 1, nfc2
       read(1, *) ipol, jpol, iat, jat
       do j = 1, ph%scell(1)*ph%scell(2)*ph%scell(3)
          read(1, *) t1, t2, t3, &
               ph%ifc2(ipol, jpol, iat, jat, t1, t2, t3)
       end do
    end do
    close(1)
    
    !Enforce the conservation of momentum in the simplest way possible.
    do i = 1, 3
       do j = 1, 3
          do iat = 1, nat
             ph%ifc2(i, j, iat, iat, 1, 1, 1) = ph%ifc2(i, j, iat, iat, 1, 1, 1) - &
                  sum(ph%ifc2(i, j, iat, :, :, :, :))
          end do
       end do
    end do

    ph%cell_r(:, 1:3) = transpose(crys%lattvecs)/bohr2nm
    do i = 1, 3
       ph%cell_r(i, 0) = dnrm2(3, ph%cell_r(i, 1:3), 1)
    end do
    ph%cell_g(:, 1:3) = transpose(crys%reclattvecs)*bohr2nm
    do i = 1, 3
       ph%cell_g(i, 0) = dnrm2(3, ph%cell_g(i, 1:3), 1)
    end do

    wscell(1,1:3) = ph%cell_r(1,1:3)*ph%scell(1)
    wscell(2,1:3) = ph%cell_r(2,1:3)*ph%scell(2)
    wscell(3,1:3) = ph%cell_r(3,1:3)*ph%scell(3)

    j = 1
    do m1 = -2, 2
       do m2 = -2, 2
          do m3 = -2, 2
             if(all((/m1, m2, m3/).eq.0)) then
                cycle
             end if
             do i = 1, 3
                ph%rws(j, i) = wscell(1, i)*m1 + wscell(2, i)*m2 + wscell(3, i)*m3
             end do
             ph%rws(j, 0) = 0.5*dot_product(ph%rws(j, 1:3), ph%rws(j, 1:3))
             j = j + 1
          end do
       end do
    end do

    do i = 1, nat
       ph%mm(i, i) = mass(tipo(i))
       ph%rr(i, i, :) = 0
       do j = i + 1, nat
          ph%mm(i, j) = sqrt(mass(tipo(i))*mass(tipo(j)))
          ph%rr(i, j, 1:3) = r(i, 1:3) - r(j, 1:3)
          ph%mm(j, i) = ph%mm(i, j)
          ph%rr(j, i, 1:3) = -ph%rr(i, j, 1:3)
       end do
    end do
    
  end subroutine read_ifc2
  
  subroutine read_ifc3(ph, crys)
    !! Read the 3rd order force constants in the thirdorder.py format.
    !! This subroutine is adapted from ShengBTE.

    class(phonon), intent(inout) :: ph
    type(crystal), intent(in) :: crys
    
    !Local variables
    real(dp) :: tmp(3,3), r(crys%numatoms, 3), wscell(3,0:3), celldm(6), at(3,3), &
         mass(crys%numelements), zeff(crys%numatoms, 3, 3), eps(3, 3)
    real(dp), allocatable :: fc(:,:,:,:)
    integer(k8) :: ii, jj, ll, mm, nn, ltem, mtem, ntem, info, P(3), &
         na1, na2, na3, j1, j2, j3, na1_, na2_, na3_, j1_, j2_, j3_, sc_nat, &
         triplet_counter, nR, nR_, qscell(3), tipo(crys%numatoms), t1, t2, t3, i, j, &
         iat, jat, ibrav, m1, m2, m3, ntype, nat, jn1, jn2, jn3
    integer(k8), allocatable :: R2(:,:), R3(:,:)
    character(len = 1) :: polar_key
    character(len = 6) :: label(crys%numelements)
    logical :: sheng_file_exists, d3q_file_exists, save_nR

    !Check what force constants files have been provided
    sheng_file_exists = .False.
    d3q_file_exists = .False.
    inquire(file = 'FORCE_CONSTANTS_3RD', exist = sheng_file_exists)
    inquire(file = 'mat3R', exist = d3q_file_exists)

    if(sheng_file_exists .and. .not. d3q_file_exists) then
       call print_message('ShengBTE format third order force constants provided.')
    else if(d3q_file_exists .and. .not. sheng_file_exists) then
       call print_message('d3q format third order force constants provided.')
    else if(d3q_file_exists .and. sheng_file_exists) then
       call print_message(&
            'Both ShengBTE and d3q format third order force constants provided. Defaulting to ShengBTE format.')
       d3q_file_exists = .False.
    else
       call exit_with_message('Third order force constant file not provided. Exiting.')
    end if

    if(sheng_file_exists) then
       !The file is in a simple sparse format, described in detail in
       !the user documentation. See Doc/ShengBTE.pdf.
       open(1, file = 'FORCE_CONSTANTS_3RD', status = "old")
       read(1, *) ph%numtriplets
       allocate(ph%Index_i(ph%numtriplets), ph%Index_j(ph%numtriplets), ph%Index_k(ph%numtriplets))
       allocate(ph%ifc3(3, 3, 3, ph%numtriplets), ph%R_j(3, ph%numtriplets), ph%R_k(3,ph%numtriplets))
       do ii = 1, ph%numtriplets
          read(1, *) jj
          read(1, *) ph%R_j(:, ii) !Ang
          read(1, *) ph%R_k(:, ii) !Ang
          read(1, *) ph%Index_i(ii), ph%Index_j(ii), ph%Index_k(ii)
          do ll = 1, 3
             do mm = 1, 3
                do nn = 1, 3
                   read(1, *) ltem, mtem, ntem, ph%ifc3(ll, mm, nn, ii)
                end do
             end do
          end do
       end do
       close(1)
       !IFC3 units are eV/Ang^3
    else
       !See SUBROUTINE read_fc3_grid of fc3_interp.f90 of the d3q code
       !for more information about the format.
       open(1, file = 'mat3R', status = "old")
       
       !Read some stuff that will not be used in the code.
       read(1,*) ntype, nat, ibrav, celldm(1:6)
       if (ibrav==0) then
          read(1,*) ((at(i,j),i=1,3),j=1,3)
       end if

       do i = 1, ntype
          read(1, *) j, label(i), mass(i)
       end do

       do i = 1, nat
          read(1, *) j, tipo(i), r(i, 1:3)
       end do

       read(1, *) polar_key
       if(polar_key == "T") then
          do i = 1, 3
             read(1, *) eps(i, 1:3)
          end do
          do i = 1, nat
             read(1, *)
             do j = 1, 3
                read(1, *) zeff(i, j, 1:3)
             end do
          end do
       end if
       read(1,*) qscell(1:3)

       !Number of atoms in the supercell
       sc_nat = product(qscell)*crys%numatoms

       save_nR = .true.
       do na1 = 1, sc_nat
          do na2 = 1, sc_nat
             do na3 = 1, sc_nat
                do j1 =1, 3
                   jn1 = j1 + (na1 - 1)*3
                   do j2 =1, 3
                      jn2 = j2 + (na2 - 1)*3
                      do j3 =1, 3
                         jn3 = j3 + (na3 - 1)*3
                         !Read tensor elements location and triplet atoms.
                         read(1, *) j1_, j2_, j3_, na1_, na2_, na3_
                         if(any((/na1, na2, na3, j1, j2, j3/) /= &
                              (/na1_, na2_, na3_, j1_, j2_, j3_/))) then
                            call exit_with_message(&
                                 "Wrong Triplet indices and/or tensor element location in mat3R file. Exiting.")
                         end if

                         !Read number of unit cells in file.
                         read(1, *) nR_
                         !Save this number only the first time it is read.
                         !Also, allocate the various quantities.
                         if(save_nR) then
                            nR = nR_
                            allocate(R2(3, nR), R3(3, nR), fc(3*sc_nat, 3*sc_nat, 3*sc_nat, nR))
                            save_nR = .false.
                         end if

                         if(nR_ /= nR) call exit_with_message(&
                              "Wrong number of unit cells in mat3R file. Exiting.")

                         do ii = 1, nR
                            !R2 and R3 seem to be the same for every nR chunk. Not
                            !sure why the file is written in this format.
                            read(1, *) R2(:, ii), R3(:, ii), fc(jn1, jn2, jn3, ii)
                            !At this point the fc units are Ry/Bohr^3
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do

       !Convert to the standard format.
       ph%numtriplets = nR*sc_nat**3 !Number of triplets
       triplet_counter = 0
       do ii = 1, nR
          do na1 = 1, sc_nat
             do na2 = 1, sc_nat
                do na3 = 1, sc_nat
                   triplet_counter = triplet_counter + 1

                   !Triplet
                   ph%Index_i(triplet_counter) = na1
                   ph%Index_j(triplet_counter) = na2
                   ph%Index_k(triplet_counter) = na3

                   !Positions of the 2nd and 3rd atom in the triplet
                   !converted to Cartesian coordinates (Ang).
                   ph%R_j(:, triplet_counter) = matmul(crys%lattvecs, R2(:, ii)/dble(ph%qmesh))*10.0_dp
                   ph%R_k(:, triplet_counter) = matmul(crys%lattvecs, R3(:, ii)/dble(ph%qmesh))*10.0_dp
                   
                   do j1 =1, 3
                      jn1 = j1 + (na1 - 1)*3
                      do j2 =1, 3
                         jn2 = j2 + (na2 - 1)*3
                         do j3 =1, 3
                            jn3 = j3 + (na3 - 1)*3
                            ph%ifc3(j1, j2, j3, triplet_counter) = fc(jn1, jn2, jn3, ii)
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
       ph%ifc3 = ph%ifc3*Ryd2eV/(bohr2nm*10.0_dp)**3 !eV/Ang^3
    end if

    !Each vector is rounded to the nearest lattice vector.
    tmp = crys%lattvecs
    call dgesv(3, ph%numtriplets, tmp, 3, P, ph%R_j, 3, info)
    ph%R_j = matmul(crys%lattvecs, anint(ph%R_j/10.0_dp)) !nm
    tmp = crys%lattvecs
    call dgesv(3, ph%numtriplets, tmp, 3, P, ph%R_k, 3, info)
    ph%R_k = matmul(crys%lattvecs, anint(ph%R_k/10.0_dp)) !nm
  end subroutine read_ifc3

  subroutine phonon_espresso(ph, crys, nk, kpoints, omegas, eigenvect, velocities)
    !! Subroutine to calculate phonons from the 2nd order force constants.
    !! This is adapted from Quantum Espresso and ShengBTE.

    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    integer(k8), intent(in) :: nk
    real(dp), intent(in) :: kpoints(nk, 3)
    real(dp), intent(out) :: omegas(nk, ph%numbranches)
    real(dp), optional, intent(out) :: velocities(nk, ph%numbranches, 3)
    complex(kind=8), optional, intent(out) :: eigenvect(nk, ph%numbranches, ph%numbranches)

    ! QE's 2nd-order files are in Ryd units.
    real(kind=8),parameter :: toTHz=20670.687,&
         massfactor=1.8218779*6.022e-4

    integer(k8) :: ir,nreq,ntype,nat,nbranches
    integer(k8) :: i,j,ipol,jpol,iat,jat,idim,jdim,t1,t2,t3,m1,m2,m3,ik
    integer(k8) :: ndim,nwork,ncell_g(3)
    integer(k8),allocatable :: tipo(:)
    character(len=5),allocatable :: label(:)
    real(dp) :: weight,total_weight,exp_g,ck
    real(dp) :: r_ws(3), at(3,3)
    real(dp) :: alpha,geg,gmax,kt,gr,volume_r,dnrm2
    real(dp) :: zig(3),zjg(3),dgeg(3),t(0:3),g(0:3),g_old(0:3)
    real(dp), allocatable :: omega2(:),rwork(:)
    real(dp),allocatable :: k(:,:),mass(:),r(:,:),eps(:,:)
    real(dp),allocatable :: eival(:,:),vels(:,:,:),zeff(:,:,:)
    complex(dp) :: auxi(3)
    complex(dp),allocatable :: cauxiliar(:),eigenvectors(:,:),work(:)
    complex(dp),allocatable :: dyn(:,:),dyn_s(:,:,:),dyn_g(:,:,:)
    complex(dp),allocatable :: ddyn(:,:,:),ddyn_s(:,:,:,:),ddyn_g(:,:,:,:)

    ! Quantum Espresso's 2nd-order format contains information about
    ! lattice vectors, atomic positions, Born effective charges and so
    ! forth in its header. The information is read but completely
    ! ignored. It is the user's responsibility to ensure that
    ! it is consistent with the CONTROL file.
    nwork = 1
    ntype = crys%numelements
    nat = crys%numatoms
    ndim = 3*nat
    nbranches = ndim
    
    allocate(omega2(nbranches))
    allocate(work(nwork))
    allocate(rwork(max(1,9*nat-2)))
    allocate(k(nk,3))
    allocate(mass(ntype))
    allocate(tipo(nat))
    allocate(eps(3,3))
    allocate(zeff(nat,3,3))    
    allocate(dyn(ndim,ndim))
    allocate(dyn_s(nk,ndim,ndim))
    allocate(dyn_g(nk,ndim,ndim))
    allocate(eival(ndim,nk))
    allocate(eigenvectors(ndim,ndim))
    allocate(cauxiliar(ndim))
    if(present(velocities)) then
       allocate(ddyn(ndim,ndim,3))
       allocate(ddyn_s(nk,ndim,ndim,3))
       allocate(ddyn_g(nk,ndim,ndim,3))
       allocate(vels(ndim,nk,3))
    end if
    
    mass = crys%masses/massfactor
    tipo = crys%atomtypes
    eps = transpose(crys%epsilon)
    do i = 1, nat
       zeff(i, :, :) = transpose(crys%born(:, :, i))
    end do
    
    ! Make sure operations are performed in consistent units.
    do ik = 1, nk
       k(ik, :) = matmul(crys%reclattvecs, kpoints(ik, :))
    end do
    k = k*bohr2nm
    
    volume_r = crys%volume/bohr2nm**3

    gmax=14.
    alpha=(twopi*bohr2nm/dnrm2(3,crys%lattvecs(:,1),1))**2
    geg=gmax*4.*alpha
    ncell_g=int(sqrt(geg)/ph%cell_g(:,0))+1
    
    dyn_s = 0.0_dp
    if(present(velocities)) ddyn_s = 0.0_dp
    
    do iat=1,nat
       do jat=1,nat
          total_weight=0.0d0
          do m1=-2*ph%scell(1),2*ph%scell(1)
             do m2=-2*ph%scell(2),2*ph%scell(2)
                do m3=-2*ph%scell(3),2*ph%scell(3)
                   do i=1,3
                      t(i)=m1*ph%cell_r(1,i)+m2*ph%cell_r(2,i)+m3*ph%cell_r(3,i)
                      r_ws(i)=t(i)+ph%rr(iat,jat,i)
                   end do
                   weight=0.d0
                   nreq=1
                   j=0
                   Do ir=1,124
                      ck=dot_product(r_ws,ph%rws(ir,1:3))-ph%rws(ir,0)
                      if(ck.gt.1e-6) then
                         j=1
                         cycle
                      end if
                      if(abs(ck).lt.1e-6) then
                         nreq=nreq+1
                      end if
                   end do
                   if(j.eq.0) then
                      weight=1.d0/dble(nreq)
                   end if
                   if(weight.gt.0.d0) then
                      t1=mod(m1+1,ph%scell(1))
                      if(t1.le.0) then
                         t1=t1+ph%scell(1)
                      end if
                      t2=mod(m2+1,ph%scell(2))
                      if(t2.Le.0) then
                         t2=t2+ph%scell(2)
                      end if
                      t3=mod(m3+1,ph%scell(3))
                      if(t3.le.0) then
                         t3=t3+ph%scell(3)
                      end if
                      do ik=1,nk
                         kt=dot_product(k(ik,1:3),t(1:3))
                         do ipol=1,3
                            idim = (iat-1)*3+ipol
                            do jpol=1,3
                               jdim = (jat-1)*3+jpol
                               dyn_s(ik,idim,jdim)=dyn_s(ik,idim,jdim)+&
                                    ph%ifc2(ipol,jpol,iat,jat,t1,t2,t3)*&
                                    expi(-kt)*weight
                               if(present(velocities)) then
                                  ddyn_s(ik,idim,jdim,1:3)=ddyn_s(ik,idim,jdim,1:3)-&
                                       oneI*t(1:3)*&
                                       ph%ifc2(ipol,jpol,iat,jat,t1,t2,t3)*&
                                       expi(-kt)*weight
                               end if
                            end do
                         end do
                      end do
                   end if
                   total_weight=total_weight+weight
                end do
             end do
          end do
       end do
    end do
    ! The nonanalytic correction has two components in this
    ! approximation. Results may differ slightly between this method
    ! and the one implemented in the previous subroutine.
    dyn_g = 0.0_dp
    if(present(velocities)) ddyn_g = 0.0_dp
    if(crys%polar) then
       do m1=-ncell_g(1),ncell_g(1)
          do m2=-ncell_g(2),ncell_g(2)
             do m3=-ncell_g(3),ncell_g(3)
                g(1:3)=m1*ph%cell_g(1,1:3)+&
                     m2*ph%cell_g(2,1:3)+m3*ph%cell_g(3,1:3)
                geg=dot_product(g(1:3),matmul(eps,g(1:3)))
                if(geg.gt.0.0d0.and.geg/alpha/4.0d0.lt.gmax) then
                   exp_g=exp(-geg/alpha/4.0d0)/geg
                   do iat=1,nat
                      zig(1:3)=matmul(g(1:3),zeff(iat,1:3,1:3))
                      auxi(1:3)=0.
                      do jat=1,nat
                         gr=dot_product(g(1:3),ph%rr(iat,jat,1:3))
                         zjg(1:3)=matmul(g(1:3),zeff(jat,1:3,1:3))
                         auxi(1:3)=auxi(1:3)+zjg(1:3)*expi(gr)
                      end do
                      do ipol=1,3
                         idim=(iat-1)*3+ipol
                         do jpol=1,3
                            jdim=(iat-1)*3+jpol
                            dyn_g(1:nk,idim,jdim)=dyn_g(1:nk,idim,jdim)-&
                                 exp_g*zig(ipol)*auxi(jpol)
                         end do
                      end do
                   end do
                end if
                g_old(0:3)=g(0:3)
                do ik=1,nk
                   g(1:3)=g_old(1:3)+k(ik,1:3)
                   geg=dot_product(g(1:3),matmul(eps,g(1:3)))
                   if (geg.gt.0.0d0.and.geg/alpha/4.0d0.lt.gmax) then
                      exp_g=exp(-geg/alpha/4.0d0)/geg
                      dgeg=matmul(eps+transpose(eps),g(1:3))
                      do iat=1,nat
                         zig(1:3)=matmul(g(1:3),zeff(iat,1:3,1:3))
                         do jat=1,nat
                            gr=dot_product(g(1:3),ph%rr(iat,jat,1:3))
                            zjg(1:3)=matmul(g(1:3),zeff(jat,1:3,1:3))
                            do ipol=1,3
                               idim=(iat-1)*3+ipol
                               do jpol=1,3
                                  jdim=(jat-1)*3+jpol
                                  dyn_g(ik,idim,jdim)=dyn_g(ik,idim,jdim)+&
                                       exp_g*zig(ipol)*zjg(jpol)*expi(gr)
                                  if(present(velocities)) then
                                     do i=1,3
                                        ddyn_g(ik,idim,jdim,i)=ddyn_g(ik,idim,jdim,i)+&
                                             exp_g*expi(gr)*&
                                             (zjg(jpol)*zeff(iat,i,ipol)+zig(ipol)*zeff(jat,i,jpol)+&
                                             zig(ipol)*zjg(jpol)*oneI*ph%rr(iat,jat,i)-&
                                             zig(ipol)*zjg(jpol)*(dgeg(i)/alpha/4.0+dgeg(i)/geg))
                                     end do
                                  end if
                               end do
                            end do
                         end do
                      end do
                   end if
                end do
             end do
          end do
       end do
       dyn_g = dyn_g*8.0_dp*pi/volume_r
       if(present(velocities)) ddyn_g = ddyn_g*8.0_dp*pi/volume_r
    end if
    
    ! Once the dynamical matrix has been built, the frequencies and
    ! group velocities are extracted exactly like in the previous
    ! subroutine.
    do ik = 1, nk
       dyn(:,:) = dyn_s(ik,:,:) + dyn_g(ik,:,:)
       if(present(velocities)) then
          ddyn(:,:,:) = ddyn_s(ik,:,:,:) + ddyn_g(ik,:,:,:)
       end if

       do ipol=1,3
          do jpol=1,3
             do iat=1,nat
                do jat=1,nat
                   idim=(iat-1)*3+ipol
                   jdim=(jat-1)*3+jpol
                   dyn(idim,jdim)=dyn(idim,jdim)/ph%mm(iat,jat)
                   if(present(velocities)) then
                      ddyn(idim,jdim,1:3)=ddyn(idim,jdim,1:3)/ph%mm(iat,jat)
                   end if
                end do
             end do
          end do
       end do
       
       call zheev("V","U",nbranches,dyn(:,:),nbranches,omega2,work,-1_k8,rwork,i)
       if(real(work(1)).gt.nwork) then
          nwork=nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V","U",nbranches,dyn(:,:),nbranches,omega2,work,nwork,rwork,i)
       
       if(present(eigenvect)) then
          eigenvect(ik,:,:) = transpose(dyn(:,:))
       end if
       
       omegas(ik,:)=sign(sqrt(abs(omega2)),omega2)

       if(present(velocities)) then
          do i=1,nbranches
             do j=1,3
                velocities(ik,i,j)=real(dot_product(dyn(:,i),&
                     matmul(ddyn(:,:,j),dyn(:,i))))
             end do
             velocities(ik,i,:)=velocities(ik,i,:)/(2.*omegas(ik,i))
          end do
       end if

       !Take care of gamma point.
       if(all(k(ik,1:3) == 0)) then
          omegas(ik, 1:3) = 0.0_dp
          if(present(velocities)) velocities(ik, :, :) = 0.0_dp
       end if
    end do

    !Units conversion
    omegas=omegas*Ryd2eV !eV
    if(present(velocities)) velocities=velocities*toTHz*bohr2nm !Km/s
  end subroutine phonon_espresso
end module phonon_module
