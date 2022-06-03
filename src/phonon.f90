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
!
! 20220322 : Added support for reading d3q's sparse files (XC)

module phonon_module
  !! Module containing type and procedures related to the phononic properties.

  use params, only: dp, k8, bohr2nm, pi, twopi, Ryd2eV, oneI
  use particle_module, only: particle
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
  public phonon
  
  type, extends(particle) :: phonon
     !! Data and procedures related to phonons.
     
     character(len = 2) :: prefix = 'ph'
     !! Prefix idenitfying particle type.
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

     !Data read from ifc2 file. These will be used in the phonon calculation.
     real(dp) :: rws(124, 0:3), cell_r(1:3, 0:3), cell_g(1:3, 0:3)
     real(dp), allocatable :: mm(:,:), rr(:,:,:)
      
   contains

     procedure, public :: initialize, deallocate_phonon_quantities
     procedure, private :: calculate_phonons, read_ifc2, read_ifc3, &
          phonon_espresso
     
  end type phonon

contains

  subroutine initialize(self, crys, sym, num)
    !! Initialize the phonon data type, calculate ground state phonon properties,
    !! and read 3rd order force constants data. 

    class(phonon), intent(out) :: self
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(numerics), intent(in) :: num

    call subtitle("Setting up phonons...")

    !Set phonon branches
    self%numbands = crys%numatoms*3
    !Set wave vector mesh
    self%wvmesh = num%qmesh
    !Set number of phonon wave vectors
    self%nwv = product(self%wvmesh(:))

    !Read ifc2 and related quantities
    call read_ifc2(self, crys)
    
    !Calculate harmonic properties
    call calculate_phonons(self, crys, sym, num)

    if(.not. num%onlyebte) then
       !Read ifc3s and related quantities
       call read_ifc3(self, crys)
    end if
  end subroutine initialize

  subroutine deallocate_phonon_quantities(self)
    !! Deallocate the electron eigenvectors

    class(phonon), intent(inout) :: self

    deallocate(self%evecs, self%ifc2, self%ifc3, self%Index_i, self%Index_j, self%Index_k, &
         self%mm, self%rr)
  end subroutine deallocate_phonon_quantities
  
  subroutine calculate_phonons(self, crys, sym, num)
    !! Calculate phonon quantities on the FBZ and IBZ meshes.

    class(phonon), intent(inout) :: self
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
    call distribute_points(self%nwv, chunk, start, end, num_active_images)

    !Allocate small work variable chunk for each image
    allocate(ens_chunk(chunk, self%numbands)[*])
    allocate(vels_chunk(chunk, self%numbands, 3)[*])
    allocate(evecs_chunk(chunk, self%numbands, self%numbands)[*])

    !Calculate FBZ mesh
    call calculate_wavevectors_full(self%wvmesh, self%wavevecs, blocks)

    !Print phonon FBZ mesh
    call write2file_rank2_real("ph.wavevecs_fbz", self%wavevecs)
    
    !Calculate FBZ phonon quantities
    call phonon_espresso(self, crys, chunk, self%wavevecs(start:end, :), &
         ens_chunk, evecs_chunk, vels_chunk)

    !Gather the chunks from the images and broadcast to all
    allocate(self%ens(self%nwv, self%numbands))
    allocate(self%vels(self%nwv, self%numbands, 3))
    allocate(self%evecs(self%nwv, self%numbands, self%numbands))
    sync all
    if(this_image() == 1) then
       do im = 1, num_active_images
          self%ens(start[im]:end[im], :) = ens_chunk(:,:)[im]
          self%vels(start[im]:end[im], :, :) = vels_chunk(:,:,:)[im]
          self%evecs(start[im]:end[im], :, :) = evecs_chunk(:,:,:)[im]
       end do
    end if
    call co_broadcast(self%ens, 1)
    call co_broadcast(self%vels, 1)
    call co_broadcast(self%evecs, 1)
    
    deallocate(ens_chunk, vels_chunk, evecs_chunk)
    
    !Calculate IBZ mesh
    call print_message("Calculating IBZ and IBZ -> FBZ mappings...")
    call find_irred_wedge(self%wvmesh, self%nwv_irred, self%wavevecs_irred, &
         self%indexlist_irred, self%nequiv, sym%nsymm_rot, sym%qrotations, &
         self%ibz2fbz_map, self%equiv_map, blocks)

    !Print phonon IBZ mesh
    call write2file_rank2_real("ph.wavevecs_ibz", self%wavevecs_irred)
    
    !Create symmetrizers of wave vector dependent vectors ShengBTE style
    allocate(symmetrizers_chunk(3, 3, chunk)[*])
    symmetrizers_chunk = 0.0_dp
    do iq = start, end
       kk = 0
       do jj = 1, sym%nsymm
          if(self%equiv_map(jj, iq) == iq) then
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

    !Gather from images and broadcast to all
    allocate(self%symmetrizers(3, 3, self%nwv))
    sync all
    if(this_image() == 1) then
       do im = 1, num_active_images
          self%symmetrizers(:, :, start[im]:end[im]) = symmetrizers_chunk(:,:,:)[im]
       end do
    end if
    call co_broadcast(self%symmetrizers, 1)
    
    deallocate(symmetrizers_chunk)
    
    !Symmetrize phonon energies and velocities.
    do i = 1, self%nwv_irred !an irreducible point
       ii = self%indexlist_irred(i)
       self%vels(ii,:,:)=transpose(&
            matmul(self%symmetrizers(:,:,ii),transpose(self%vels(ii,:,:))))
       do l = 1, self%nequiv(i) !number of equivalent points of i
          il = self%ibz2fbz_map(l, i, 2) ! (i, l) -> il
          s = self%ibz2fbz_map(l, i, 1) ! mapping rotation

          !energy
          self%ens(il,:) = self%ens(ii,:)

          !velocity
          do ib = 1, self%numbands
             !here use real space (Cartesian) rotations
             self%vels(il, ib, :) = matmul(sym%crotations(:, :, s), self%vels(ii, ib, :))
          end do
       end do
    end do
    
    !Print out irreducible phonon energies and velocities
    if(this_image() == 1) then
       write(numcols, "(I0)") self%numbands
       open(1, file = "ph.ens_ibz", status = "replace")
       do iq = 1, self%nwv_irred
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               self%ens(self%indexlist_irred(iq), :)
       end do
       close(1)

       write(numcols, "(I0)") 3*self%numbands
       open(1, file = "ph.vels_ibz", status = "replace")
       do iq = 1, self%nwv_irred
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               self%vels(self%indexlist_irred(iq), :, :)
       end do
       close(1)
    end if
    
    !Calculate phonon tetrahedra
    if(num%tetrahedra) then
       call print_message("Calculating phonon mesh tetrahedra...")
       call form_tetrahedra_3d(self%nwv, self%wvmesh, self%tetra, self%tetracount, &
            self%tetramap, .false.)
       call fill_tetrahedra_3d(self%tetra, self%ens, self%tetra_evals)
    else
       call print_message("Calculating phonon mesh triangles...")
       call form_triangles(self%nwv, self%wvmesh, self%triang, self%triangcount, &
            self%triangmap, .false.)
       call fill_triangles(self%triang, self%ens, self%triang_evals)
    end if
  end subroutine calculate_phonons
  
  subroutine read_ifc2(self, crys)
    !! Read the 2nd order force constants from the Quantum Espresso format.
    !! This is adapted from ShengBTE.
    
    class(phonon), intent(inout) :: self
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
    
    allocate(self%mm(crys%numatoms, crys%numatoms))
    allocate(self%rr(crys%numatoms, crys%numatoms, 3))
    
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

    self%scell = qscell
    
    !Read the force constants.
    allocate(self%ifc2(3, 3, nat, nat, self%scell(1), self%scell(2), self%scell(3)))
    nfc2 = 3*3*nat*nat
    do i = 1, nfc2
       read(1, *) ipol, jpol, iat, jat
       do j = 1, self%scell(1)*self%scell(2)*self%scell(3)
          read(1, *) t1, t2, t3, &
               self%ifc2(ipol, jpol, iat, jat, t1, t2, t3)
       end do
    end do
    close(1)
    
    !Enforce the conservation of momentum in the simplest way possible.
    do i = 1, 3
       do j = 1, 3
          do iat = 1, nat
             self%ifc2(i, j, iat, iat, 1, 1, 1) = self%ifc2(i, j, iat, iat, 1, 1, 1) - &
                  sum(self%ifc2(i, j, iat, :, :, :, :))
          end do
       end do
    end do

    self%cell_r(:, 1:3) = transpose(crys%lattvecs)/bohr2nm
    do i = 1, 3
       self%cell_r(i, 0) = dnrm2(3, self%cell_r(i, 1:3), 1)
    end do
    self%cell_g(:, 1:3) = transpose(crys%reclattvecs)*bohr2nm
    do i = 1, 3
       self%cell_g(i, 0) = dnrm2(3, self%cell_g(i, 1:3), 1)
    end do

    wscell(1,1:3) = self%cell_r(1,1:3)*self%scell(1)
    wscell(2,1:3) = self%cell_r(2,1:3)*self%scell(2)
    wscell(3,1:3) = self%cell_r(3,1:3)*self%scell(3)

    j = 1
    do m1 = -2, 2
       do m2 = -2, 2
          do m3 = -2, 2
             if(all((/m1, m2, m3/).eq.0)) then
                cycle
             end if
             do i = 1, 3
                self%rws(j, i) = wscell(1, i)*m1 + wscell(2, i)*m2 + wscell(3, i)*m3
             end do
             self%rws(j, 0) = 0.5*dot_product(self%rws(j, 1:3), self%rws(j, 1:3))
             j = j + 1
          end do
       end do
    end do

    do i = 1, nat
       self%mm(i, i) = mass(tipo(i))
       self%rr(i, i, :) = 0
       do j = i + 1, nat
          self%mm(i, j) = sqrt(mass(tipo(i))*mass(tipo(j)))
          self%rr(i, j, 1:3) = r(i, 1:3) - r(j, 1:3)
          self%mm(j, i) = self%mm(i, j)
          self%rr(j, i, 1:3) = -self%rr(i, j, 1:3)
       end do
    end do
    
  end subroutine read_ifc2
  
  subroutine read_ifc3(self, crys)
    !! Read the 3rd order force constants in the thirdorder.py format.
    !! This subroutine is adapted from ShengBTE.

    class(phonon), intent(inout) :: self
    type(crystal), intent(in) :: crys
    
    !Local variables
    real(dp) :: tmp(3,3), r(crys%numatoms, 3), celldm(6), at(3,3), &
         mass(crys%numelements), zeff(crys%numatoms, 3, 3), eps(3, 3), fc_
    real(dp), allocatable :: fc(:,:,:,:)
    integer(k8) :: ii, jj, ll, mm, nn, ltem, mtem, ntem, info, P(3), &
         na1, na2, na3, j1, j2, j3, na1_, na2_, na3_, j1_, j2_, j3_, &
         triplet_counter, nR, nR_, qscell(3), tipo(crys%numatoms), i, j, &
         ibrav, ntype, nat, jn1, jn2, jn3, ind
    integer(k8), allocatable :: nind(:), R2tmp(:,:), R3tmp(:,:), &
         R2(:,:), R3(:,:), triplet_map(:,:,:,:)
    character(len = 1) :: polar_key
    character(len = 6) :: label(crys%numelements), sparse_header
    logical :: sheng_file_exists, d3q_file_exists, &
         d3q_sparse_file_exists, save_nR

    !External procedures
    external :: dgesv

    !Check what force constants files have been provided
    sheng_file_exists = .False.
    d3q_file_exists = .False.
    d3q_sparse_file_exists = .False.
    inquire(file = 'FORCE_CONSTANTS_3RD', exist = sheng_file_exists)
    inquire(file = 'mat3R', exist = d3q_file_exists)
    inquire(file = 'mat3R.sparse', exist = d3q_sparse_file_exists)

    if(sheng_file_exists .and. (.not. d3q_file_exists) .and. (.not. d3q_sparse_file_exists)) then
       call print_message('Reading ShengBTE format third order force constants...')
    else if(d3q_file_exists .and. (.not. sheng_file_exists) .and. (.not. d3q_sparse_file_exists)) then
       call print_message('Reading d3q format third order force constants...')
    else if(d3q_sparse_file_exists .and. (.not. sheng_file_exists) .and. (.not. d3q_file_exists)) then
       call print_message('Reading d3q_sparse format third order force constants...')
    else if(sheng_file_exists .and. d3q_file_exists .and. (.not. d3q_sparse_file_exists)) then
       call print_message(&
            'Both ShengBTE and d3q format third order force constants provided. Defaulting to ShengBTE format.')
       d3q_file_exists = .False.
    else if(sheng_file_exists .and. (.not. d3q_file_exists) .and. d3q_sparse_file_exists) then
       call print_message(&
            'Both ShengBTE and d3q_sparse format third order force constants provided. Defaulting to ShengBTE format.')
       d3q_sparse_file_exists = .False.
    else if((.not. sheng_file_exists) .and. d3q_file_exists .and. d3q_sparse_file_exists) then
       call print_message(&
            'Both d3q and d3q_sparse format third order force constants provided. Defaulting to d3q_sparse format.')
       d3q_file_exists = .False.
    else if(sheng_file_exists .and. d3q_file_exists .and. d3q_sparse_file_exists) then
       call print_message(&
            'All ShengBTE, d3q and d3q_sparse format third order force constants provided. Defaulting to ShengBTE format.')
       d3q_file_exists = .False.
       d3q_sparse_file_exists = .False.
    else
       call exit_with_message('Third order force constant file not provided. Exiting.')
    end if

    if(sheng_file_exists) then
       !The file is in a simple sparse format, described in detail in
       !the user documentation. See Doc/ShengBTE.pdf.
       open(1, file = 'FORCE_CONSTANTS_3RD', status = "old")
       read(1, *) self%numtriplets
       allocate(self%Index_i(self%numtriplets), self%Index_j(self%numtriplets), &
            self%Index_k(self%numtriplets))
       allocate(self%ifc3(3, 3, 3, self%numtriplets), self%R_j(3, self%numtriplets), &
            self%R_k(3,self%numtriplets))
       do ii = 1, self%numtriplets
          read(1, *) jj
          read(1, *) self%R_j(:, ii) !Ang
          read(1, *) self%R_k(:, ii) !Ang
          read(1, *) self%Index_i(ii), self%Index_j(ii), self%Index_k(ii)
          do ll = 1, 3
             do mm = 1, 3
                do nn = 1, 3
                   read(1, *) ltem, mtem, ntem, self%ifc3(ll, mm, nn, ii)
                end do
             end do
          end do
       end do
       close(1)
       !IFC3 units are eV/Ang^3
    else if(d3q_file_exists) then
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
       
       save_nR = .true.
       do na1 = 1, crys%numatoms
          do na2 = 1, crys%numatoms
             do na3 = 1, crys%numatoms
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
                                 "Wrong triplet indices and/or tensor element location in mat3R file. Exiting.")
                         end if

                         !Read number of unit cells in file.
                         read(1, *) nR_
                         !Save this number only the first time it is read.
                         !Also, allocate the various quantities.
                         if(save_nR) then
                            nR = nR_
                            allocate(R2(3, nR), R3(3, nR), &
                                 fc(self%numbands, self%numbands, self%numbands, nR))
                            save_nR = .false.
                         end if

                         if(nR_ /= nR) call exit_with_message(&
                              "Wrong number of unit cells in mat3R file. Exiting.")

                         do ii = 1, nR
                            !R2 and R3 are the same for every nR chunk.
                            read(1, *) R2(:, ii), R3(:, ii), fc(jn1, jn2, jn3, ii)
                            !At this point the fc units are Ry/Bohr^3
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
       
       !Number of triplets
       self%numtriplets = nR*crys%numatoms**3
       
       !Allocate quantities
       allocate(self%Index_i(self%numtriplets), self%Index_j(self%numtriplets), self%Index_k(self%numtriplets))
       allocate(self%ifc3(3, 3, 3, self%numtriplets), self%R_j(3, self%numtriplets), self%R_k(3,self%numtriplets))
       
       !Convert to the standard format.
       triplet_counter = 0
       do ii = 1, nR
          do na1 = 1, crys%numatoms
             do na2 = 1, crys%numatoms
                do na3 = 1, crys%numatoms
                   triplet_counter = triplet_counter + 1

                   !Triplet
                   self%Index_i(triplet_counter) = na1
                   self%Index_j(triplet_counter) = na2
                   self%Index_k(triplet_counter) = na3

                   !Positions of the 2nd and 3rd atom in the triplet
                   !converted to Cartesian coordinates (Ang).
                   self%R_j(:, triplet_counter) = &
                        matmul(crys%lattvecs, R2(:, ii)*10.0_dp) !Ang
                   self%R_k(:, triplet_counter) = &
                        matmul(crys%lattvecs, R3(:, ii)*10.0_dp) !Ang
                   
                   do j1 =1, 3
                      jn1 = j1 + (na1 - 1)*3
                      do j2 =1, 3
                         jn2 = j2 + (na2 - 1)*3
                         do j3 =1, 3
                            jn3 = j3 + (na3 - 1)*3
                            self%ifc3(j1, j2, j3, triplet_counter) = fc(jn1, jn2, jn3, ii)
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
       self%ifc3 = self%ifc3*Ryd2eV/(bohr2nm*10.0_dp)**3 !eV/Ang^3
    else if(d3q_sparse_file_exists) then
       !See SUBROUTINE read_fc3_sparse of fc3_interp.f90 of the d3q code
       !for more information about the format.
       open(1, file = 'mat3R.sparse', status = "old")

       !Check file really in sparse mode
       read(1,*) sparse_header
       if (sparse_header /= "sparse") then
         call exit_with_message('Not really a d3q sparse file. Exiting.')
       end if
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

       !Read number of unit cells in file.
       read(1, *) nR
       !Also, allocate the various quantities.
       allocate(nind(nR), R2tmp(3,nR), R3tmp(3,nR))
       allocate(triplet_map(crys%numatoms,crys%numatoms,crys%numatoms,nR))

       do ii = 1, nR
          ! Read how many (atom indices + components) at that position, and atom coordinates
          read(1, *) nind(ii), R2tmp(:, ii), R3tmp(:, ii)
       end do

       ! Will need to make two passes on the rest of file
       !  1st pass will give number of triplets
       !  2nd pass will do the reading proper

       !Number of triplets
       triplet_map = 0
       self%numtriplets = 0
       do ii = 1, nR
          do ind = 1, nind(ii)
             read(1,*) jn1, jn2, jn3, fc_

             ! Unwrap the indices
             na1 = (jn1-1)/3 + 1
             na2 = (jn2-1)/3 + 1
             na3 = (jn3-1)/3 + 1

             ! Compute atom_indices_visited index
             if ( triplet_map(na1,na2,na3,ii)==0 ) then
                self%numtriplets = self%numtriplets + 1
                triplet_map(na1,na2,na3,ii) = self%numtriplets
             end if

          end do
          read(1, *)
       end do
       
       !Allocate quantities
       allocate(self%Index_i(self%numtriplets), self%Index_j(self%numtriplets), self%Index_k(self%numtriplets))
       allocate(self%ifc3(3, 3, 3, self%numtriplets), self%R_j(3, self%numtriplets), self%R_k(3,self%numtriplets))
       allocate(R2(3, self%numtriplets), R3(3, self%numtriplets))
       
       ! Read the FC3s proper
       rewind(1)
       read(1,*) sparse_header
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

       !Read number of unit cells in file.
       read(1, *) nR

       do ii = 1, nR
          ! Read how many (atom indices + components) at that position, and atom coordinates
          read(1, *) nind(ii), tmp(:, 2), tmp(:, 3)
       end do

       self%ifc3 = 0._dp
       do ii = 1, nR
          do ind = 1, nind(ii)
             read(1,*) jn1, jn2, jn3, fc_

             ! Unwrap the indices
             na1 = (jn1-1)/3 + 1
             na2 = (jn2-1)/3 + 1
             na3 = (jn3-1)/3 + 1
             j1  = mod(jn1-1, 3) + 1
             j2  = mod(jn2-1, 3) + 1
             j3  = mod(jn3-1, 3) + 1

             !Triplet
             triplet_counter = triplet_map(na1,na2,na3,ii)
             self%Index_i(triplet_counter) = na1
             self%Index_j(triplet_counter) = na2
             self%Index_k(triplet_counter) = na3
             self%ifc3(j1, j2, j3, triplet_counter) = fc_
             R2(:,triplet_counter) = R2tmp(:,ii)
             R3(:,triplet_counter) = R3tmp(:,ii)

          end do
          read(1, *)
       end do

       ! Finish converting to the standard format.
       do ii = 1, self%numtriplets
          !Positions of the 2nd and 3rd atom in the triplet
          !converted to Cartesian coordinates (Ang).
          self%R_j(:, ii) = &
               matmul(crys%lattvecs, R2(:, ii)*10.0_dp) !Ang
          self%R_k(:, ii) = &
               matmul(crys%lattvecs, R3(:, ii)*10.0_dp) !Ang
       end do
       self%ifc3 = self%ifc3*Ryd2eV/(bohr2nm*10.0_dp)**3 !eV/Ang^3
    end if

    !Each vector is rounded to the nearest lattice vector.
    tmp = crys%lattvecs
    call dgesv(3, self%numtriplets, tmp, 3, P, self%R_j, 3, info)
    self%R_j = matmul(crys%lattvecs, anint(self%R_j/10.0_dp)) !nm
    tmp = crys%lattvecs
    call dgesv(3, self%numtriplets, tmp, 3, P, self%R_k, 3, info)
    self%R_k = matmul(crys%lattvecs, anint(self%R_k/10.0_dp)) !nm

    if(this_image() == 1) &
       write(*, "(A, I10)") " Number triplets read in = ", self%numtriplets
  end subroutine read_ifc3

  subroutine phonon_espresso(self, crys, nk, kpoints, omegas, eigenvect, velocities)
    !! Subroutine to calculate phonons from the 2nd order force constants.
    !
    ! This is adapted from ShengBTE's subroutine of the same name.
    ! ShengBTE is distributed under GPL v3 or later.
    
    class(phonon), intent(in) :: self
    type(crystal), intent(in) :: crys
    integer(k8), intent(in) :: nk
    real(dp), intent(in) :: kpoints(nk, 3)
    real(dp), intent(out) :: omegas(nk, self%numbands)
    real(dp), optional, intent(out) :: velocities(nk, self%numbands, 3)
    complex(dp), optional, intent(out) :: eigenvect(nk, self%numbands, self%numbands)

    ! QE's 2nd-order files are in Ryd units.
    real(dp),parameter :: toTHz=20670.687,&
         massfactor=1.8218779*6.022e-4

    integer(k8) :: ir,nreq,ntype,nat,nbranches
    integer(k8) :: i,j,ipol,jpol,iat,jat,idim,jdim,t1,t2,t3,m1,m2,m3,ik
    integer(k8) :: ndim,nwork,ncell_g(3)
    integer(k8),allocatable :: tipo(:)
    real(dp) :: weight,total_weight,exp_g,ck
    real(dp) :: r_ws(3)
    real(dp) :: alpha,geg,gmax,kt,gr,volume_r,dnrm2
    real(dp) :: zig(3),zjg(3),dgeg(3),t(0:3),g(0:3),g_old(0:3)
    real(dp), allocatable :: omega2(:),rwork(:)
    real(dp),allocatable :: k(:,:),mass(:),eps(:,:)
    real(dp),allocatable :: eival(:,:),vels(:,:,:),zeff(:,:,:)
    complex(dp) :: auxi(3)
    complex(dp),allocatable :: cauxiliar(:),eigenvectors(:,:),work(:)
    complex(dp),allocatable :: dyn(:,:),dyn_s(:,:,:),dyn_g(:,:,:)
    complex(dp),allocatable :: ddyn(:,:,:),ddyn_s(:,:,:,:),ddyn_g(:,:,:,:)

    !External procedures
    external :: zheev
    
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
    ncell_g=int(sqrt(geg)/self%cell_g(:,0))+1
    
    dyn_s = 0.0_dp
    if(present(velocities)) ddyn_s = 0.0_dp
    
    do iat=1,nat
       do jat=1,nat
          total_weight=0.0d0
          do m1=-2*self%scell(1),2*self%scell(1)
             do m2=-2*self%scell(2),2*self%scell(2)
                do m3=-2*self%scell(3),2*self%scell(3)
                   do i=1,3
                      t(i)=m1*self%cell_r(1,i)+m2*self%cell_r(2,i)+m3*self%cell_r(3,i)
                      r_ws(i)=t(i)+self%rr(iat,jat,i)
                   end do
                   weight=0.d0
                   nreq=1
                   j=0
                   Do ir=1,124
                      ck=dot_product(r_ws,self%rws(ir,1:3))-self%rws(ir,0)
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
                      t1=mod(m1+1,self%scell(1))
                      if(t1.le.0) then
                         t1=t1+self%scell(1)
                      end if
                      t2=mod(m2+1,self%scell(2))
                      if(t2.Le.0) then
                         t2=t2+self%scell(2)
                      end if
                      t3=mod(m3+1,self%scell(3))
                      if(t3.le.0) then
                         t3=t3+self%scell(3)
                      end if
                      do ik=1,nk
                         kt=dot_product(k(ik,1:3),t(1:3))
                         do ipol=1,3
                            idim = (iat-1)*3+ipol
                            do jpol=1,3
                               jdim = (jat-1)*3+jpol
                               dyn_s(ik,idim,jdim)=dyn_s(ik,idim,jdim)+&
                                    self%ifc2(ipol,jpol,iat,jat,t1,t2,t3)*&
                                    expi(-kt)*weight
                               if(present(velocities)) then
                                  ddyn_s(ik,idim,jdim,1:3)=ddyn_s(ik,idim,jdim,1:3)-&
                                       oneI*t(1:3)*&
                                       self%ifc2(ipol,jpol,iat,jat,t1,t2,t3)*&
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

    !The nonanalytic correction has two components in this approximation.
    dyn_g = 0.0_dp
    if(present(velocities)) ddyn_g = 0.0_dp
    if(crys%polar) then
       do m1=-ncell_g(1),ncell_g(1)
          do m2=-ncell_g(2),ncell_g(2)
             do m3=-ncell_g(3),ncell_g(3)
                g(1:3)=m1*self%cell_g(1,1:3)+&
                     m2*self%cell_g(2,1:3)+m3*self%cell_g(3,1:3)
                geg=dot_product(g(1:3),matmul(eps,g(1:3)))
                if(geg.gt.0.0d0.and.geg/alpha/4.0d0.lt.gmax) then
                   exp_g=exp(-geg/alpha/4.0d0)/geg
                   do iat=1,nat
                      zig(1:3)=matmul(g(1:3),zeff(iat,1:3,1:3))
                      auxi(1:3)=0.
                      do jat=1,nat
                         gr=dot_product(g(1:3),self%rr(iat,jat,1:3))
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
                            gr=dot_product(g(1:3),self%rr(iat,jat,1:3))
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
                                             zig(ipol)*zjg(jpol)*oneI*self%rr(iat,jat,i)-&
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
                   dyn(idim,jdim)=dyn(idim,jdim)/self%mm(iat,jat)
                   if(present(velocities)) then
                      ddyn(idim,jdim,1:3)=ddyn(idim,jdim,1:3)/self%mm(iat,jat)
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
             velocities(ik,i,:)=velocities(ik,i,:)/(2.0_dp*omegas(ik,i))
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
