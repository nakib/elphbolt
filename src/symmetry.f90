! Copyright 2020 elphbolt contributors.
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

module symmetry_module
  !! Module containing type and procedures related to crystal and
  !! Brillouin zone symmetries.

  use precision, only: r64, i64
  use misc, only: mux_vector, demux_mesh, demux_vector, &
       exit_with_message, subtitle, distribute_points, shrink
  use crystal_module, only : crystal

  use iso_c_binding
  use spglib_f08, only: SpglibDataset, spg_get_dataset, &
       SpglibMagneticDataset, spg_get_magnetic_dataset

  implicit none

  private
  public symmetry, find_equiv_map, find_irred_wedge, create_fbz2ibz_map, &
       fbz2ibz, symmetrize_3x3_tensor, symmetrize_3x3_tensor_noTR

  type symmetry
     !! Data and procedure related to symmetries.

     integer(i64) :: nsymm
     !! Number of spacegroup symmetries.
     integer(i64) :: nsymm_rot
     !! Number of rotations (including time-reversed ones if use_TR is true).
     logical :: use_TR = .true.
     !! Is time-reversal symmetry included in the rotations list?
     integer(i64), allocatable :: rotations_orig(:,:,:)
     !! Rotations without time-reversal, real space, crystal coordinates.
     real(r64), allocatable :: crotations_orig(:,:,:)
     !! Rotations without time-reversal, real space, Cartesian coordinates.
     real(r64), allocatable :: qrotations_orig(:,:,:)
     !! Rotations without time-reversal, reciprocal space, crystal coordinates.
     !And with time-reversal, if use_TR is true (time reversed sector is the 2nd half of last axis).
     !If use_TR is false, the following are identical to the above "_orig" arrays:
     integer(i64), allocatable :: rotations(:,:,:) 
     !! Rotations with time-reversal, real space, crystal coordinates.
     real(r64), allocatable :: crotations(:,:,:)
     !! Rotations with time-reversal, real space, Cartesian coordinates.
     real(r64), allocatable :: qrotations(:,:,:)
     !! Rotations with time-reversal, reciprocal space, crystal coordinates.
     logical, allocatable :: is_TR(:)
     !! Is the rotation (in the above list) a time-reversed one?
     character(len=10) :: international
     !! Spacegroup in Hermann–Mauguin (or international) notation.

   contains

     procedure :: calculate_symmetries

  end type symmetry

contains

  subroutine calculate_symmetries(self, crys, mesh, use_TR, magmoms)
    !! Subroutine to generate the symmetry related data for a given crystal.
    !!
    !! This subroutine closely follows parts of config.f90 of the ShengBTE code.
    !!
    !! If use_TR is true (default), the space group of the structure is used and
    !! each operation is doubled by its time-reversed partner.
    !! If use_TR is false and magmoms is not provided, only the space group
    !! operations are used.
    !! If use_TR is false and magmoms is provided, the magnetic space group
    !! from spglib is used. Each of its operations carries its own
    !! time-reversal flag: an operation with the flag set acts on wave vectors
    !! as k -> -R k and is stored with the negated rotation, just like the
    !! time-reversed partners above.

    class(symmetry), intent(out) :: self
    type(crystal), intent(in) :: crys
    integer(i64), intent(in) :: mesh(3)
    logical, intent(in), optional :: use_TR
    !! Include time-reversal symmetry? Default is true.
    real(r64), intent(in), optional :: magmoms(:, :)
    !! Cartesian axial vector (magnetic moment) on each basis atom. Shape (3, numatoms).
    !! Only used if use_TR is false.

    !Internal variables:
    integer(i64) :: i, ii, jj, kk, ll, info, nq, nlen
    integer(i64) :: P(3)
    integer(i64), allocatable :: rtmp(:,:,:), local_equiv_map(:,:)
    logical, allocatable :: valid(:), is_TR_tmp(:)
    real(r64), allocatable :: crtmp(:,:,:), qrtmp(:,:,:)
    real(r64), allocatable :: translations(:,:), ctranslations(:,:)
    real(r64) :: tmp1(3, 3), tmp2(3, 3), tmp3(3, 3)
    type(spglibdataset) :: symdataset
    type(SpglibMagneticDataset) :: magdataset
    logical :: use_magnetic
    logical, allocatable :: op_TR(:), TR_flag_spglib(:)
    integer(i64), allocatable :: order(:)
    integer(i64) :: iorder, nops
    real(r64), allocatable :: all_translations(:, :)
    integer(kind = C_INT) :: numatoms_cint
    integer(kind = C_INT) :: atomtypes_cint(crys%numatoms)

    !External procedures
    external :: dgesv

    call subtitle("Analyzing symmetry...")

    self%use_TR = .true.
    if(present(use_TR)) self%use_TR = use_TR

    use_magnetic = .false.
    if(present(magmoms)) then
       if(self%use_TR) then
          if(this_image() == 1) print*, &
               "Magnetic moments provided but use_TR is true. Ignoring the moments."
       else
          use_magnetic = .true.
       end if
    end if

    !Number of points in wave vector mesh
    nq = product(mesh)

    numatoms_cint = crys%numatoms
    atomtypes_cint = crys%atomtypes
    symdataset = spg_get_dataset(transpose(crys%lattvecs), crys%basis, &
         atomtypes_cint, numatoms_cint, 1.0e-5_r64)

    !Grab subset of symmetry info from spglib data set
    ! Number of symmetry operations
    nops = symdataset%n_operations

    !Operations of the magnetic space group, if requested
    allocate(TR_flag_spglib(nops))
    TR_flag_spglib = .false.
    if(use_magnetic) then
       if(size(magmoms, 1) /= 3 .or. size(magmoms, 2) /= crys%numatoms) &
            call exit_with_message("magmoms should have shape (3, numatoms). Exiting.")
       
       magdataset = spg_get_magnetic_dataset(transpose(crys%lattvecs), crys%basis, &
            atomtypes_cint, reshape(magmoms, [3*crys%numatoms]), 1_C_INT, &
            numatoms_cint, 1_C_INT, 1.0e-5_r64)

       nops = magdataset%n_operations

       deallocate(TR_flag_spglib)

       allocate(TR_flag_spglib(nops))
       TR_flag_spglib = (magdataset%time_reversals == 1)
    end if

    !Order the operations so that the unitary ones come before the time-reversed
    !ones. order(i) is the index, in spglib's list, of the i-th operation in our list.
    allocate(order(nops))

    iorder = 0

    do i = 1, nops
       if(.not. TR_flag_spglib(i)) then
          iorder = iorder + 1
          order(iorder) = i
       end if
    end do

    do i = 1, nops
       if(TR_flag_spglib(i)) then
          iorder = iorder + 1
          order(iorder) = i
       end if
    end do

    !op_TR(i) is the time-reversal flag of the i-th operation in our list.
    allocate(op_TR(nops))
    do i = 1, nops
       op_TR(i) = TR_flag_spglib(order(i))
    end do

    self%nsymm = nops

    ! Double the above to take time reversal symetry (TRS) into account, if requested.
    if(self%use_TR) then
       self%nsymm_rot = 2*self%nsymm
    else
       self%nsymm_rot = self%nsymm
    end if

    ! Allocate internal symmetry related data
    allocate(self%rotations(3,3,self%nsymm_rot),self%crotations(3,3,self%nsymm_rot),&
         self%qrotations(3,3,self%nsymm_rot),self%rotations_orig(3,3,self%nsymm),&
         self%crotations_orig(3,3,self%nsymm),self%qrotations_orig(3,3,self%nsymm),&
         translations(3,self%nsymm),ctranslations(3,self%nsymm))

    ! Mind the difference in convention between spglib and this code
    do i = 1, self%nsymm
       if(use_magnetic) then
          self%rotations_orig(:, :, i) = transpose(magdataset%rotations(:, :, order(i)))
       else
          self%rotations_orig(:, :, i) = transpose(symdataset%rotations(:, :, order(i)))
       end if
    end do

    ! Set just the front half
    self%rotations(:, :, 1:self%nsymm) = self%rotations_orig

    ! Set the translations
    if(use_magnetic) then
       all_translations = magdataset%translations
    else
       all_translations = symdataset%translations
    end if
    translations = all_translations(:, order)

    ! Set the international symbol
    self%international = trim(adjustl(symdataset%international_symbol))
    !!

    if(this_image() == 1) then
       write(*, "(A, A)") "Crystal symmetry group = ", self%international
       write(*, "(A, I3)") "Spacegroup number = ", symdataset%spacegroup_number
       write(*, "(A, I5)") "Number of crystal symmetries (without time-reversal) = ", self%nsymm
       write(*, "(A, L2)") "Time-reversal symmetry included = ", self%use_TR
       if(use_magnetic) write(*, "(A, I5, A, I5)") &
            "Magnetic space group: number of operations = ", self%nsymm, &
            ", of which time-reversed = ", count(op_TR)
    end if

    !Set symmetry operations in Cartesian basis.
    ctranslations = matmul(crys%lattvecs, translations)

    do i = 1, self%nsymm
       tmp1 = transpose(crys%lattvecs)
       tmp2 = transpose(matmul(crys%lattvecs, self%rotations_orig(:, :, i)))
       call dgesv(3, 3, tmp1, 3, P, tmp2, 3, info)
       self%crotations_orig(:, :, i) = transpose(tmp2)
    end do
    self%crotations(:, :, 1:self%nsymm) = self%crotations_orig
    !!

    !Transform the rotation matrices to the reciprocal-space basis.
    do i = 1, self%nsymm
       tmp1 = matmul(transpose(crys%lattvecs), crys%lattvecs)
       tmp2 = transpose(self%rotations_orig(:, :, i))
       tmp3 = tmp1
       call dgesv(3, 3, tmp1, 3, P, tmp2, 3, info)
       self%qrotations_orig(:, :, i) = transpose(matmul(tmp2, tmp3))
    end do
    self%qrotations(:,:,1:self%nsymm) = self%qrotations_orig

    !Fill the second half of the rotation matrix list using TRS.
    if(self%use_TR) then
       self%rotations(:,:,self%nsymm+1:2*self%nsymm) = -self%rotations_orig(:,:,1:self%nsymm)
       self%qrotations(:,:,self%nsymm+1:2*self%nsymm) = -self%qrotations_orig(:,:,1:self%nsymm)
       self%crotations(:,:,self%nsymm+1:2*self%nsymm) = -self%crotations_orig(:,:,1:self%nsymm)
    end if

    !Magnetic space group operations that come with time reversal act as -R on k.
    if(use_magnetic) then
       do i = 1, self%nsymm
          if(op_TR(i)) then
             self%rotations(:, :, i) = -self%rotations_orig(:, :, i)
             self%qrotations(:, :, i) = -self%qrotations_orig(:, :, i)
             self%crotations(:, :, i) = -self%crotations_orig(:, :, i)
          end if
       end do
    end if

    !Flag the time-reversed rotations.
    allocate(self%is_TR(self%nsymm_rot))
    self%is_TR = .false.
    if(self%use_TR) self%is_TR(self%nsymm + 1:2*self%nsymm) = .true.
    if(use_magnetic) self%is_TR = op_TR

    !Find rotations that are either duplicated or incompatible with mesh.
    allocate(local_equiv_map(self%nsymm_rot,nq))
    call find_equiv_map(self%nsymm_rot,local_equiv_map,mesh,self%qrotations)
    allocate(valid(self%nsymm_rot))
    valid = .true.
    jj = 0
    do ii = 1,self%nsymm_rot
       if(valid(ii) .and. any(local_equiv_map(ii,:) == -1)) then
          valid(ii) = .false.
          jj = jj + 1
       end if
    end do
    if(this_image() == 1 .and. jj /= 0) then
       print*, jj, &
            "Rotations are incompatible with the wave vector mesh and will be discarded."
    end if
    ll = 0
    do ii = 2,self%nsymm_rot
       do i = 1,ii - 1
          if(.not. valid(i)) cycle
          if(all(self%rotations(:,:,ii) == self%rotations(:,:,i))) then
             valid(ii) = .false.
             ll = ll + 1
             exit
          end if
       end do
    end do
    if(this_image() == 1 .and. ll /= 0) then
       write(*, "(A, I5)") "Number of duplicated rotations to be discarded = ", ll
    end if

    !Filter out those rotations through a series of move_alloc calls.
    !Arrays to take into account: rotations,crotations,qrotations.
    if(ll + jj /= 0) then
       allocate(rtmp(3,3,self%nsymm_rot - ll - jj))
       allocate(crtmp(3,3,self%nsymm_rot - ll - jj))
       allocate(qrtmp(3,3,self%nsymm_rot - ll - jj))
       allocate(is_TR_tmp(self%nsymm_rot - ll - jj))
       kk = 0
       do ii = 1,self%nsymm_rot
          if(valid(ii)) then
             kk = kk + 1
             rtmp(:,:,kk) = self%rotations(:,:,ii)
             crtmp(:,:,kk) = self%crotations(:,:,ii)
             qrtmp(:,:,kk) = self%qrotations(:,:,ii)
             is_TR_tmp(kk) = self%is_TR(ii)
          end if
       end do
       self%nsymm_rot = self%nsymm_rot - ll - jj
       call move_alloc(rtmp,self%rotations)
       call move_alloc(crtmp,self%crotations)
       call move_alloc(qrtmp,self%qrotations)
       call move_alloc(is_TR_tmp,self%is_TR)
    end if
  end subroutine calculate_symmetries

  subroutine find_star(q_in,q_out,mesh,qrotations)
    !! Compute all images of a wave vector (crystal coords.) under the
    !! rotational symmetry operations.

    integer(i64), intent(in) :: q_in(3), mesh(3)
    real(r64), intent(in) :: qrotations(:,:,:)
    real(r64), intent(out) :: q_out(:,:)

    integer(i64) :: ii, nsymm_rot

    nsymm_rot = size(qrotations(1,1,:))

    do ii = 1, nsymm_rot
       q_out(:, ii) = mesh*matmul(qrotations(:, :, ii),dble(q_in)/mesh)
    end do
  end subroutine find_star

  subroutine find_equiv_map(nsymm_rot,equiv_map,mesh,qrotations,indexlist)
    !! Subroutine to create the map of equivalent wave vectors.

    integer(i64), intent(in) :: nsymm_rot, mesh(3)
    real(r64), intent(in) :: qrotations(:,:,:)
    integer(i64), optional, intent(in) :: indexlist(:)
    integer(i64), intent(out) :: equiv_map(:,:)

    !Local variables
    integer(i64) :: nmesh, chunk, counter, im, num_active_images
    integer(i64), allocatable :: index_mesh(:,:)
    integer(i64) :: i, isym, ivec(3), base
    real(r64) :: vec(3), vec_star(3, nsymm_rot), dnrm2
    integer(i64), allocatable :: start[:], end[:], equiv_map_chunk(:,:)[:]

    if(present(indexlist)) then
       nmesh = size(indexlist)
    else
       nmesh = product(mesh)
    end if

    allocate(index_mesh(3,nmesh))

    !Create mesh of demuxed 0-based indices.
    base = 0
    if(present(indexlist)) then
       call demux_mesh(index_mesh, mesh, base, indexlist)
    else
       call demux_mesh(index_mesh, mesh, base)
    end if

    !Allocate start and end coarrays
    allocate(start[*], end[*])

    !Divide wave vectors among images
    call distribute_points(nmesh, chunk, start, end, num_active_images)

    !Allocate small work variable chunk for each image
    allocate(equiv_map_chunk(nsymm_rot, chunk)[*])

    !Only work with the active images
    if(this_image() <= num_active_images) then
       counter = 0
       do i = start, end !Run over total number of wave vectors.
          !Increase counter
          counter = counter + 1
          call find_star(index_mesh(:,i),vec_star,mesh,qrotations) !Find star of wave vector.
          do isym = 1, nsymm_rot !Run over all rotational symmetries of system.
             vec = vec_star(:,isym) !Pick image.
             ivec = nint(vec) !Snap to nearest integer grid.
             !Check norm and save mapping:
             if(dnrm2(3,abs(vec - dble(ivec)),1) >= 1e-2_r64) then
                equiv_map_chunk(isym, counter) = -1
             else
                equiv_map_chunk(isym, counter) = mux_vector(modulo(ivec,mesh),mesh,base)
             end if
          end do
       end do
    end if

    !Gather equiv_map_chunks in equiv_map and broadcast to all
    sync all
    if(this_image() == 1) then
       do im = 1, num_active_images
          equiv_map(:, start[im]:end[im]) = equiv_map_chunk(:,:)[im]
       end do
    end if
    sync all
    call co_broadcast(equiv_map, 1)
    sync all
  end subroutine find_equiv_map

  subroutine find_irred_wedge(mesh,nwavevecs_irred,wavevecs_irred, &
       indexlist_irred,nequivalent,nsymm_rot,qrotations,ibz2fbz_map,equivalence_map,blocks,indexlist)
    !! Find the irreducible wedge of the FBZ and other quantities.
    !! Wedge finding algorithm is inspired by ShengBTE.
    !!
    !! mesh is the array of number of points along the reciprocal lattice vectors
    !! nwavevecs_irred is the number of irreducible vectors
    !! wavevecs_irred are the irreducible vectors
    !! indexlist_irred is the list containing the muxed indices of the irreducible wave vectors
    !! nequivalent lists the number of equivalent points of each irreducible point
    !! ibz2fbz_map contains the map from an irreducible BZ (IBZ) vector to its FBZ images
    !!   The third axis contains the pair (symmetry index, image)
    !! equivalence_map is a map of the equivalent points under rotations
    !! blocks chooses whether the FBZ is energy restricted or not
    !! indexlist is the sorted list of indices of the wavevectors
    !!   in the energy restricted FBZ which must be present if blocks is true

    integer(i64), intent(in) :: mesh(3)
    logical, intent(in) :: blocks
    integer(i64), intent(in) :: nsymm_rot
    real(r64), intent(in) :: qrotations(:,:,:)
    integer(i64), optional, intent(in) :: indexlist(:)
    integer(i64), intent(out) :: nwavevecs_irred
    integer(i64), allocatable, intent(out) :: indexlist_irred(:), &
         nequivalent(:), ibz2fbz_map(:,:,:), equivalence_map(:,:)
    real(r64), allocatable, intent(out) :: wavevecs_irred(:,:)

    !Local variables
    integer(i64) :: nwavevecs, i, imux, s, image, imagelist(nsymm_rot), &
         nrunninglist, counter, ijk(3), aux, num_active_images, chunk, check
    integer(i64), allocatable :: runninglist(:), ibz2fbz_map_tmp(:,:,:), &
         start[:], end[:]
    logical :: in_list

    if(blocks .and. .not. present(indexlist)) &
         call exit_with_message("If blocks is true then indexlist must be present")

    if(blocks) then
       nwavevecs = size(indexlist)
    else
       nwavevecs = product(mesh)
    end if

    allocate(equivalence_map(nsymm_rot, nwavevecs))

    if(blocks) then
       call find_equiv_map(nsymm_rot, equivalence_map, mesh, qrotations, indexlist)
    else
       call find_equiv_map(nsymm_rot, equivalence_map, mesh, qrotations)
    end if

    allocate(indexlist_irred(nwavevecs), nequivalent(nwavevecs), &
         runninglist(nwavevecs), ibz2fbz_map_tmp(nsymm_rot, nwavevecs, 2))

    !Allocate coarrays
    allocate(start[*], end[*])

    runninglist = 0
    nrunninglist = 0
    nwavevecs_irred = 0
    nequivalent = 0
    counter = 0

    do i = 1, nwavevecs !Take a point from the FBZ
       !Get the muxed index of the wave vector
       if(blocks) then
          imux = indexlist(i)
       else
          imux = i
       end if

       !Check if point is not already in the running list of points
       in_list = .false.
       if(nrunninglist > 0) then
          !Divide wave vectors among images
          call distribute_points(nrunninglist, chunk, start, end, num_active_images)

          check = 0
          if(start > 0) then
             if(any(runninglist(start:end) == imux)) check = 1
          end if
          sync all
          call co_sum(check)
          sync all

          if(check > 0) in_list = .true.
       end if
       sync all

       if(.not. in_list) then
          !Increment irreducible point counter
          nwavevecs_irred = nwavevecs_irred + 1
          !Save point to irreducible wedge list
          indexlist_irred(nwavevecs_irred) = imux

          !Generate images of this irreducible point
          do s = 1, nsymm_rot !Take a rotation
             image = equivalence_map(s, i) !This is the image

             !Check if image is not already in the list of images
             if(.not. any(imagelist(1:nequivalent(nwavevecs_irred)) == image)) then
                !Increment equivalent image counter
                nequivalent(nwavevecs_irred) = & 
                     nequivalent(nwavevecs_irred) + 1
                aux = nequivalent(nwavevecs_irred)
                !Save image to list of images and running list of
                !points that have already been considered
                imagelist(aux) = image
                nrunninglist = nrunninglist + 1
                runninglist(nrunninglist) = image
                !Save mapping of the irreducible point to its FBZ image
                ibz2fbz_map_tmp(aux, &
                     nwavevecs_irred, :) = [s, image]
             end if
          end do
          counter = counter + nequivalent(nwavevecs_irred)
       end if
    end do

    !Check for error
    if(nwavevecs /= counter) call exit_with_message("Severe error: Could not find irreducible wedge.")

    if(this_image() == 1) then
       write(*, "(A, I10)") " Number of FBZ wave vectors = ", counter
       write(*, "(A, I10)") " Number IBZ wave vectors = ", nwavevecs_irred
    end if

    !Deallocate some internal data
    deallocate(runninglist)

    !Shrink indexlist_irred to the correct size
    call shrink(indexlist_irred, nwavevecs_irred)

    !Shrink nequivalent to the correct size
    call shrink(nequivalent, nwavevecs_irred)

    allocate(ibz2fbz_map(nsymm_rot, nwavevecs_irred, 2))
    ibz2fbz_map(:, 1:nwavevecs_irred, :) = ibz2fbz_map_tmp(:, 1:nwavevecs_irred, :)
    deallocate(ibz2fbz_map_tmp)

    !Create crystal coords IBZ wave vectors
    allocate(wavevecs_irred(nwavevecs_irred,3))
    do i = 1, nwavevecs_irred !run over total number of vectors
       imux = indexlist_irred(i)
       call demux_vector(imux, ijk, mesh, 0_i64) !get 0-based (i,j,k) indices

       wavevecs_irred(i,:) = dble(ijk)/mesh !wave vectors in crystal coordinates
    end do
  end subroutine find_irred_wedge

  function fbz2ibz(iwvmux,nwv_irred,nequiv,ibz2fbz_map)
    !! Find index in IBZ blocks list for a given FBZ blocks muxed vector index

    integer(i64), intent(in) :: iwvmux, nwv_irred, nequiv(nwv_irred), ibz2fbz_map(:,:,:)
    integer(i64) :: i, l, il, fbz2ibz

    fbz2ibz = -1

    !Sum over all ibz points
    do i = 1,nwv_irred !an irreducible point
       do l = 1,nequiv(i) !number of equivalent points of i
          !Get (i, l) -> il, the muxed vector index of image
          il = ibz2fbz_map(l,i,2)

          if(il == iwvmux) then 
             fbz2ibz = i
             exit
          end if

       end do
    end do

    if(fbz2ibz == -1) then 
       print*, 'Error in fbz2ibz for input index = ', iwvmux
    end if
  end function fbz2ibz

  subroutine create_fbz2ibz_map(fbz2ibz_map, nwv, nwv_irred, indexlist, nequiv, ibz2fbz_map)
    !! Subroutine to create map of FBZ blocks to IBZ blocks

    integer(i64), intent(in) :: nwv, nwv_irred, indexlist(nwv), &
         nequiv(nwv_irred),ibz2fbz_map(:,:,:)
    integer(i64), intent(out), allocatable :: fbz2ibz_map(:)

    integer(i64) :: iwv

    allocate(fbz2ibz_map(nwv))

    do iwv = 1,nwv !Run over all wave vectors
       fbz2ibz_map(iwv) = fbz2ibz(indexlist(iwv),nwv_irred,nequiv,ibz2fbz_map)
    end do
  end subroutine create_fbz2ibz_map

  subroutine symmetrize_3x3_tensor(tensor, crotations)
    !! Symmetrize a 3x3 tensor.

    real(r64), intent(inout) :: tensor(3,3)
    real(r64), intent(in) :: crotations(:,:,:)
    integer(i64) :: irot, nrots
    real(r64) :: aux(3,3)

    nrots = size(crotations(1, 1, :))

    aux(:,:) = 0.0_r64
    do irot = 1, nrots
       aux(:,:) = aux(:,:) + matmul(crotations(:, :, irot),&
            matmul(tensor, transpose(crotations(:, :, irot))))
    end do

    tensor(:,:) = aux(:,:)/nrots
  end subroutine symmetrize_3x3_tensor

  subroutine symmetrize_3x3_tensor_noTR(tensor, crotations, Bfield)
    !! Symmetrize a 3x3 tensor in the presence of a B-field.
    !! Note: Only for B-field of the form [0 0 Bz]
    !! The crotations passed in must not contain time-reversed operations,
    !! e.g. crotations_orig of the symmetry object.

    real(r64), intent(inout) :: tensor(3, 3)
    real(r64), intent(in) :: crotations(:, :, :)
    real(r64), intent(in) :: Bfield(3) !At the moment not doing anything with it.
    integer :: irot, nrots, i
    real(r64) :: aux(3,3)

    nrots = size(crotations(1, 1, :))

    aux = 0.0_r64
    do irot = 1, nrots
       aux(:, :) = aux(:, :) + matmul(crotations(:, :, irot),&
            matmul(tensor, transpose(crotations(:, :, irot))))
    end do

    !Symmetrize along z:
    tensor(3, 1:2) = aux(3, 1:2)/nrots
    tensor(1:2, 3) = aux(1:2, 3)/nrots

    !Symmetrize diagonals
    do i = 1, 2
       tensor(i, i) = aux(i, i)/nrots
    end do

    !Enforce Onsager in the simplest way. Will refine later.
    tensor(2, 1) = -tensor(1, 2)
  end subroutine symmetrize_3x3_tensor_noTR
end module symmetry_module
