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

module symmetry_module
  !! Module containing type and procedures related crystal and
  !! Brillouin zone symmetries.

  use params, only: dp, k8
  use misc, only: mux_vector, demux_mesh, demux_vector, &
       exit_with_message, subtitle, distribute_points
  use crystal_module !, only :: crystal
  use spglib_wrapper, only: get_operations, get_cartesian_operations, get_num_operations
  
  implicit none

  private
  public symmetry, find_equiv_map, find_irred_wedge, create_fbz2ibz_map, &
       fbz2ibz, symmetrize_3x3_tensor

  type symmetry
     !! Data and procedure related to symmetries.
     
     integer(k8) :: nsymm
     !! Number of spacegroup symmetries.
     integer(k8) :: nsymm_rot
     !! Number of rotations.
     integer(k8), allocatable :: rotations_orig(:,:,:)
     !! Rotations without time-reversal, real space, crystal coordinates.
     real(dp), allocatable :: crotations_orig(:,:,:)
     !! Rotations without time-reversal, real space, Cartesian coordinates.
     real(dp), allocatable :: qrotations_orig(:,:,:)
     !! Rotations without time-reversal, reciprocal space, crystal coordinates.
     !And with time-reversal (time reversed sector is the 2nd half of last axis):
     integer(k8), allocatable :: rotations(:,:,:) 
     !! Rotations with time-reversal, real space, crystal coordinates.
     real(dp), allocatable :: crotations(:,:,:)
     !! Rotations with time-reversal, real space, Cartesian coordinates.
     real(dp), allocatable :: qrotations(:,:,:)
     !! Rotations with time-reversal, reciprocal space, crystal coordinates.
     character(len=10) :: international
     !! Spacegroup in Hermann–Mauguin (or international) notation.

   contains

     procedure :: calculate_symmetries
     
  end type symmetry

contains

  subroutine calculate_symmetries(sym, crys, mesh)
    !! Subroutine to generate the symmetry related data for a given crystal.
    !!
    !! This subroutine closely follows parts of config.f90 of the ShengBTE code.

    class(symmetry), intent(out) :: sym
    type(crystal), intent(in) :: crys
    integer(k8), intent(in) :: mesh(3)

    !Internal variables:
    integer(k8) :: i, j, k, ii, jj, kk, ll, info, nq, nlen
    integer(k8) :: P(3)
    integer(k8), allocatable :: rtmp(:,:,:), local_equiv_map(:,:)
    logical, allocatable :: valid(:)
    real(dp), allocatable :: crtmp(:,:,:), qrtmp(:,:,:)
    real(dp), allocatable :: translations(:,:), ctranslations(:,:)
    real(dp) :: tmp1(3,3), tmp2(3,3), tmp3(3,3)

    call subtitle("Analyzing symmetry...")
    
    !Number of points in wave vector mesh
    nq = product(mesh)
    
    !Number of crystal symmetries.
    sym%nsymm = get_num_operations(crys%lattvecs,crys%numatoms,crys%atomtypes,crys%basis)
    !Double the above to take time reversal symetry (TRS) into account.
    sym%nsymm_rot = 2*sym%nsymm

    allocate(sym%rotations(3,3,sym%nsymm_rot),sym%crotations(3,3,sym%nsymm_rot),&
         sym%qrotations(3,3,sym%nsymm_rot),sym%rotations_orig(3,3,sym%nsymm),&
         sym%crotations_orig(3,3,sym%nsymm),sym%qrotations_orig(3,3,sym%nsymm),&
         translations(3,sym%nsymm),ctranslations(3,sym%nsymm))

    !Get symmetry operations.
    call get_operations(crys%lattvecs,crys%numatoms,crys%atomtypes,&
         crys%basis,sym%nsymm,sym%rotations_orig,translations,sym%international)
    sym%rotations(:,:,1:sym%nsymm) = sym%rotations_orig

    if(this_image() == 1) then
       !This is a hacky fix to the problem of a trailing binary character
       !printing that happens on some machines.
       nlen = len(trim(sym%international)) - 1
       write(*, "(A, A)") "Crystal symmetry group = ", sym%international(1:nlen)
       write(*, "(A, I5)") "Number of crystal symmetries (without time-reversal) = ", sym%nsymm
    end if

    !Get symmertry operations in Cartesian basis.
    call get_cartesian_operations(crys%lattvecs,sym%nsymm,&
         sym%rotations_orig,translations,&
         sym%crotations_orig,ctranslations)
    sym%crotations(:,:,1:sym%nsymm) = sym%crotations_orig

    !Transform the rotation matrices to the reciprocal-space basis.
    do i = 1,sym%nsymm
       tmp1 = matmul(transpose(crys%lattvecs),crys%lattvecs)
       tmp2 = transpose(sym%rotations_orig(:, :, i))
       tmp3 = tmp1
       call dgesv(3,3,tmp1,3,P,tmp2,3,info)
       sym%qrotations_orig(:,:,i) = transpose(matmul(tmp2,tmp3))
    end do
    sym%qrotations(:,:,1:sym%nsymm) = sym%qrotations_orig

    !Fill the second half of the rotation matrix list using TRS.
    sym%rotations(:,:,sym%nsymm+1:2*sym%nsymm) = -sym%rotations_orig(:,:,1:sym%nsymm)
    sym%qrotations(:,:,sym%nsymm+1:2*sym%nsymm) = -sym%qrotations_orig(:,:,1:sym%nsymm)
    sym%crotations(:,:,sym%nsymm+1:2*sym%nsymm) = -sym%crotations_orig(:,:,1:sym%nsymm)

    !Find rotations that are either duplicated or incompatible with mesh.
    allocate(local_equiv_map(sym%nsymm_rot,nq))
    call find_equiv_map(sym%nsymm_rot,local_equiv_map,mesh,sym%qrotations)
    allocate(valid(sym%nsymm_rot))
    valid = .true.
    jj = 0
    do ii = 1,sym%nsymm_rot
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
    do ii = 2,sym%nsymm_rot
       do i = 1,ii - 1
          if(.not. valid(i)) cycle
          if(all(sym%rotations(:,:,ii) == sym%rotations(:,:,i))) then
             valid(ii) = .false.
             ll = ll + 1
             exit
          end if
       end do
    end do
    if(this_image() == 1 .and. ll == 0) then
       write(*, "(A, I5)") "Number of duplicated rotations to be discarded = ", ll
    end if

    !Filter out those rotations through a series of move_alloc calls.
    !Arrays to take into account: rotations,crotations,qrotations.
    if(ll + jj /= 0) then
       allocate(rtmp(3,3,sym%nsymm_rot - ll - jj))
       allocate(crtmp(3,3,sym%nsymm_rot - ll - jj))
       allocate(qrtmp(3,3,sym%nsymm_rot - ll - jj))
       kk = 0
       do ii = 1,sym%nsymm_rot
          if(valid(ii)) then
             kk = kk + 1
             rtmp(:,:,kk) = sym%rotations(:,:,ii)
             crtmp(:,:,kk) = sym%crotations(:,:,ii)
             qrtmp(:,:,kk) = sym%qrotations(:,:,ii)
          end if
       end do
       sym%nsymm_rot = sym%nsymm_rot - ll - jj
       call move_alloc(rtmp,sym%rotations)
       call move_alloc(crtmp,sym%crotations)
       call move_alloc(qrtmp,sym%qrotations)
    end if
  end subroutine calculate_symmetries

  subroutine find_star(q_in,q_out,mesh,qrotations)
    !! Compute all images of a wave vector (crystal coords.) under the
    !! rotational symmetry operations.

    integer(k8), intent(in) :: q_in(3), mesh(3)
    real(dp), intent(in) :: qrotations(:,:,:)
    real(dp), intent(out) :: q_out(:,:)

    integer(k8) :: ii, nsymm_rot

    nsymm_rot = size(qrotations(1,1,:))

    do ii = 1, nsymm_rot
       q_out(:, ii) = mesh*matmul(qrotations(:, :, ii),dble(q_in)/mesh)
    end do
  end subroutine find_star

  subroutine find_equiv_map(nsymm_rot,equiv_map,mesh,qrotations,indexlist)
    !! Subroutine to create the map of equivalent wave vectors.

    integer(k8), intent(in) :: nsymm_rot, mesh(3)
    real(dp), intent(in) :: qrotations(:,:,:)
    integer(k8), optional, intent(in) :: indexlist(:)
    integer(k8), intent(out) :: equiv_map(:,:)

    integer(k8) :: nmesh
    integer(k8), allocatable :: index_mesh(:,:)
    integer(k8) :: i, isym, ivec(3), base
    real(dp) :: vec(3), vec_star(3, nsymm_rot), dnrm2

    if(present(indexlist)) then
       nmesh = size(indexlist)
    else
       nmesh = product(mesh)
    end if

    allocate(index_mesh(3,nmesh))

    !Create mesh of demuxed 0-based indices.
    base = 0
    if(present(indexlist)) then
       call demux_mesh(index_mesh,nmesh,mesh,base,indexlist)
    else
       call demux_mesh(index_mesh,nmesh,mesh,base)
    end if

    do i = 1,nmesh !Run over total number of wave vectors.
       call find_star(index_mesh(:,i),vec_star,mesh,qrotations) !Find star of wave vector.
       do isym = 1,nsymm_rot !Run over all rotational symmetries of system.
          vec = vec_star(:,isym) !Pick image.
          ivec = nint(vec) !Snap to nearest integer grid.
          !Check norm and save mapping:
          if(dnrm2(3,abs(vec - dble(ivec)),1) >= 1e-2_dp) then
             equiv_map(isym,i) = -1
          else
             equiv_map(isym,i) = mux_vector(modulo(ivec,mesh),mesh,base)
          end if
       end do
    end do
  end subroutine find_equiv_map

  subroutine find_equiv_map_parallel(nsymm_rot,equiv_map,mesh,qrotations,indexlist)
    !! Subroutine to create the map of equivalent wave vectors.

    integer(k8), intent(in) :: nsymm_rot, mesh(3)
    real(dp), intent(in) :: qrotations(:,:,:)
    integer(k8), optional, intent(in) :: indexlist(:)
    integer(k8), intent(out) :: equiv_map(:,:)

    integer(k8) :: nmesh, chunk, counter, im, num_active_images
    integer(k8), allocatable :: index_mesh(:,:), start[:], end[:]
    integer(k8) :: i, isym, ivec(3), base
    real(dp) :: vec(3), vec_star(3, nsymm_rot), dnrm2
    integer(k8), allocatable :: equiv_map_chunk(:,:)[:]

    if(present(indexlist)) then
       nmesh = size(indexlist)
    else
       nmesh = product(mesh)
    end if

    allocate(index_mesh(3,nmesh))

    !Create mesh of demuxed 0-based indices.
    base = 0
    if(present(indexlist)) then
       call demux_mesh(index_mesh,nmesh,mesh,base,indexlist)
    else
       call demux_mesh(index_mesh,nmesh,mesh,base)
    end if

    !Allocate start and end coarrays
    allocate(start[*], end[*])
    
    !Divide wave vectors among images
    call distribute_points(nmesh, chunk, start, end, num_active_images)

    !Allocate small work variable chunk for each image
    allocate(equiv_map_chunk(nsymm_rot, chunk)[*])

    counter = 0
    do i = start, end !Run over total number of wave vectors.
       !Increase counter
       counter = counter + 1
       call find_star(index_mesh(:,i),vec_star,mesh,qrotations) !Find star of wave vector.
       do isym = 1, nsymm_rot !Run over all rotational symmetries of system.
          vec = vec_star(:,isym) !Pick image.
          ivec = nint(vec) !Snap to nearest integer grid.
          !Check norm and save mapping:
          if(dnrm2(3,abs(vec - dble(ivec)),1) >= 1e-2_dp) then
             equiv_map_chunk(isym, counter) = -1
          else
             equiv_map_chunk(isym, counter) = mux_vector(modulo(ivec,mesh),mesh,base)
          end if
       end do
    end do
    
    !Collect equiv_map_chunks in equiv_map
    sync all
    do im = 1, num_active_images
       equiv_map(:, start[im]:end[im]) = equiv_map_chunk(:,:)[im]
    end do
    sync all
  end subroutine find_equiv_map_parallel

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

    integer(k8), intent(in) :: mesh(3)
    logical, intent(in) :: blocks
    integer(k8), intent(in) :: nsymm_rot
    real(dp), intent(in) :: qrotations(:,:,:)
    integer(k8), optional, intent(in) :: indexlist(:)
    integer(k8), intent(out) :: nwavevecs_irred
    integer(k8), allocatable, intent(out) :: indexlist_irred(:), &
         nequivalent(:), ibz2fbz_map(:,:,:), equivalence_map(:,:)
    real(dp), allocatable, intent(out) :: wavevecs_irred(:,:)

    !Local variables
    integer(k8) :: nwavevecs, i, imux, s, image, imagelist(nsymm_rot), &
         nrunninglist, counter, ijk(3)
    integer(k8), allocatable :: runninglist(:), &
         indexlist_irred_tmp(:), nequivalent_tmp(:), ibz2fbz_map_tmp(:,:,:)
    logical :: proceed

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
       !call find_equiv_map_parallel(nsymm_rot, equivalence_map, mesh, qrotations, indexlist)
    else
       call find_equiv_map(nsymm_rot, equivalence_map, mesh, qrotations)
       !call find_equiv_map_parallel(nsymm_rot, equivalence_map, mesh, qrotations)
    end if

    allocate(indexlist_irred_tmp(nwavevecs), nequivalent_tmp(nwavevecs), &
         runninglist(nwavevecs), ibz2fbz_map_tmp(nsymm_rot, nwavevecs, 2))

    nwavevecs_irred = 0
    nequivalent_tmp = 0
    nrunninglist = 0
    runninglist = 0
    counter = 0
    do i = 1,nwavevecs !Take a point from the FBZ
       !Get the muxed index of the wave vector
       if(blocks) then
          imux = indexlist(i)
       else
          imux = i
       end if

       !Check if point is not already in the running list of points
       proceed = .not. any(runninglist(1:nrunninglist) == imux)
       if(proceed) then
          !Increment irreducible point counter
          nwavevecs_irred = nwavevecs_irred + 1
          !Save point to irreducible wedge list
          indexlist_irred_tmp(nwavevecs_irred) = imux

          !Generate images of this irreducible point
          do s = 1,nsymm_rot !Take a rotation
             image = equivalence_map(s, i) !This is the image

             !Check if image is not already in the list of images
             if(.not. any(imagelist(1:nequivalent_tmp(nwavevecs_irred)) == image)) then
                !Increment equivalent image counter
                nequivalent_tmp(nwavevecs_irred) = & 
                     nequivalent_tmp(nwavevecs_irred) + 1
                !Save image to list of images and running list of
                !points that have already been considered
                imagelist(nequivalent_tmp(nwavevecs_irred)) = image
                nrunninglist = nrunninglist + 1
                runninglist(nrunninglist) = image
                !Save mapping of the irreducible point to its FBZ image
                ibz2fbz_map_tmp(nequivalent_tmp(nwavevecs_irred), &
                     nwavevecs_irred, :) = (/s, image/)
             end if
          end do
          counter = counter + nequivalent_tmp(nwavevecs_irred)
       end if
    end do

    !Check for error
    if(nwavevecs /= counter) call exit_with_message("Severe error: Could not find irreducible wedge.")

    if(this_image() == 1) write(*, "(A, I10)") " Number of FBZ wave vectors = ", counter
    if(this_image() == 1) write(*, "(A, I10)") " Number IBZ wave vectors = ", nwavevecs_irred

    !Deallocate some internal data
    deallocate(runninglist)

    !Copy the tmp data into (much) smaller sized global data holders
    allocate(indexlist_irred(nwavevecs_irred), nequivalent(nwavevecs_irred), &
         ibz2fbz_map(nsymm_rot, nwavevecs_irred, 2))
    indexlist_irred(1:nwavevecs_irred) = indexlist_irred_tmp(1:nwavevecs_irred)
    nequivalent(1:nwavevecs_irred) = nequivalent_tmp(1:nwavevecs_irred)
    ibz2fbz_map(:, 1:nwavevecs_irred, :) = ibz2fbz_map_tmp(:, 1:nwavevecs_irred, :)

    !Deallocate some internal data
    deallocate(indexlist_irred_tmp, nequivalent_tmp, ibz2fbz_map_tmp)

    !Create crystal coords IBZ wave vectors
    allocate(wavevecs_irred(nwavevecs_irred,3))
    do i = 1,nwavevecs_irred !run over total number of vectors
       imux = indexlist_irred(i)
       call demux_vector(imux, ijk, mesh, 0_k8) !get 0-based (i,j,k) indices

       wavevecs_irred(i,:) = dble(ijk)/mesh !wave vectors in crystal coordinates
    end do
  end subroutine find_irred_wedge

  function fbz2ibz(iwvmux,nwv_irred,nequiv,ibz2fbz_map)
    !! Find index in IBZ blocks list for a given FBZ blocks muxed vector index

    integer(k8), intent(in) :: iwvmux, nwv_irred, nequiv(nwv_irred), ibz2fbz_map(:,:,:)
    integer(k8) :: i, l, il, fbz2ibz

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

    integer(k8), intent(in) :: nwv, nwv_irred, indexlist(nwv), &
         nequiv(nwv_irred),ibz2fbz_map(:,:,:)
    integer(k8), intent(out), allocatable :: fbz2ibz_map(:)

    integer(k8) :: iwv

    allocate(fbz2ibz_map(nwv))

    do iwv = 1,nwv !Run over all wave vectors
       fbz2ibz_map(iwv) = fbz2ibz(indexlist(iwv),nwv_irred,nequiv,ibz2fbz_map)
    end do
  end subroutine create_fbz2ibz_map
  
  subroutine symmetrize_3x3_tensor(tensor, crotations)
    !! Symmetrize a 3x3 tensor.

    real(dp), intent(inout) :: tensor(3,3)
    real(dp), intent(in) :: crotations(:,:,:)
    integer(k8) :: irot, nrots
    real(dp) :: aux(3,3)

    nrots = size(crotations(1, 1, :))
    
    aux(:,:) = 0.0_dp
    do irot = 1, nrots
       aux(:,:) = aux(:,:) + matmul(crotations(:, :, irot),&
            matmul(tensor, transpose(crotations(:, :, irot))))
    end do

    tensor(:,:) = aux(:,:)/nrots
  end subroutine symmetrize_3x3_tensor
  
end module symmetry_module
