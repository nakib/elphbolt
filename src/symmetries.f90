module symmetry_module
  !! Module containing type and procedures related crystal and
  !! Brillouin zone symmetries.

  use params, only: dp, k4
  use misc, only: mux_vector, demux_mesh
  use crystal_module !, only :: crystal
  use spglib_wrapper, only: get_operations, get_cartesian_operations, get_num_operations
  
  implicit none

  private
  public symmetry

  type symmetry
     !! Data and procedure related to symmetries.
     
     integer(k4) :: nsymm
     !! Number of spacegroup symmetries.
     integer(k4) :: nsymm_rot
     !! Number of rotations.
     integer(k4), allocatable :: rotations_orig(:,:,:)
     !! Rotations without time-reversal, real space, crystal coordinates.
     real(dp), allocatable :: crotations_orig(:,:,:)
     !! Rotations without time-reversal, real space, Cartesian coordinates.
     real(dp), allocatable :: qrotations_orig(:,:,:)
     !! Rotations without time-reversal, reciprocal space, crystal coordinates.
     !And with time-reversal (time reversed sector is the 2nd half of last axis):
     integer(k4), allocatable :: rotations(:,:,:) 
     !! Rotations with time-reversal, real space, crystal coordinates.
     real(dp), allocatable :: crotations(:,:,:)
     !! Rotations with time-reversal, real space, Cartesian coordinates.
     real(dp), allocatable :: qrotations(:,:,:)
     !! Rotations with time-reversal, reciprocal space, crystal coordinates.
     integer(k4), allocatable :: equiv_map(:,:)
     !! Map of equivalent points under rotations.
     !! Axis 1 runs over rotations.
     !! Axis 2 runs over wave vectors (full Brillouin zone).
     integer(k4), allocatable :: equiv_map_blocks(:,:)
     !! Map of equivalent points under rotations.
     !! Axis 1 runs over rotations.
     !! Axis 2 runs over wave vectors (energy windowed blocks of full Brillouin zone).
     character(len=10) :: international
     !! Spacegroup in Hermann–Mauguin (or international) notation.

   contains

     procedure :: calculate_symmetries
     
  end type symmetry

contains

  subroutine calculate_symmetries(s, c, mesh)
    !! Subroutine to generate the symmetry related data for a given crystal.
    !!
    !! This subroutine closely follows parts of config.f90 of the ShengBTE code.

    class(symmetry), intent(out) :: s
    type(crystal), intent(in) :: c
    integer(k4), intent(in) :: mesh(3)

    !Internal variables:
    integer(k4) :: i, j, k, ii, jj, kk, ll, info, nq
    integer(k4) :: P(3)
    integer(k4), allocatable :: rtmp(:,:,:), local_equiv_map(:,:)
    logical, allocatable :: valid(:)
    real(dp), allocatable :: crtmp(:,:,:), qrtmp(:,:,:)
    real(dp), allocatable :: translations(:,:), ctranslations(:,:)
    real(dp) :: tmp1(3,3), tmp2(3,3), tmp3(3,3)

    if(this_image() == 1) print*, 'Analyzing symmetry...'
    
    !Number of points in wave vector mesh
    nq = product(mesh)
    
    !Number of crystal symmetries.
    s%nsymm = get_num_operations(c%lattvecs,c%numatoms,c%atomtypes,c%basis)
    !Double the above to take time reversal symetry (TRS) into account.
    s%nsymm_rot = 2*s%nsymm

    allocate(s%rotations(3,3,s%nsymm_rot),s%crotations(3,3,s%nsymm_rot),&
         s%qrotations(3,3,s%nsymm_rot),s%rotations_orig(3,3,s%nsymm),&
         s%crotations_orig(3,3,s%nsymm),s%qrotations_orig(3,3,s%nsymm),&
         translations(3,s%nsymm),ctranslations(3,s%nsymm))

    !Get symmetry operations.
    call get_operations(c%lattvecs,c%numatoms,c%atomtypes,&
         c%basis,s%nsymm,s%rotations_orig,translations,s%international)
    s%rotations(:,:,1:s%nsymm) = s%rotations_orig

    if(this_image() == 1) then
       print*, "Crystal symmetry group = ", trim(s%international)
       print*, "Number of crystal symmetries (without time-reversal) = ", s%nsymm
    end if

    !Get symmertry operations in Cartesian basis.
    call get_cartesian_operations(c%lattvecs,s%nsymm,&
         s%rotations_orig,translations,&
         s%crotations_orig,ctranslations)
    s%crotations(:,:,1:s%nsymm) = s%crotations_orig

    !Transform the rotation matrices to the reciprocal-space basis.
    do i = 1,s%nsymm
       tmp1 = matmul(transpose(c%lattvecs),c%lattvecs)
       tmp2 = transpose(s%rotations_orig(:, :, i))
       tmp3 = tmp1
       call dgesv(3,3,tmp1,3,P,tmp2,3,info)
       s%qrotations_orig(:,:,i) = transpose(matmul(tmp2,tmp3))
    end do
    s%qrotations(:,:,1:s%nsymm) = s%qrotations_orig

    !Fill the second half of the rotation matrix list using TRS.
    s%rotations(:,:,s%nsymm+1:2*s%nsymm) = -s%rotations_orig(:,:,1:s%nsymm)
    s%qrotations(:,:,s%nsymm+1:2*s%nsymm) = -s%qrotations_orig(:,:,1:s%nsymm)
    s%crotations(:,:,s%nsymm+1:2*s%nsymm) = -s%crotations_orig(:,:,1:s%nsymm)

    !Find rotations that are either duplicated or incompatible with mesh.
    allocate(local_equiv_map(s%nsymm_rot,nq))
    call find_equiv_map(s%nsymm_rot,local_equiv_map,mesh,s%qrotations)
    allocate(valid(s%nsymm_rot))
    valid = .true.
    jj = 0
    do ii = 1,s%nsymm_rot
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
    do ii = 2,s%nsymm_rot
       do i = 1,ii - 1
          if(.not. valid(i)) cycle
          if(all(s%rotations(:,:,ii) == s%rotations(:,:,i))) then
             valid(ii) = .false.
             ll = ll + 1
             exit
          end if
       end do
    end do
    if(this_image() == 1 .and. ll == 0) then
       print*, "Number of duplicated rotations to be discarded = ", ll
    end if

    !Filter out those rotations through a series of move_alloc calls.
    !Arrays to take into account: rotations,crotations,qrotations.
    if(ll + jj /= 0) then
       allocate(rtmp(3,3,s%nsymm_rot - ll - jj))
       allocate(crtmp(3,3,s%nsymm_rot - ll - jj))
       allocate(qrtmp(3,3,s%nsymm_rot - ll - jj))
       kk = 0
       do ii = 1,s%nsymm_rot
          if(valid(ii)) then
             kk = kk + 1
             rtmp(:,:,kk) = s%rotations(:,:,ii)
             crtmp(:,:,kk) = s%crotations(:,:,ii)
             qrtmp(:,:,kk) = s%qrotations(:,:,ii)
          end if
       end do
       s%nsymm_rot = s%nsymm_rot - ll - jj
       call move_alloc(rtmp,s%rotations)
       call move_alloc(crtmp,s%crotations)
       call move_alloc(qrtmp,s%qrotations)
    end if
  end subroutine calculate_symmetries

  subroutine find_star(q_in,q_out,mesh,qrotations)
    !! Compute all images of a wave vector (crystal coords.) under the
    !! rotational symmetry operations.

    integer(k4), intent(in) :: q_in(3), mesh(3)
    real(dp), intent(in) :: qrotations(:,:,:)
    real(dp), intent(out) :: q_out(:,:)

    integer(k4) :: ii, nsymm_rot

    nsymm_rot = size(qrotations(1,1,:))

    do ii = 1, nsymm_rot
       q_out(:, ii) = mesh*matmul(qrotations(:, :, ii),dble(q_in)/mesh)
    end do
  end subroutine find_star

  subroutine find_equiv_map(nsymm_rot,equiv_map,mesh,qrotations,indexlist)
    !! Subroutine to create the map of equivalent wave vectors.

    integer(k4), intent(in) :: nsymm_rot, mesh(3)
    real(dp), intent(in) :: qrotations(:,:,:)
    integer(k4), optional, intent(in) :: indexlist(:)
    integer(k4), intent(out) :: equiv_map(:,:)

    integer(k4) :: nmesh
    integer(k4), allocatable :: index_mesh(:,:)
    integer(k4) :: i, isym, ivec(3), base
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

    deallocate(index_mesh)
  end subroutine find_equiv_map
  
end module symmetry_module
