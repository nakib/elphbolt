module config
  !! Module containing crystal information and subroutines related to the calculation.

  use params, only: dp, k4, twopi
  use derived_types, only: crystal_data, reciprocal_lattice_data, symmetry_data, electron_data
  use misc, only: cross_product, mux_vector, demux_mesh
  use spglib_wrapper, only: get_operations, get_cartesian_operations, get_num_operations
  
  implicit none

  private
  public crystal_data, reciprocal_lattice_data, initialize_system
    
contains

  subroutine initialize_system(crys,reclat,sym,el)
    !! Subroutine to read input.nml file and set up the calculation environment.
    !! Outputs crystal of type crystal_data and reciprocal_lattice of type
    !! reciprocal_lattice_date.

    type(crystal_data), intent(out) :: crys
    type(reciprocal_lattice_data), intent(out) :: reclat
    type(symmetry_data), intent(out) :: sym
    type(electron_data), intent(out) :: el
    !type(phonon_data), intent(out) :: ph

    !Local variables:
    integer(k4) :: i, j, k, numelements, numatoms, qmesh(3)
    integer(k4), allocatable :: atomtypes(:)
    character(len=3), allocatable :: elements(:)
    real(dp), allocatable :: masses(:), basis(:,:), basis_cart(:,:)
    real(dp) :: lattvecs(3,3)

    !Namelists:
    namelist /allocations/ numelements, numatoms
    namelist /crystal_info/ elements, atomtypes, basis, lattvecs
    namelist /numerics/ qmesh

    if(this_image() == 1) print*, 'Reading input.nml...'
    
    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Set values from input:
    ! Read allocations.
    read(1, nml = allocations)
    if(numelements < 1 .or. numatoms < 1 .or. numatoms > numelements) then
       !TODO Catch error
       print*, 'Bad input(s) in allocations.'
    end if

    crys%numelements = numelements
    crys%numatoms = numatoms
        
    allocate(elements(numelements), atomtypes(numatoms), &
         masses(numatoms), basis(3,numatoms))
    
    ! Read crystal_info.
    read(1, nml = crystal_info)
    if(any(atomtypes < 1) .or. any(masses < 0)) then
       !TODO Catch error
       print*, 'Bad input(s) in crystal_info.'
    end if
    
    allocate(crys%elements(numelements), crys%atomtypes(numatoms), &
         crys%masses(numatoms), crys%basis(3,numatoms), &
         crys%basis_cart(3,numatoms))
    crys%elements(:) = elements(:)
    crys%atomtypes(:) = atomtypes(:)
    crys%masses(:) = masses(:)
    crys%lattvecs(:,:) = lattvecs(:,:)
    crys%basis(:,:) = basis(:,:)
    crys%basis_cart(:,:) = matmul(lattvecs,basis)

    ! Calculate reciprocal lattice vectors and real and reciprocal cell volumes.
    do i = 1,3
       j = mod(i,3) + 1
       k = mod(j,3) + 1
       reclat%reclattvecs(:,i) = &
            cross_product(crys%lattvecs(:,j), crys%lattvecs(:,k))
    end do
    crys%volume = abs(dot_product(crys%lattvecs(:,1), &
         reclat%reclattvecs(:,1)))
    reclat%volume_bz = twopi/crys%volume
    reclat%reclattvecs(:,:) = &
         reclat%volume_bz*reclat%reclattvecs(:,:)

    ! Read mesh information.
    qmesh = (/0, 0, 0/)
    read(1, nml = numerics)
    if(any(qmesh <= 0)) then
       !TODO Catch error
       print*, 'Bad input(s) in numerics.'
    end if
    el%coarse_mesh = qmesh !Coarse electron mesh is the same as the phonon mesh.
    
    !Close input file.
    close(1) 

    !Calculate all symmetry related information.
    call calculate_symmetries(crys,sym,qmesh)
    
    !!!
  end subroutine initialize_system

  subroutine calculate_symmetries(crys,sym,qmesh)
    !! Subroutine to generate the symmetry related data
    !! for a given crystal.
    !!
    !! This subroutine closely follows parts of config.f90 of the ShengBTE code.

    type(crystal_data), intent(in) :: crys
    integer(k4), intent(in) :: qmesh(3)
    type(symmetry_data), intent(out) :: sym
    
    !Some internal variables:
    integer(k4) :: nq, i, j, k, ii, jj, kk, ll, info
    integer(k4) :: P(3)
    integer(k4), allocatable :: rtmp(:,:,:), equiv_map(:,:)
    logical, allocatable :: valid(:)
    real(dp), allocatable :: crtmp(:,:,:), qrtmp(:,:,:)
    real(dp), allocatable :: translations(:,:), ctranslations(:,:)
    real(dp) :: tmp1(3,3), tmp2(3,3), tmp3(3,3)

    !Number of wave vectors.
    nq = product(qmesh)
    
    !Number of crystal symmetries.
    sym%nsymm = get_num_operations(crys%lattvecs,crys%numatoms,&
         crys%atomtypes,crys%basis)
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
       print*, "Symmetry group = ", trim(sym%international)
       write(*,*) "Number of symmetries (without time-reversal) = ", sym%nsymm
    end if

    !Get symmertry operations in Cartesian basis.
    call get_cartesian_operations(crys%lattvecs,sym%nsymm,&
         sym%rotations_orig,translations,&
         sym%crotations_orig,ctranslations)
    sym%crotations(:,:,1:sym%nsymm) = sym%crotations_orig

    !Transform the rotation matrices to the reciprocal-space basis.
    do i = 1, sym%nsymm
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
    allocate(equiv_map(sym%nsymm_rot,nq))
    call find_equiv_map(sym%nsymm_rot,equiv_map,qmesh,sym%qrotations)
    allocate(valid(sym%nsymm_rot))
    valid = .true.
    jj = 0
    do ii = 1,sym%nsymm_rot
       if(valid(ii) .and. any(equiv_map(ii,:) == -1)) then
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
       print*, ll, "Duplicated rotations will be discarded."
    end if

    ! Filter out those rotations through a series of move_alloc calls.
    ! Arrays to take into account: rotations,crotations,qrotations.
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
    !!
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
  
end module config
