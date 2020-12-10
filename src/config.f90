module config
  !! Module containing crystal information and subroutines related to the calculation.

  use data
  use params, only: dp, k4, twopi
  use misc, only: cross_product, mux_vector, demux_mesh, exit_with_message
  use spglib_wrapper, only: get_operations, get_cartesian_operations, get_num_operations
  use wannier, only: read_EPW_Wannier
  
  implicit none

  private
  public initialize_system, find_equiv_map
    
contains

  subroutine initialize_system
    !! Subroutine to read input.nml file and set up the calculation environment.
    
    !Local variables:
    integer(k4) :: i, j, k

    !Namelists:
    namelist /allocations/ numelements, numatoms
    namelist /crystal_info/ elements, atomtypes, basis, lattvecs, polar, &
         born, epsilon
    namelist /numerics/ qmesh, mesh_ref, fsthick
    namelist /control/ datadumpdir, read_g2, read_V
    namelist /electrons/ enref
    
    if(this_image() == 1) print*, 'Reading input.nml...'
    
    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Set values from input:
    ! Read allocations.
    read(1, nml = allocations)
    if(numelements < 1 .or. numatoms < 1 .or. numatoms > numelements) then
       call exit_with_message('Bad input(s) in allocations.')
    end if

    ! Set number of phonon branches.
    numbranches = 3*numatoms
        
    allocate(elements(numelements), atomtypes(numatoms), born(3,3,numatoms), &
         masses(numatoms), basis(3,numatoms), basis_cart(3,numatoms))
    
    ! Read crystal_info.
    polar = .false.
    epsilon = 0.0_dp
    born = 0.0_dp
    read(1, nml = crystal_info)
    if(any(atomtypes < 1) .or. any(masses < 0)) then
       call exit_with_message('Bad input(s) in crystal_info.')
    end if
    
    basis_cart(:,:) = matmul(lattvecs,basis)

    ! Calculate reciprocal lattice vectors and real and reciprocal cell volumes.
    do i = 1,3
       j = mod(i,3) + 1
       k = mod(j,3) + 1
       reclattvecs(:,i) = &
            cross_product(lattvecs(:,j), lattvecs(:,k))
    end do
    volume = abs(dot_product(lattvecs(:,1),reclattvecs(:,1)))
    volume_bz = twopi/volume
    reclattvecs(:,:) = volume_bz*reclattvecs(:,:)

    ! Read electrons information.
    read(1, nml = electrons)
    
    ! Read numerics information.
    qmesh = (/0, 0, 0/)
    mesh_ref = 1
    fsthick = 0.0_dp
    read(1, nml = numerics)
    if(any(qmesh <= 0) .or. mesh_ref < 1 .or. fsthick < 0) then
       call exit_with_message('Bad input(s) in numerics.')
    end if
    nq = product(qmesh)
    kmesh = qmesh*mesh_ref !Electron mesh is commensurate to and refined from the phonon mesh
    nk = product(kmesh)
    
    ! Read control information.
    datadumpdir = './'
    read_g2 = .false.
    read_V = .false.
    read(1, nml = control)
    datadumpdir = trim(datadumpdir)
    ! Create data dump directory
    if(this_image() == 1) call system('mkdir ' // trim(adjustl(datadumpdir)))

    ! Create matrix elements data directories
    g2dir = trim(adjustl(datadumpdir))//'g2'
    if(this_image() == 1) call system('mkdir ' // trim(adjustl(g2dir)))
    Vdir = trim(adjustl(datadumpdir))//'V'
    if(this_image() == 1) call system('mkdir ' // trim(adjustl(Vdir)))
    
    !Close input file.
    close(1) 

    !Set current work directory.
    call getcwd(cwd)
    cwd = trim(cwd)

    !Print out run-specific information.
    if(this_image() == 1) then
       print*, 'Working directory: ', cwd
       print*, 'Data dump directory: ', datadumpdir
       print*, 'e-ph vertex data directory: ', g2dir
       print*, 'ph-ph vertex data directory: ', Vdir
       if(read_g2) then
          print*, 'Previous e-ph vertex will be read from disk.'
       else
          print*, 'e-ph vertex will be calculated.'
       end if
       if(read_V) then
          print*, 'Previous ph-ph vertex will be read from disk.'
       else
          print*, 'ph-ph vertex will be calculated.'
       end if
    end if
    
    !Print out crystal and reciprocal lattice information.
    if(this_image() == 1) then
       print*, 'Lattice vectors [nm]:'
       print*, lattvecs(:,1)
       print*, lattvecs(:,2)
       print*, lattvecs(:,3)
       print*, 'Primitive cell volume =', volume, 'nm^3'

       print*, 'Reciprocal lattice vectors [1/nm]:'
       print*, reclattvecs(:,1)
       print*, reclattvecs(:,2)
       print*, reclattvecs(:,3)
       print*, 'Brillouin zone volume =', volume_bz, '1/nm^3'

       print*, 'Phonon wave vector mesh = ', qmesh
       print*, 'Number of phonon wave vectors = ', nq
       
       print*, 'Electron wave vector mesh = ', kmesh
       print*, 'Number of electron wave vectors = ', nk
    end if
    
    !Calculate all symmetry related information.
    call calculate_symmetries

    !Read EPW data
    call read_EPW_Wannier()
  end subroutine initialize_system

  subroutine calculate_symmetries
    !! Subroutine to generate the symmetry related data for a given crystal.
    !!
    !! This subroutine closely follows parts of config.f90 of the ShengBTE code.
    
    !Some internal variables:
    integer(k4) :: i, j, k, ii, jj, kk, ll, info
    integer(k4) :: P(3)
    integer(k4), allocatable :: rtmp(:,:,:), equiv_map(:,:)
    logical, allocatable :: valid(:)
    real(dp), allocatable :: crtmp(:,:,:), qrtmp(:,:,:)
    real(dp), allocatable :: translations(:,:), ctranslations(:,:)
    real(dp) :: tmp1(3,3), tmp2(3,3), tmp3(3,3)

    !Number of crystal symmetries.
    nsymm = get_num_operations(lattvecs,numatoms,atomtypes,basis)
    !Double the above to take time reversal symetry (TRS) into account.
    nsymm_rot = 2*nsymm

    allocate(rotations(3,3,nsymm_rot),crotations(3,3,nsymm_rot),&
         qrotations(3,3,nsymm_rot),rotations_orig(3,3,nsymm),&
         crotations_orig(3,3,nsymm),qrotations_orig(3,3,nsymm),&
         translations(3,nsymm),ctranslations(3,nsymm))

    !Get symmetry operations.
    call get_operations(lattvecs,numatoms,atomtypes,&
         basis,nsymm,rotations_orig,translations,international)
    rotations(:,:,1:nsymm) = rotations_orig

    if(this_image() == 1) then
       print*, "Symmetry group = ", trim(international)
       print*, "Number of symmetries (without time-reversal) = ", nsymm
    end if

    !Get symmertry operations in Cartesian basis.
    call get_cartesian_operations(lattvecs,nsymm,&
         rotations_orig,translations,&
         crotations_orig,ctranslations)
    crotations(:,:,1:nsymm) = crotations_orig

    !Transform the rotation matrices to the reciprocal-space basis.
    do i = 1,nsymm
       tmp1 = matmul(transpose(lattvecs),lattvecs)
       tmp2 = transpose(rotations_orig(:, :, i))
       tmp3 = tmp1
       call dgesv(3,3,tmp1,3,P,tmp2,3,info)
       qrotations_orig(:,:,i) = transpose(matmul(tmp2,tmp3))
    end do
    qrotations(:,:,1:nsymm) = qrotations_orig

    !Fill the second half of the rotation matrix list using TRS.
    rotations(:,:,nsymm+1:2*nsymm) = -rotations_orig(:,:,1:nsymm)
    qrotations(:,:,nsymm+1:2*nsymm) = -qrotations_orig(:,:,1:nsymm)
    crotations(:,:,nsymm+1:2*nsymm) = -crotations_orig(:,:,1:nsymm)

    !Find rotations that are either duplicated or incompatible with mesh.
    allocate(equiv_map(nsymm_rot,nq))
    call find_equiv_map(nsymm_rot,equiv_map,qmesh,qrotations)
    allocate(valid(nsymm_rot))
    valid = .true.
    jj = 0
    do ii = 1,nsymm_rot
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
    do ii = 2,nsymm_rot
       do i = 1,ii - 1
          if(.not. valid(i)) cycle
          if(all(rotations(:,:,ii) == rotations(:,:,i))) then
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
       allocate(rtmp(3,3,nsymm_rot - ll - jj))
       allocate(crtmp(3,3,nsymm_rot - ll - jj))
       allocate(qrtmp(3,3,nsymm_rot - ll - jj))
       kk = 0
       do ii = 1,nsymm_rot
          if(valid(ii)) then
             kk = kk + 1
             rtmp(:,:,kk) = rotations(:,:,ii)
             crtmp(:,:,kk) = crotations(:,:,ii)
             qrtmp(:,:,kk) = qrotations(:,:,ii)
          end if
       end do
       nsymm_rot = nsymm_rot - ll - jj
       call move_alloc(rtmp,rotations)
       call move_alloc(crtmp,crotations)
       call move_alloc(qrtmp,qrotations)
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
  
end module config
