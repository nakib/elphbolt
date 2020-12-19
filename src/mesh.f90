module mesh
  !! Module containing subroutines and function for setting up the
  !! full and irreducible wave vector meshes, electron and phonon quantities
  !! on these, and the various index mappings.
  
  use params, only: dp, k4
  use misc, only: demux_vector, sort_int, exit_with_message, print_message, demux_state
  use config, only: find_equiv_map
  use wannier, only: el_wann_epw, ph_wann_epw
  use data, only: el_ibz2fbz_map, el_indexlist_irred, el_nequiv, el_wavevecs, &
       el_wavevecs_irred, el_ens_irred, el_evecs_irred, el_vels_irred, kmesh, &
       nk_irred, nsymm_rot, qrotations, numwannbands, el_indexlist, el_ens, &
       el_vels, el_evecs, nk, crotations, enref, fsthick, nstates_inwindow, &
       el_fbz2ibz_map, nstates_irred_inwindow, IBZ_inwindow_states, cwd, qmesh, &
       ph_wavevecs, ph_ens, ph_vels, ph_evecs, nq, nq_irred, numbranches, &
       ph_ibz2fbz_map, ph_indexlist_irred, ph_nequiv, ph_wavevecs_irred, ph_indexlist, &
       ph_fbz2ibz_map
  
  implicit none

  private 
  public calculate_electrons, calculate_phonons
  
contains

  subroutine calc_wavevectors_full(mesh,wavevecs,blocks,indexlist)
    !! Calculate wave vectors (crystal coords.) of the full Brillouin zone (FBZ)
    !!
    !! mesh is the array of number of points along the reciprocal lattice vectors
    !! wavevecs is the list of all the wave vectors

    integer(k4), intent(in) :: mesh(3)
    logical, intent(in) :: blocks
    integer(k4), optional, intent(in) :: indexlist(:)
    real(dp), allocatable, intent(out) :: wavevecs(:,:)
    integer(k4) :: nwavevecs, ijk(3), i, imux
    
    if(blocks .and. .not. present(indexlist)) &
         call exit_with_message("If blocks is true then indexlist must be present")

    if(blocks) then
       nwavevecs = size(indexlist)
    else
       nwavevecs = product(mesh)
    end if

    allocate(wavevecs(nwavevecs, 3))
    do i = 1, nwavevecs !run over total number of vectors
       if(blocks) then
          imux = indexlist(i)
       else
          imux = i
       end if
       call demux_vector(imux, ijk, mesh, 0_k4) !get 0-based (i,j,k) indices
       wavevecs(i,:) = dble(ijk)/mesh !wave vectors in crystal coordinates
    end do
  end subroutine calc_wavevectors_full

  subroutine find_irred_wedge(mesh,nwavevecs_irred,wavevecs_irred, &
       indexlist_irred,nequivalent,nsymm_rot,qrotations,ibz2fbz_map,blocks,indexlist)
    !! Find the irreducible wedge of the FBZ and other quantities
    !! Wedge finding algorithm is inspired by ShengBTE
    !!
    !! mesh is the array of number of points along the reciprocal lattice vectors
    !! nwavevecs_irred is the number of irreducible vectors
    !! wavevecs_irred are the irreducible vectors
    !! indexlist_irred is the list containing the muxed indices of the irreducible wave vectors
    !! nequivalent lists the number of equivalent points of each irreducible point
    !! ibz2fbz_map contains the map from an irreducible BZ (IBZ) vector to its FBZ images
    !!   The third axis contains the pair (symmetry index, image)
    !! blocks chooses whether the FBZ is energy restricted or not
    !! indexlist is the sorted list of indices of the wavevectors
    !!   in the energy restricted FBZ which must be present if blocks is true

    integer(k4), intent(in) :: mesh(3)
    logical, intent(in) :: blocks
    integer(k4), intent(in) :: nsymm_rot
    real(dp), intent(in) :: qrotations(:,:,:)
    integer(k4), optional, intent(in) :: indexlist(:)
    integer(k4), intent(out) :: nwavevecs_irred
    integer(k4), allocatable, intent(out) :: indexlist_irred(:), &
         nequivalent(:), ibz2fbz_map(:,:,:)
    real(dp), allocatable, intent(out) :: wavevecs_irred(:,:)
    integer(k4) :: nwavevecs, i, imux, s, image, imagelist(nsymm_rot), &
         nrunninglist, counter, ijk(3)
    integer(k4), allocatable :: equivalence_map(:,:), runninglist(:), &
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
    else
       call find_equiv_map(nsymm_rot, equivalence_map, mesh, qrotations)
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

    if(this_image() == 1) write(*, *) "Number of FBZ wave vectors = ", counter
    if(this_image() == 1) write(*, *) "Number IBZ wave vectors = ", nwavevecs_irred

    !Deallocate some internal data
    deallocate(runninglist, equivalence_map)

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
       call demux_vector(imux, ijk, mesh, 0_k4) !get 0-based (i,j,k) indices
       
       wavevecs_irred(i,:) = dble(ijk)/mesh !wave vectors in crystal coordinates
    end do
  end subroutine find_irred_wedge

!!$  function fbz2ibz(ikmux)
!!$    !! Find index in IBZ blocks list for a given FBZ blocks muxed vector index
!!$
!!$    integer(k4), intent(in) :: ikmux
!!$    integer(k4) :: i, l, il, fbz2ibz
!!$
!!$    fbz2ibz = -1
!!$
!!$    !Sum over all ibz points
!!$    do i = 1,nk_irred !an irreducible point
!!$       do l = 1,el_nequiv(i) !number of equivalent points of i
!!$          !Get (i, l) -> il, the muxed vector index of image
!!$          il = el_ibz2fbz_map(l,i,2)
!!$
!!$          if(il == ikmux) then 
!!$             fbz2ibz = i
!!$             exit
!!$          end if
!!$          
!!$       end do
!!$    end do
!!$
!!$    if(fbz2ibz == -1) then 
!!$       print*, 'Error in fbz2ibz for input index = ', ikmux
!!$    end if
!!$  end function fbz2ibz
!!$
!!$  subroutine create_fbz2ibz_map
!!$    !! Subroutine to create map of FBZ blocks to IBZ blocks
!!$
!!$    integer(k4) :: ik
!!$
!!$    allocate(el_fbz2ibz_map(nk))
!!$
!!$    do ik = 1,nk
!!$       el_fbz2ibz_map(ik) = fbz2ibz(el_indexlist(ik))
!!$    end do
!!$  end subroutine create_fbz2ibz_map

  function fbz2ibz(iwvmux,nwv_irred,nequiv,ibz2fbz_map)
    !! Find index in IBZ blocks list for a given FBZ blocks muxed vector index

    integer(k4), intent(in) :: iwvmux, nwv_irred, nequiv(nwv_irred), ibz2fbz_map(:,:,:)
    integer(k4) :: i, l, il, fbz2ibz

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
    
    integer(k4), intent(in) :: nwv, nwv_irred, indexlist(nwv), &
         nequiv(nwv_irred),ibz2fbz_map(:,:,:)
    integer(k4), intent(out), allocatable :: fbz2ibz_map(:)

    integer(k4) :: iwv

    allocate(fbz2ibz_map(nwv))

    do iwv = 1,nwv !Run over all wave vectors
       fbz2ibz_map(iwv) = fbz2ibz(indexlist(iwv),nwv_irred,nequiv,ibz2fbz_map)
    end do
  end subroutine create_fbz2ibz_map

  subroutine apply_energy_window(nk, indexlist, energies)
    !! Subroutine to find the Fermi window restricted blocks of BZ.
    !! This could be used for FBZ and IBZ.
    !!
    !! nk is the number of mesh points - will be updated to
    !! the number of mesh points in energy restricted blocks
    !!
    !! indexlist is the list of wave vector indices - will be updated
    !! to the list of indices in the blocks

    integer(k4), intent(inout) :: nk
    integer(k4), allocatable, intent(inout) :: indexlist(:)
    real(dp), intent(in) :: energies(nk, numwannbands)

    integer(k4) :: ik, count, inwindow(nk)
    real(dp) :: aux(numwannbands)

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
    !! Subroutine to find FBZ blocks quanties by eliminating points that
    !! lie outside the energy window.

    integer(k4), intent(in) :: indexlist(:)
    real(dp), allocatable, intent(inout) :: energies(:,:), velocities(:,:,:)
    integer(k4) :: i, nk
    real(dp), allocatable :: energies_tmp(:,:), velocities_tmp(:,:,:)

    nk = size(indexlist)

    allocate(energies_tmp(nk, numwannbands), velocities_tmp(nk, numwannbands, 3))

    do i = 1, nk
       energies_tmp(i, :) = energies(indexlist(i), :)
       velocities_tmp(i, :, :) = velocities(indexlist(i), :, :)
    end do

    deallocate(energies, velocities)
    allocate(energies(nk, numwannbands), velocities(nk, numwannbands, 3))

    energies(1:nk, :) = energies_tmp(1:nk, :)
    velocities(1:nk, :, :) = velocities_tmp(1:nk, :, :)
  end subroutine fbz_blocks_quantities
  
  subroutine calculate_electrons
    !! Calculate electron energy window restricted wave vector meshes
    !! and the electronic properties on them

    !Some utitlity variables
    integer(k4) :: i, l, s, il, ib, count, istate
    real(dp), allocatable :: el_ens_tmp(:, :), el_vels_tmp(:, :, :)

    !Switch for mesh utilites with or without energy restriction
    logical :: blocks

    !I/O related
    character(len = 1024) :: filename, numcols

    !The electronic mesh setup proceeds in multiple steps:
    ! 1. Calculate full electron wave vector mesh
    call print_message("Calculating FBZ...")
    blocks = .false.
    call calc_wavevectors_full(kmesh, el_wavevecs, blocks)

    ! 2. Calculate the IBZ
    call print_message("Calculating IBZ and IBZ -> FBZ mappings...")
    call find_irred_wedge(kmesh, nk_irred, el_wavevecs_irred, &
         el_indexlist_irred, el_nequiv, nsymm_rot, qrotations, el_ibz2fbz_map, blocks)

    ! 3. Calculate IBZ quantities
    call print_message("Calculating IBZ energies...")
    allocate(el_ens_irred(nk_irred, numwannbands), &
         el_vels_irred(nk_irred, numwannbands, 3), &
         el_evecs_irred(nk_irred, numwannbands, numwannbands))
    call el_wann_epw(nk_irred, el_wavevecs_irred, el_ens_irred, &
         el_vels_irred, el_evecs_irred)

    ! 4. Map out FBZ quantities from IBZ ones
    call print_message("Mapping out FBZ energies...")
    allocate(el_indexlist(nk), el_ens(nk, numwannbands), el_vels(nk, numwannbands, 3))

    do i = 1,nk_irred !an irreducible point
       do l = 1,el_nequiv(i) !number of equivalent points of i
          il = el_ibz2fbz_map(l, i, 2) ! (i, l) -> il
          s = el_ibz2fbz_map(l, i, 1) ! mapping rotation

          !index list
          el_indexlist(il) = il

          !energy
          el_ens(il,:) = el_ens_irred(i,:)

          !velocity
          do ib = 1,numwannbands
             !here use real space (Cartesian) rotations
             el_vels(il, ib, :) = matmul(crotations(:, :, s), el_vels_irred(i, ib, :))
          end do
       end do
    end do

    ! 5. Find energy window restricted FBZ blocks
    !    After this step, nk, el_indexlist will refer
    !    to the energy restricted mesh
    call print_message("Calculating Fermi window restricted FBZ blocks...")
    call apply_energy_window(nk, el_indexlist, el_ens)

    ! 6. Sort index list and related quanties of FBZ blocks
    call print_message("Sorting FBZ blocks index list...")
    call sort_int(el_indexlist)

    ! 7. Get FBZ blocks wave vectors, energies, velocities and eigenvectors
    !    After this step, el%wavevecs, el%ens, el%vels, and el%evecs
    !    will refer to the energy restricted mesh
    call print_message("Calcutating FBZ blocks quantities...")
    !wave vectors
    deallocate(el_wavevecs)
    blocks = .true.
    call calc_wavevectors_full(kmesh, el_wavevecs, blocks, el_indexlist) !wave vectors

    !energies and velocities
    call fbz_blocks_quantities(el_indexlist, el_ens, el_vels)

    !Get FBZ blocks eigenvectors from direct calculations since we are
    !not getting these from IBZ quantities via symmetry rotations
    allocate(el_evecs(nk, numwannbands, numwannbands))
    allocate(el_ens_tmp(nk, numwannbands), el_vels_tmp(nk, numwannbands, 3))
    call el_wann_epw(nk, el_wavevecs, el_ens_tmp, el_vels_tmp, el_evecs)
    deallocate(el_ens_tmp, el_vels_tmp) !free up memory

    ! 8. Find IBZ of energy window restricted blocks
    !    After this step, el%nk_irred, el%indexlist_irred, 
    !    el%wavevecs_irred, el%nequiv, and el%ibz2fbz_map
    !    will refer to the energy restricted mesh 
    call print_message("Calculating IBZ blocks...")
    deallocate(el_wavevecs_irred, el_indexlist_irred, el_nequiv, el_ibz2fbz_map)
    blocks = .true.
    call find_irred_wedge(kmesh, nk_irred, el_wavevecs_irred, &
         el_indexlist_irred, el_nequiv, nsymm_rot, qrotations, &
         el_ibz2fbz_map, blocks, el_indexlist)
    
    ! 9. Get IBZ blocks energies, velocities, and eigen vectors
    !    energies and velocities
    call print_message("Calcutating IBZ blocks quantities...")
    deallocate(el_ens_irred, el_vels_irred, el_evecs_irred)
    allocate(el_ens_irred(nk_irred, numwannbands), &
         el_vels_irred(nk_irred, numwannbands, 3), &
         el_evecs_irred(nk_irred, numwannbands, numwannbands))
    call el_wann_epw(nk_irred, el_wavevecs_irred, el_ens_irred, &
         el_vels_irred, el_evecs_irred)

    ! 10. Calculate the number of FBZ blocks electronic states
    !     available for scattering
    nstates_inwindow = 0
    do i = 1,nk !over FBZ blocks
       do ib = 1,numwannbands !bands
          if(abs(el_ens(i, ib) - enref) <= fsthick) &
               nstates_inwindow = nstates_inwindow + 1
       end do
    end do
    if(this_image() == 1) write(*, *) "Number of energy restricted FBZ blocks states = ", nstates_inwindow

    ! 11. Create FBZ blocks to IBZ blocks map
    call print_message("Calculating FBZ -> IBZ mappings...")
    !call create_fbz2ibz_map
    call create_fbz2ibz_map(el_fbz2ibz_map,nk,nk_irred,el_indexlist,el_nequiv,el_ibz2fbz_map)

    ! 12. Calculate the number of IBZ electronic states available for scattering
    nstates_irred_inwindow = 0
    do istate = 1,nk_irred*numwannbands
       !Demux state index into band (ib) and wave vector (i) indices
       call demux_state(istate, numwannbands, ib, i)
       if(abs(el_ens_irred(i, ib) - enref) <= fsthick) then 
          nstates_irred_inwindow = nstates_irred_inwindow + 1
       end if
    end do
    if(this_image() == 1) write(*, *) "Number of energy restricted IBZ blocks states = ", &
         nstates_irred_inwindow

    !Calculate list of IBZ in-window states = (wave vector index, band index)
    allocate(IBZ_inwindow_states(nstates_irred_inwindow,2))
    count = 0
    do istate = 1, nk_irred*numwannbands
       !Demux state index into band (ib) and wave vector (i) indices
       call demux_state(istate, numwannbands, ib, i)
       if(abs(el_ens_irred(i, ib) - enref) <= fsthick) then
          count = count + 1
          IBZ_inwindow_states(count,:) = (/i, ib/)
       end if
    end do

    !Write IBZ in-window states as text data to file
    if(this_image() == 1) then
       call chdir(cwd)
       filename = 'IBZ_inwindow_states'
       write(numcols,"(I0)") 2
       open(1,file=trim(filename),status='replace')
       write(1,*) "#k-vec index     band index"
       do i = 1, nstates_irred_inwindow
          write(1,"("//trim(adjustl(numcols))//"I10)") IBZ_inwindow_states(i,:)
       end do
       close(1)
    end if
  end subroutine calculate_electrons

!!$  subroutine calculate_phonons
!!$    !! Calculate phonon quantities on the FBZ and IBZ meshes.
!!$
!!$    integer(k4) :: iq
!!$    !Switch for mesh utilites with or without energy restriction
!!$    logical :: blocks
!!$    blocks = .false.
!!$
!!$    call print_message("Calculating phonon FBZ quantities...")
!!$
!!$    allocate(ph_indexlist(nq))
!!$    do iq = 1, nq
!!$       ph_indexlist(iq) = iq
!!$    end do
!!$    
!!$    !Calculate FBZ mesh
!!$    call calc_wavevectors_full(qmesh, ph_wavevecs, blocks)
!!$
!!$    !Calculate FBZ phonon quantities
!!$    allocate(ph_ens(nq, numbranches))
!!$    allocate(ph_vels(nq, numbranches, 3))
!!$    allocate(ph_evecs(nq, numbranches, numbranches))    
!!$    call ph_wann_epw(nq, ph_wavevecs, ph_ens, ph_vels, ph_evecs)
!!$
!!$    !Calculate IBZ mesh
!!$    call print_message("Calculating IBZ and IBZ -> FBZ mappings...")
!!$    call find_irred_wedge(qmesh, nq_irred, ph_wavevecs_irred, &
!!$         ph_indexlist_irred, ph_nequiv, nsymm_rot, qrotations, ph_ibz2fbz_map, blocks)
!!$
!!$    !Create FBZ to IBZ map
!!$    call print_message("Calculating FBZ -> IBZ mappings...")
!!$    call create_fbz2ibz_map(ph_fbz2ibz_map, nq, nq_irred, ph_indexlist, ph_nequiv, ph_ibz2fbz_map)
!!$  end subroutine calculate_phonons
end module mesh
