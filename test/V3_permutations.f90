program V3_permutations
  !! generation and classification of irreducible triplets using permutation symmetry
  !!
  !! Returns a mapping of full triplets to a reduced consisting of (ilambda2, ilambda1) pairs
  !! lambda is defined as (iband, i * k), where iband is the band index and ik wave-vector index

  use precision, only: i64, r64
  use misc, only: print_message, subtitle, timer, exit_with_message, mux_vector, demux_state, &
       mux_state, twonorm, permutations, demux_vector, lex_less_1d
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use phonon_module, only: phonon

  implicit none

  abstract interface
     real(r64) function Vm2_3ph(ev1_s1, ev2_s2, ev3_s3, &
          Index_i, Index_j, Index_k, ifc3, phases_q2q3, ntrip, nb)
       import r64, i64

       integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
       complex(r64), intent(in) :: phases_q2q3(ntrip), ev1_s1(nb), ev2_s2(nb), ev3_s3(nb)
       real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)
     end function Vm2_3ph
  end interface

  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(phonon) :: ph
  type(timer) :: t_event
  procedure(Vm2_3ph), pointer :: Vm2_calculator
  complex(r64), allocatable :: R1(:, :, :), R3(:)
  complex(r64) :: ev1(3, 2), ev2(3, 2), ev3(3, 2)

  real(r64), allocatable :: V2(:, :, :, :, :)
  real(r64), allocatable :: V2_minimal_set(:, :, :, :, :)
  integer(i64), allocatable :: M(:, :, :)

  !1st axis: canonical triplet (istate1, istate2, istate3)
  !2nd axis: istate2
  !3rd axis: istate1
  !integer(i64), allocatable :: triplet_permutation_maps(:, :, :)

  integer :: count_full, count_minimal
  integer(i64) :: lambda1, lambda2, lambda3
  integer(i64) :: istate1, iq2, iq3_minus, s2, s3, nstates_irred, s1, iq1_ibz, &
       iq1
  real(r64) :: val !value

  !call triplet_test()

  if(this_image() == 1) then
     write(*, '(A)')  'V3 permutations playground'
     write(*, '(A, I5)') 'Number of coarray images = ', num_images()
  end if

  !Set up crystal
  call crys%initialize

  !Set up numerics data
  call num%initialize(crys)

  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)

  !Calculate phonons
  call ph%initialize(crys, sym, num)

  !Calculate ph-ph vertex (cpu, original)
  Vm2_calculator => Vm2_3ph_reference
  call t_event%start_timer('Full V set calculation')
  call calculate_3ph_interaction(ph, crys, num, V2, Vm2_calculator)
  print*, 'value = ', twonorm(pack(V2, .true.))
  call t_event%end_timer('Full V set calculation')
  !print*, 'value = ', twonorm(pack(V2_minimal_set, .true.))
  ! print*, V2(3, 1, 2, 3, 3)
  ! print*, V2(3, 2, 4, 3, 2)
  ! print*, V2(4, 1, 4, 1, 4)
  ! print*, V2(5, 1, 5, 1, 5)
  ! print*, V2(6, 1, 6, 1, 6)

  !Minimal set calculations
  call t_event%start_timer('Minimal V set calculation')
!!$  call map_triplet_full_to_reduced_new(ph%numbands, ph%wvmesh, &
!!$       ph%nwv_irred*ph%numbands, ph%nwv*ph%numbands, &
!!$       triplet_permutation_maps)
  call calculate_3ph_interaction_minimalset_new(ph, crys, num, &
       V2, Vm2_calculator)
  print*, 'value = ', twonorm(pack(V2, .true.))
  call t_event%end_timer('Minimal V set calculation')

!!$  count_full = 0
!!$  print *, '---------------------------------------------------------------'
!!$  print *, '   lambda1    lambda2    lambda3      Value'
!!$  print *, '---------------------------------------------------------------'
!!$
!!$  do lambda1 = 1, 6
!!$     call demux_state(lambda1, ph%numbands, s1, iq1)
!!$
!!$     do lambda2 = 1, 10
!!$        call demux_state(lambda2, ph%numbands, s2, iq2)
!!$
!!$        do lambda3 = 1, 10
!!$           call demux_state(lambda3, ph%numbands, s3, iq3_minus)
!!$
!!$           if(s3 <= size(V2, 1) .and. iq3_minus <= size(V2, 2) .and. &
!!$                s2 <= size(V2, 3) .and. iq2 <= size(V2, 4) .and. &
!!$                lambda1 <= size(V2, 5)) then
!!$
!!$              val = V2(s3, iq3_minus, s2, iq2, lambda1)
!!$
!!$              !if(abs(val) > 1.0e-7_r64) then
!!$              write(*,'(3I10, 2X, F16.12)') lambda1, lambda2, lambda3, val
!!$              count_full = count_full + 1
!!$              !end if
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  print *, 'Number of V2 elements:', count_full
!!$
!!$  !Calculate V2_minimal_set (reference)
!!$  Vm2_calculator => Vm2_3ph_reference
!!$  call t_event%start_timer('reference V2 minimal set')
!!$  call calculate_3ph_interaction_minimalset(ph, crys, num, V2_minimal_set, M, Vm2_calculator)
!!$  call t_event%end_timer('reference V2 minimal set')
!!$  print*, 'value = ', twonorm(pack(V2_minimal_set, .true.))
!!$  !print*, 'value = ', twonorm(pack(V2_minimal_set, .true.))
!!$  ! print*, V2(3, 1, 2, 3, 3)
!!$  ! print*, V2(3, 2, 4, 3, 2)
!!$  ! print*, V2(4, 1, 4, 1, 4)
!!$  ! print*, V2(5, 1, 5, 1, 5)
!!$  ! print*, V2(6, 1, 6, 1, 6)
!!$
!!$  count_minimal = 0
!!$  print *, '---------------------------------------------------------------'
!!$  print *, '   lambda1    lambda2    lambda3      Value'
!!$  print *, '---------------------------------------------------------------'
!!$
!!$  do lambda1 = 1, 6
!!$     call demux_state(lambda1, ph%numbands, s1, iq1)
!!$
!!$     do lambda2 = 1, 10
!!$        call demux_state(lambda2, ph%numbands, s2, iq2)
!!$
!!$        do lambda3 = 1, 10
!!$           call demux_state(lambda3, ph%numbands, s3, iq3_minus)
!!$
!!$           !canonical value from V2_minimal_set
!!$           if(s3 <= size(V2_minimal_set, 1) .and. iq3_minus <= size(V2_minimal_set, 2) .and. &
!!$                s2 <= size(V2_minimal_set, 3) .and. iq2 <= size(V2_minimal_set, 4) .and. &
!!$                lambda1 <= size(V2_minimal_set, 5)) then
!!$
!!$              val = V2_minimal_set(s3, iq3_minus, s2, iq2, lambda1)
!!$
!!$              if(abs(val) /= 1.0_r64) then
!!$                 write(*,'(3I10, 2X, F16.12)') lambda1, lambda2, lambda3, val
!!$                 count_minimal = count_minimal + 1
!!$              end if
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  print *, 'Number of V2 minimal set elements:', count_minimal

!!$  ! Reduction factor between V2 and V2 minimal set
!!$  print *, 'Reduction factor (minimal/full):', real(count_minimal)/real(count_full)
!!$  print *, 'Symmetry saving (in %):', (1.0 - real(count_minimal)/real(count_full))*100.0

  !call triplet_test

contains

  subroutine triplet_test()
    integer(i64) :: nbands, mesh_size(3)
    integer(i64), allocatable :: lambda1_list(:), lambda2_list(:)
    integer(i64), allocatable :: M(:, :, :)
    integer :: ilambda1, ilambda2

    if (this_image() == 1) then
       print *, "V3 permutations to canonical (lambda1, lambda2, lambda3) form"
       print *, "Number of coarray images =", num_images()
    end if

    nbands = 2
    mesh_size = [2, 2, 2]*1_i64

    allocate(lambda1_list(2)); lambda1_list = [1, 2]*1_i64
    allocate(lambda2_list(16)); lambda2_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]*1_i64

    call map_triplet_full_to_reduced(nbands, mesh_size, lambda1_list, lambda2_list, M)

    print *, "M(ilambda2, ilambda1) → canonical triplet: (lambda1, lambda2, lambda3)"
    do ilambda1 = 1, size(lambda1_list)
       do ilambda2 = 1, size(lambda2_list)
          write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'M(', ilambda2, ',', ilambda1, ') = (', &
               M(1, ilambda2, ilambda1), ',', M(2, ilambda2, ilambda1), ',', M(3, ilambda2, ilambda1), ')'
       end do
    end do
  end subroutine triplet_test

  subroutine map_triplet_full_to_reduced(nbands, mesh_size, lambda1_list, lambda2_list, M)
    !! This subroutine maps the canonical triplet back to (ilambda2, ilambda1) and stores it in M(:, ilambda2, ilambda1).
    !!
    !! nbands Number of bands
    !! mesh_size Wavevector discretization
    !! lambda1_list The smaller list of states
    !! lambda2_list The larger list of states
    !! M Mapping of every interaction triplet of states to its canonical representative

    integer(i64), intent(in) :: nbands, mesh_size(3)
    integer(i64), intent(in) :: lambda1_list(:), lambda2_list(:)
    integer(i64), allocatable, intent(out) :: M(:, :, :)

    integer(i64) :: ilambda1, ilambda2, iband1, iband2, iband3, ik1, ik2, ik3
    integer(i64) :: q1(3), q2(3), q3(3)
    integer(i64), allocatable :: all_perms(:, :)
    integer(i64) :: triplet_full(3), permuted_triplet(3), canonical_triplet(3)
    integer(i64) :: i, j

    allocate(M(3, size(lambda2_list), size(lambda1_list)))

    ! Initialize all element of M to a known integer value -1: mean that the triplet has not been assigned yet
    M = -1_i64

    ! Get all permutations of 3 indices
    all_perms = permutations(3_i64)
    do ilambda1 = 1, size(lambda1_list)
       ! Demux state for lambda1_list: for each lambda1_list value, convert its index to (iband1, ik1)
       call demux_state(lambda1_list(ilambda1), nbands, iband1, ik1)

       ! Demux vector to get q1
       call demux_vector(ik1, q1, mesh_size, base = 0_i64)
       do ilambda2 = 1, size(lambda2_list)
          ! Demux state for lambda2_list: for each lambda2_list value, convert its index to (iband2, ik2)
          call demux_state(lambda2_list(ilambda2), nbands, iband2, ik2)

          ! Demux vector to get q1
          call demux_vector(ik2, q2, mesh_size, base = 0_i64)

          ! Compute q3 using modular arithmetic: such that momentum is conserved: q1 - q2 - q3 = 0 mod G
          q3 = modulo(q1 - q2, mesh_size)

          ! Compute ik3 using mux_vector
          ik3 = mux_vector(q3, mesh_size, base = 0_i64)

          ! Iterate over all possible bands for third state: to have (iband3, ik3)
          do iband3 = 1, nbands
             triplet_full = [mux_state(nbands, iband1, ik1), &
                  mux_state(nbands, iband2, ik2), &
                  mux_state(nbands, iband3, ik3)]

             ! Initialize the canonical form of the triplet to the original (unpermuted) triplet, before trying other triplet
             ! triplet_full is a 3x2 array: ((iband1, ik1), (iband2, ik2), (iband3, ik3))
             canonical_triplet = triplet_full

             do i = 1, size(all_perms, 2)
                do j = 1, 3
                   permuted_triplet(j) = triplet_full(all_perms(j, i))
                end do

                ! Use a lexicographic comparison function to keep the smallest permutation
                if(lex_less_1d(permuted_triplet, canonical_triplet)) canonical_triplet = permuted_triplet
             end do

             ! after elaborating the triplet (lambda1, lambda2, lambda3) as linear state indices
             ! and after generating all 6 permutations
             ! keep only the lexicographically smallest permutation (the irreducible triplet)
             ! Store it in M(:, ilambda2, ilambda1)

             M(:, ilambda2, ilambda1) = canonical_triplet

             ! once the first valid canonical triplet is found (for a given lambda1 and lambda2, we can store it and skip checking (with exit) other values of iband3
             ! symmetry will handle all others, we only interested in one representative triplet identified using lex_less_1d.
             exit
          end do
       end do
    end do
  end subroutine map_triplet_full_to_reduced

  subroutine map_triplet_full_to_reduced_new(nbands, mesh_size, &
       nstates_irred, nstates, M)
    !! This subroutine maps the canonical triplet back to (ilambda2, ilambda1) and stores it in M(:, ilambda2, ilambda1).
    !!
    !! nbands Number of bands
    !! mesh_size Wavevector discretization
    !! nstates_irred Number of irreducible states
    !! nstates Number of states
    !! M Mapping of every interaction triplet of states to its canonical representative

    integer(i64), intent(in) :: nbands, mesh_size(3), nstates_irred, nstates
    integer(i64), allocatable, intent(out) :: M(:, :, :)

    integer(i64) :: ilambda1, ilambda2, iband1, iband2, iband3, ik1, ik2, ik3
    integer(i64) :: q1(3), q2(3), q3(3)
    integer(i64), allocatable :: all_perms(:, :)
    integer(i64) :: triplet_full(3), permuted_triplet(3), canonical_triplet(3)
    integer(i64) :: i, j

    allocate(M(3, nstates, nstates_irred))

    ! Initialize all element of M to a known integer value -1: mean that the triplet has not been assigned yet
    M = -1_i64

    ! Get all permutations of 3 indices
    all_perms = permutations(3_i64)
    do ilambda1 = 1, nstates_irred
       ! Demux state for lambda1_list: for each lambda1_list value, convert its index to (iband1, ik1)
       call demux_state(ilambda1, nbands, iband1, ik1)

       ! Demux vector to get q1
       call demux_vector(ik1, q1, mesh_size, base = 0_i64)

       do ilambda2 = 1, nstates
          ! Demux state for lambda2_list: for each lambda2_list value, convert its index to (iband2, ik2)
          call demux_state(ilambda2, nbands, iband2, ik2)

          ! Demux vector to get q1
          call demux_vector(ik2, q2, mesh_size, base = 0_i64)

          ! Compute q3 using modular arithmetic: such that momentum is conserved: q1 - q2 - q3 = 0 mod G
          q3 = modulo(q1 - q2, mesh_size)

          ! Compute ik3 using mux_vector
          ik3 = mux_vector(q3, mesh_size, base = 0_i64)

          ! Iterate over all possible bands for third state: to have (iband3, ik3)
          do iband3 = 1, nbands
             triplet_full = [mux_state(nbands, iband1, ik1), &
                  mux_state(nbands, iband2, ik2), &
                  mux_state(nbands, iband3, ik3)]

             ! Initialize the canonical form of the triplet to the original (unpermuted) triplet, before trying other triplet
             ! triplet_full is a 3x2 array: ((iband1, ik1), (iband2, ik2), (iband3, ik3))
             canonical_triplet = triplet_full

             do i = 1, size(all_perms, 2)
                do j = 1, 3
                   permuted_triplet(j) = triplet_full(all_perms(j, i))
                end do

                ! Use a lexicographic comparison function to keep the smallest permutation
                if(lex_less_1d(permuted_triplet, canonical_triplet)) canonical_triplet = permuted_triplet
             end do

             ! after elaborating the triplet (lambda1, lambda2, lambda3) as linear state indices
             ! and after generating all 6 permutations
             ! keep only the lexicographically smallest permutation (the irreducible triplet)
             ! Store it in M(:, ilambda2, ilambda1)

             M(:, ilambda2, ilambda1) = canonical_triplet

             ! once the first valid canonical triplet is found (for a given lambda1 and lambda2, we can store it and skip checking (with exit) other values of iband3
             ! symmetry will handle all others, we only interested in one representative triplet identified using lex_less_1d.
             exit
          end do
       end do
    end do
  end subroutine map_triplet_full_to_reduced_new

  subroutine map_triplet_full_to_reduced_ilambda1(nbands, mesh_size, ilambda1, nstates, M)
    !! This subroutine maps the canonical triplet back to (ilambda2, ilambda1) and stores it in M_ilambda1(:, ilambda2).
    !!
    !! nbands Number of bands
    !! mesh_size Wavevector discretization
    !! ilambda1 Initial phonon state
    !! nstates Number of states
    !! M Mapping of every interaction triplet of states to its canonical representative

    integer(i64), intent(in) :: nbands, mesh_size(3), ilambda1, nstates
    integer(i64), allocatable, intent(out) :: M(:, :)

    integer(i64) :: ilambda2, iband1, iband2, iband3, ik1, ik2, ik3
    integer(i64) :: q1(3), q2(3), q3(3)
    integer(i64), allocatable :: all_perms(:, :)
    integer(i64) :: triplet_full(3), permuted_triplet(3), canonical_triplet(3)
    integer(i64) :: i, j

    allocate(M(3, nstates))

    ! Get all permutations of 3 indices
    all_perms = permutations(3_i64)

    ! Demux the initial phonon state
    call demux_state(ilambda1, nbands, iband1, ik1)

    ! Demux vector to get q1
    call demux_vector(ik1, q1, mesh_size, base = 0_i64)

    do ilambda2 = 1, nstates
       ! Demux state for lambda2_list: for each lambda2_list value, convert its index to (iband2, ik2)
       call demux_state(ilambda2, nbands, iband2, ik2)

       ! Demux vector to get q1
       call demux_vector(ik2, q2, mesh_size, base = 0_i64)

       ! Compute q3 using modular arithmetic: such that momentum is conserved: q1 - q2 - q3 = 0 mod G
       q3 = modulo(q1 - q2, mesh_size)

       ! Compute ik3 using mux_vector
       ik3 = mux_vector(q3, mesh_size, base = 0_i64)

       ! Iterate over all possible bands for third state: to have (iband3, ik3)
       do iband3 = 1, nbands
          triplet_full = [mux_state(nbands, iband1, ik1), &
               mux_state(nbands, iband2, ik2), &
               mux_state(nbands, iband3, ik3)]

          ! Initialize the canonical form of the triplet to the original (unpermuted) triplet, before trying other triplet
          ! triplet_full is a 3x2 array: ((iband1, ik1), (iband2, ik2), (iband3, ik3))
          canonical_triplet = triplet_full

          do i = 1, size(all_perms, 2)
             do j = 1, 3
                permuted_triplet(j) = triplet_full(all_perms(j, i))
             end do

             ! Use a lexicographic comparison function to keep the smallest permutation
             if(lex_less_1d(permuted_triplet, canonical_triplet)) canonical_triplet = permuted_triplet
          end do

          ! after elaborating the triplet (lambda1, lambda2, lambda3) as linear state indices
          ! and after generating all 6 permutations
          ! keep only the lexicographically smallest permutation (the irreducible triplet)
          ! Store it in M(:, ilambda2, ilambda1)

          M(:, ilambda2) = canonical_triplet

          ! once the first valid canonical triplet is found (for a given lambda1 and lambda2, we can store it and skip checking (with exit) other values of iband3
          ! symmetry will handle all others, we only interested in one representative triplet identified using lex_less_1d.
          exit
       end do
    end do
  end subroutine map_triplet_full_to_reduced_ilambda1

  subroutine calculate_3ph_interaction(ph, crys, num, V2, Vm2_calculator)
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    real(r64), allocatable, intent(out) :: V2(:, :, :, :, :)
    procedure(Vm2_3ph), pointer, intent(in) :: Vm2_calculator

    !Local variables
    integer(i64) :: istate1, nstates_irred, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
         idim, jdim, s2s3, counter
    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         aux
    complex(r64) :: phases(ph%numtriplets)

    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

!!$    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))
!!$
!!$    V2 = 0.0

    counter = 0

    !Run over first phonon IBZ states
    !do istate1 = 1, nstates_irred
    !Demux state index into branch (s) and wave vector (iq) indices
    !call demux_state(istate1, ph%numbands, s1, iq1_ibz)

    do iq1_ibz = 1, ph%nwv_irred

       !Muxed index of wave vector from the IBZ index list.
       !This will be used to access IBZ information from the FBZ quantities.
       iq1 = ph%indexlist_irred(iq1_ibz)

       !Initial (IBZ blocks) wave vector (crystal coords.)
       q1 = ph%wavevecs(iq1, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%wvmesh)

       do iq2 = 1, ph%nwv
          !Initial (IBZ blocks) wave vector (crystal coords.)
          q2 = ph%wavevecs(iq2, :)

          !Convert from crystal to 0-based index vector
          q2_indvec = nint(q2*ph%wvmesh)

          !Folded final phonon wave vector
          q3_minus_indvec = modulo(q1_indvec - q2_indvec, ph%wvmesh) !0-based index vector
          q3_minus = q3_minus_indvec/dble(ph%wvmesh) !crystal coords.

          !Muxed index of q3_minus
          iq3_minus = mux_vector(q3_minus_indvec, ph%wvmesh, 0_i64)

          q2_cart = matmul(crys%reclattvecs, q2)
          q3_minus_cart = matmul(crys%reclattvecs, q3_minus)

!!$          phases = exp((0.0_r64, -1.0_r64)* &
!!$               (matmul(q2_cart, ph%R_j) + matmul(q3_minus_cart, ph%R_k)))

          do s1 = 1, ph%numbands

             istate1 = mux_state(ph%numbands, s1, iq1_ibz)

             !Combined loop over the 2nd and 3rd phonon bands
             do s2s3 = 1, ph%numbands**2
                s2 = int((s2s3 - 1)/ph%numbands) + 1 !changes slow
                s3 = modulo(s2s3 - 1, ph%numbands) + 1 !changes fast

                counter = counter + 1

!!$                aux = Vm2_calculator(ph%evecs(iq1, s1, :), &
!!$                     ph%evecs(iq2, s2, :), ph%evecs(iq3_minus, s3, :), &
!!$                     ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:, :, :, :), &
!!$                     phases(:), ph%numtriplets, ph%numbands)
!!$
!!$                V2(s3, iq3_minus, s2, iq2, istate1) = aux
             end do  !s2s3
          end do     !s1
       end do        !iq2
    end do           !iq1_ibz

    print*, 'Number of matrix element computed = ', counter
  end subroutine calculate_3ph_interaction

!!$  subroutine calculate_3ph_interaction_minimalset_new(ph, crys, num, &
!!$       permutations_map, V2, Vm2_calculator)
  subroutine calculate_3ph_interaction_minimalset_new(ph, crys, num, &
       V2, Vm2_calculator)
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    !integer(i64), intent(in) :: permutations_map(:, :, :)    
    real(r64), allocatable, intent(out) :: V2(:, :, :, :, :)
    procedure(Vm2_3ph), pointer, intent(in) :: Vm2_calculator

    !Local variables
    integer(i64) :: istate1, istate2, istate3, nstates_irred, nstates, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
         idim, jdim, s2s3, counter
    integer(i64), allocatable :: permutations_map_istate1(:, :)
    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         aux
    complex(r64) :: phases(ph%numtriplets)

    !Total number of IBZ and FBZ states
    nstates_irred = ph%nwv_irred*ph%numbands
    nstates = ph%nwv*ph%numbands

    !TEST
    !print*, 'dimensions of permutations_map = ', size(permutations_map)
    !print*, 'nstates_irred, nstates = ', nstates, nstates_irred
    !!

!!$    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))
!!$
!!$    V2 = 0.0

    counter = 0

    !Run over first phonon IBZ states
    do istate1 = 1, nstates_irred
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate1, ph%numbands, s1, iq1_ibz)

       !do iq1_ibz = 1, ph%nwv_irred

       !Muxed index of wave vector from the IBZ index list.
       !This will be used to access IBZ information from the FBZ quantities.
       iq1 = ph%indexlist_irred(iq1_ibz)

       !Initial (IBZ blocks) wave vector (crystal coords.)
       q1 = ph%wavevecs(iq1, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%wvmesh)

       call map_triplet_full_to_reduced_ilambda1(ph%numbands, &
            ph%wvmesh, istate1, nstates, permutations_map_istate1)

       do iq2 = 1, ph%nwv
          !Initial (IBZ blocks) wave vector (crystal coords.)
          q2 = ph%wavevecs(iq2, :)

          !Convert from crystal to 0-based index vector
          q2_indvec = nint(q2*ph%wvmesh)

          !Folded final phonon wave vector
          q3_minus_indvec = modulo(q1_indvec - q2_indvec, ph%wvmesh) !0-based index vector
          q3_minus = q3_minus_indvec/dble(ph%wvmesh) !crystal coords.

          !Muxed index of q3_minus
          iq3_minus = mux_vector(q3_minus_indvec, ph%wvmesh, 0_i64)

          q2_cart = matmul(crys%reclattvecs, q2)
          q3_minus_cart = matmul(crys%reclattvecs, q3_minus)

!!$          phases = exp((0.0_r64, -1.0_r64)* &
!!$               (matmul(q2_cart, ph%R_j) + matmul(q3_minus_cart, ph%R_k)))

          !do s1 = 1, ph%numbands
          !istate1 = mux_state(ph%numbands, s1, iq1_ibz)

          do s2 = 1, ph%numbands
             istate2 = mux_state(ph%numbands, s2, iq2)

             do s3 = 1, ph%numbands
                istate3 = mux_state(ph%numbands, s3, iq3_minus)

                !Need only compute the matrix element for one of the permutations
                if(all([istate1, istate2, istate3] == &
                     permutations_map_istate1(:, istate2))) then
                   !Count how many processes were explicitly computed
                   counter = counter + 1

!!$                      aux = Vm2_calculator(ph%evecs(iq1, s1, :), &
!!$                           ph%evecs(iq2, s2, :), ph%evecs(iq3_minus, s3, :), &
!!$                           ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:, :, :, :), &
!!$                           phases(:), ph%numtriplets, ph%numbands)
!!$
!!$                      !Beware: Here only the minimal subset will be non-zero
!!$                      V2(s3, iq3_minus, s2, iq2, istate1) = aux
                end if
             end do
          end do
          !end do
       end do
    end do

    print*, 'Number of matrix element computed = ', counter
  end subroutine calculate_3ph_interaction_minimalset_new

  subroutine calculate_3ph_interaction_minimalset(ph, crys, num, V2_minimal_set, M, Vm2_calculator)
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    real(r64), allocatable, intent(out) :: V2_minimal_set(:, :, :, :, :)
    integer(i64), allocatable, intent(out) :: M(:, :, :) ! M(3, istate2, istate1)
    procedure(Vm2_3ph), pointer, intent(in) :: Vm2_calculator

    !Local variables
    integer(i64) :: istate1, istate2, nstates_irred, nstates_full, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
         idim, jdim, s2s3, lambda1, lambda2, lambda3
    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         aux
    integer(i64), allocatable :: all_perms(:, :)
    integer(i64) :: triplet_full(3), permuted_triplet(3), canonical_triplet(3)
    integer(i64) :: i, j
    complex(r64) :: phases(ph%numtriplets)

    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands
    nstates_full  = ph%nwv * ph%numbands

    ! Get all permutations of 3 indices
    all_perms = permutations(3_i64)

    allocate(V2_minimal_set(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))
    allocate(M(3, nstates_full, nstates_irred))

    V2_minimal_set = 0.0

    !M(3, istate2, istate1) stores the canonical triplet for each pair of states (lambda1, lambda2)
    !istate1 runs over irreducible states: nstates_irred
    !istate2 runs over all full states

    !Run over first phonon IBZ states
    !do istate1 = 1, nstates_irred
    !Demux state index into branch (s) and wave vector (iq) indices
    !call demux_state(istate1, ph%numbands, s1, iq1_ibz)

    do iq1_ibz = 1, ph%nwv_irred

       !Muxed index of wave vector from the IBZ index list.
       !This will be used to access IBZ information from the FBZ quantities.
       iq1 = ph%indexlist_irred(iq1_ibz)

       !Initial (IBZ blocks) wave vector (crystal coords.)
       q1 = ph%wavevecs(iq1, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%wvmesh)

       do iq2 = 1, ph%nwv
          !Initial (IBZ blocks) wave vector (crystal coords.)
          q2 = ph%wavevecs(iq2, :)

          !Convert from crystal to 0-based index vector
          q2_indvec = nint(q2*ph%wvmesh)

          !Folded final phonon wave vector
          q3_minus_indvec = modulo(q1_indvec - q2_indvec, ph%wvmesh) !0-based index vector
          q3_minus = q3_minus_indvec/dble(ph%wvmesh) !crystal coords.

          !Muxed index of q3_minus
          iq3_minus = mux_vector(q3_minus_indvec, ph%wvmesh, 0_i64)

          q2_cart = matmul(crys%reclattvecs, q2)
          q3_minus_cart = matmul(crys%reclattvecs, q3_minus)

          phases = exp((0.0_r64, -1.0_r64)* &
               (matmul(q2_cart, ph%R_j) + matmul(q3_minus_cart, ph%R_k)))

          do s1 = 1, ph%numbands

             istate1 = mux_state(ph%numbands, s1, iq1_ibz)

             !Combined loop over the 2nd and 3rd phonon bands
             !do s2s3 = 1, ph%numbands**2
             !s2 = int((s2s3 - 1)/ph%numbands) + 1 !changes slow
             !s3 = modulo(s2s3 - 1, ph%numbands) + 1 !changes fast
             do s2 = 1, ph%numbands
                istate2 = mux_state(ph%numbands, s2, iq2) !istate2 spans the FBZ, that makes it the long list
                do s3 = 1, ph%numbands

                   !Convert (s, iq) indices to canonical triplet (lambda1, lambda2, lambda3)
                   lambda1 = mux_state(ph%numbands, s1, iq1)
                   lambda2 = mux_state(ph%numbands, s2, iq2)
                   lambda3 = mux_state(ph%numbands, s3, iq3_minus)
                   triplet_full = [lambda1, lambda2, lambda3]

                   canonical_triplet = triplet_full

                   do i = 1, size(all_perms, 2)
                      do j = 1, 3
                         permuted_triplet(j) = triplet_full(all_perms(j, i))
                      end do

                      ! Use a lexicographic comparison function to keep the smallest permutation
                      if(lex_less_1d(permuted_triplet, canonical_triplet)) canonical_triplet = permuted_triplet
                   end do

                   ! Store the first valid canonical triplet found in M(:, istate2, istate1)
                   !M(:, istate2, istate1) = canonical_triplet

                   !only compute and store the value for the canonical triplet
                   if(all(triplet_full == canonical_triplet)) then
                      aux = Vm2_calculator(ph%evecs(iq1, s1, :), &
                           ph%evecs(iq2, s2, :), ph%evecs(iq3_minus, s3, :), &
                           ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:, :, :, :), &
                           phases(:), ph%numtriplets, ph%numbands)

                      V2_minimal_set(s3, iq3_minus, s2, iq2, istate1) = aux
                   end if
                end do
             end do
          end do
       end do
    end do
  end subroutine calculate_3ph_interaction_minimalset

  real(r64) function Vm2_3ph_reference(ev1_s1, ev2_s2, ev3_s3, &
       Index_i, Index_j, Index_k, ifc3, phases_q2q3, ntrip, nb)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.

    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
    complex(r64), intent(in) :: phases_q2q3(ntrip), ev1_s1(nb), ev2_s2(nb), ev3_s3(nb)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)

    !Local variables
    integer(i64) :: it, a, b, c, aind, bind, cind
    complex(r64) :: aux1, aux2, aux3, V0

    aux1 = (0.0_r64, 0.0_r64)
    do it = 1, ntrip
       aind = 3*(Index_k(it) - 1)
       bind = 3*(Index_j(it) - 1)
       cind = 3*(Index_i(it) - 1)
       V0 = (0.0_r64, 0.0_r64)
       do a = 1, 3
          aux2 = conjg(ev3_s3(a + aind))
          do b = 1, 3
             aux3 = aux2*conjg(ev2_s2(b + bind))
             do c = 1, 3
                if(ifc3(c, b, a, it) /= 0.0_r64) then
                   V0 = V0 + ifc3(c, b, a, it)*ev1_s1(c + cind)*aux3
                end if
             end do
          end do
       end do
       aux1 = aux1 + V0*phases_q2q3(it)
    end do

    Vm2_3ph_reference = abs(aux1)**2
  end function Vm2_3ph_reference

end program V3_permutations
