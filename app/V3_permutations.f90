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

  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(phonon) :: ph

  call triplet_test()

  if(this_image() == 1) then
     write(*, '(A)')  'V3 permutations playground'
     write(*, '(A, I5)') 'Number of coarray images = ', num_images()
  end if

  call triplet_test

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
end program V3_permutations
