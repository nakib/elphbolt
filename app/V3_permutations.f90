program V3_permutations
   !! generation and classification of irreducible triplets using permutation symmetry
   !!
   !! Returns a mapping of full triplets to a reduced consisting of (ilambda1, ilambda2) pairs
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

contains

   subroutine triplet_test()
      integer(i64) :: nbands, mesh_size(3)
      integer(i64), allocatable :: lambda_1(:), lambda_2(:)
      integer(i64), allocatable :: M(:, :, :)
      integer :: ilambda1, ilambda2

      if (this_image() == 1) then
         print *, "V3 permutations to canonical (lambda1, lambda2, lambda3) form"
         print *, "Number of coarray images =", num_images()
      end if

      nbands = 2
      mesh_size = [2, 2, 2]*1_i64

      allocate(lambda_1(2)); lambda_1 = [1, 2]*1_i64
      allocate(lambda_2(16)); lambda_2 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]*1_i64

      call map_triplet_full_to_reduced(nbands, mesh_size, lambda_1, lambda_2, M)

      print *, "M(ilambda1, ilambda2) → canonical triplet: (lambda1, lambda2, lambda3)"
      do ilambda1 = 1, size(lambda_1)
         do ilambda2 = 1, size(lambda_2)
            write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'M(', ilambda1, ',', ilambda2, ') = (', &
               M(1, ilambda1, ilambda2), ',', M(2, ilambda1, ilambda2), ',', M(3, ilambda1, ilambda2), ')'
         end do
      end do
   end subroutine triplet_test

   subroutine map_triplet_full_to_reduced(nbands, mesh_size, lambda_1, lambda_2, M)
      !! This subroutine maps the canonical triplet is converted back to (ilambda1, ilambda2) and stored in M(:, iq1, iq2).
      !!
      !! M is a 3D array: M(3, size(lambda_1), size(lambda_2))
      !! For each pair (ilambda1, ilambda2), you construct a triplet (lambda1, lambda2, lambda3)
      !! each input lambda index corresponds to a state (iband, ik), where
      !! ik3 is computed from momentum conservation q1 - q2 - q3 = 0 mod mesh
      !! for each triplet, all 3! = 6 permutations are generated using
      !! the Johnson–Trotter algorithm (via the permutations() function), but only the irreducible triplet are conserved

      integer(i64), intent(in) :: nbands, mesh_size(3)
      integer(i64), intent(in) :: lambda_1(:), lambda_2(:)
      integer(i64), allocatable, intent(out) :: M(:, :, :)

      integer(i64) :: ilambda1, ilambda2, iband1, iband2, iband3, ik1, ik2, ik3
      integer(i64) :: q1(3), q2(3), q3(3)
      integer(i64), allocatable :: all_perms(:, :)
      integer(i64) :: triplet_full(3), permuted_triplet(3), canonical_triplet(3)
      integer(i64) :: i, j

      allocate(M(3, size(lambda_1), size(lambda_2)))

      !! Initialize all element of M to a known integer value -1: mean that the triplet has not been assigned yet
      M = -1_i64

      all_perms = permutations(3_i64)

      do ilambda1 = 1, size(lambda_1)

         ! Demux state for lambda_1: for each lambda_1 value, convert its index to (iband1, ik1)
         call demux_state(lambda_1(ilambda1), nbands, iband1, ik1)

         ! Demux vector to get q1
         call demux_vector(ik1, q1, mesh_size, base=0_i64)

         do ilambda2 = 1, size(lambda_2)

            ! Demux state for lambda_2: for each lambda_2 value, convert its index to (iband2, ik2)
            call demux_state(lambda_2(ilambda2), nbands, iband2, ik2)

            ! Demux vector to get q1
            call demux_vector(ik2, q2, mesh_size, base=0_i64)

            ! Compute q3 using modular arithmetic: such that momentum is conserved: q1 - q2 - q3 = 0 mod G
            q3 = modulo(q1 - q2, mesh_size)

            ! Compute ik3 using mux_vector
            ik3 = mux_vector(q3, mesh_size, base=0_i64)

            ! Initialize the canonical triplet with a high value, so any triplet I compare against will be smaller
            canonical_triplet = [huge(0_i64), huge(0_i64), huge(0_i64)]

            ! Iterate over all possible bands for third state: to have (iband3, ik3)
            do iband3 = 1, nbands
               triplet_full = [mux_state(nbands, iband1, ik1), &
                  mux_state(nbands, iband2, ik2), &
                  mux_state(nbands, iband3, ik3)]

               !! Initialize the canonical form of the triplet to the original (unpermuted) triplet, before trying other triplet
               !! triplet_full is a 3x2 array: ((iband1, ik1), (iband2, ik2), (iband3, ik3))
               canonical_triplet = triplet_full

               do i = 1, size(all_perms, 2)
                  do j = 1, 3
                     permuted_triplet(j) = triplet_full(all_perms(j, i))
                  end do

                  !! Use a lexicographic comparison function to keep the smallest permutation
                  if (lex_less_1d(permuted_triplet, canonical_triplet)) canonical_triplet = permuted_triplet
               end do

               !! after elaborating the triplet (lambda1, lambda2, lambda3) as linear state indices
               !! and after generating all 6 permutations
               !! keep only the lexicographically smallest permutation (the irreducible triplet)
               !! Store it in M(:, ilambda1, ilambda2)

               M(:, ilambda1, ilambda2) = canonical_triplet
               exit
            end do
         end do
      end do
   end subroutine map_triplet_full_to_reduced

end program V3_permutations
