program V3_permutations
   !! generation and classification of irreducible triplets using permutation symmetry
   !!
   !! Returns triplets consisting of (iband, ik) pairs
   !! lambda is defined as (iband, i * k), where iband is the band index and ik wave-vector index
   !! calculates all triplet permutations using Johnson and Trotter algorithm
   !! and groups them under lexicographically smallest representatives

   use precision, only: i64, r64
   use misc, only: print_message, subtitle, timer, exit_with_message, mux_vector, demux_state, &
      mux_state, twonorm, permutations, demux_vector, lex_less
   use numerics_module, only: numerics
   use crystal_module, only: crystal
   use symmetry_module, only: symmetry
   use phonon_module, only: phonon

   implicit none

   type Triplet_Set
      !! Stores unique irreducible triplets and their corresponding (iq1, iq2) labels

      integer(i64), allocatable :: canonical_representative(:, :, :)
      integer(i64), allocatable :: iq_pairs(:, :)
   end type Triplet_Set

   type(numerics) :: num
   type(crystal) :: crys
   type(symmetry) :: sym
   type(phonon) :: ph
   !type(timer) :: t_event

   if(this_image() == 1) then
      write(*, '(A)')  'V3 permutations playground'
      write(*, '(A, I5)') 'Number of coarray images = ', num_images()
   end if

   call triplet_test()

contains

   subroutine triplet_test()

      integer(i64) :: nbands, test
      integer(i64), parameter :: mesh_size(3) = [2, 2, 2]
      integer(i64), allocatable :: lambda_1(:), lambda_2(:)
      type(Triplet_Set) :: result
      character(len=100) :: mesh_str

      test = 1
      nbands = 2
      allocate(lambda_1(2)); lambda_1 = [1, 2]
      allocate(lambda_2(16)); lambda_2 = [1, 2, 3, 4, 5, 6, 7, 8, &
         9, 10, 11, 12, 13, 14, 15, 16]

      write(mesh_str, '(3(I0,:,1x))') mesh_size
      write(*,'(A,I0,A,I0,A,A,A)') "Test ", test, " : nbands= ", nbands, " mesh= [", trim(adjustl(mesh_str)), "]"
      call generate_triplets(nbands, mesh_size, lambda_1, lambda_2, result)
      call triplet_result(result)
   end subroutine triplet_test

   subroutine generate_triplets(nbands, mesh_size, lambda_1, lambda_2, triplet_data)
      !! this subroutine generates all valid triplets of the form ((iband1, ik1), (iband2, ik2), (iband3, ik3))
      !!
      !! each input lambda index corresponds to a state (iband, ik), where
      !! ik3 is computed from momentum conservation q1 - q2 - q3 = 0 mod mesh
      !! for each triplet, all 3! = 6 permutations are generated using
      !! the Johnson–Trotter algorithm (via the permutations() function), but only the irreducible triplet are conserved

      integer(i64), intent(in) :: nbands, mesh_size(3)
      integer(i64), intent(in) :: lambda_1(:), lambda_2(:)
      type(Triplet_Set), intent(out) :: triplet_data

      integer(i64) :: ilambda1, ilambda2, iband1, ik1, iband2, ik2, iband3, ik3
      integer(i64) :: q1(3), q2(3), q3(3), m1, m2, base
      integer(i64) :: triplet(3, 2), permuted(3, 2), sorted(3, 2)
      integer(i64), allocatable :: all_perms(:, :), tmp(:, :, :), iq_tmp(:, :)
      integer(i64) :: i, j, k, n_stored, max_num_triplets
      integer(i64) :: perm(3)
      logical :: triplet_exists

      !Give an upper bound on the maximum number of triplets based on mesh and band count
      max_num_triplets = int(1.2_r64*size(lambda_1)*size(lambda_2)*nbands)

      allocate(tmp(3, 2, max_num_triplets))
      allocate(iq_tmp(2, max_num_triplets))

      n_stored = 0

      all_perms = permutations(3_i64)

      do ilambda1 = 1, size(lambda_1)

         ! Demux state for lambda_1: for each lambda_1 value, convert its index to (iband1, ik1)
         call demux_state(lambda_1(ilambda1), nbands, iband1, ik1)

         ! Demux vector to get q1
         call demux_vector(ik1, q1, mesh_size, base = 0_i64)

         do ilambda2 = 1, size(lambda_2)

            ! Demux state for lambda_2: for each lambda_2 value, convert its index to (iband2, ik2)
            call demux_state(lambda_2(ilambda2), nbands, iband2, ik2)

            ! Demux vector to get q1
            call demux_vector(ik2, q2, mesh_size, base = 0_i64)

            ! Compute q3 using modular arithmetic: such that momentum is conserved: q1 - q2 - q3 ≡ 0 mod G
            q3 = modulo(q1 - q2, mesh_size)

            ! Compute ik3 using mux_vector
            ik3 = mux_vector(q3, mesh_size, base = 0_i64)

            ! Iterate over all possible bands for third state: to have (iband3, ik3)
            do iband3 = 1, nbands
               triplet(:, 1) = [iband1, iband2, iband3]
               triplet(:, 2) = [ik1, ik2, ik3]

               sorted = triplet

               ! Generate all permutations using Johnson–Trotter, compare them, and finds the lexicographically smallest triplet
               do i = 1, size(all_perms, 2)
                  perm = all_perms(:, i)

                  do j = 1, 3
                     permuted(j, 1) = triplet(perm(j), 1)
                     permuted(j, 2) = triplet(perm(j), 2)
                  end do

                  if(lex_less(permuted, sorted)) sorted = permuted
               end do

               ! Check and store unique triplets to avoid redundant triplet
               triplet_exists = .false.
               do k = 1, n_stored
                  if (all(sorted == tmp(:, :, k))) then
                     triplet_exists = .true.
                     exit
                  end if
               end do

               if(.not. triplet_exists) then
                  n_stored = n_stored + 1
                  tmp(:, :, n_stored) = sorted
                  iq_tmp(:, n_stored) = [ilambda1, ilambda2]
               end if
            end do
         end do
      end do

      allocate(triplet_data%canonical_representative(3, 2, n_stored))
      allocate(triplet_data%iq_pairs(2, n_stored))
      triplet_data%canonical_representative = tmp(:, :, 1:n_stored)
      triplet_data%iq_pairs = iq_tmp(:, 1:n_stored)

   end subroutine generate_triplets

   subroutine triplet_result(result)
      type(Triplet_Set), intent(in) :: result
      integer(i64) :: k

      do k = 1, size(result%iq_pairs, 2)
         write(*,'(A,"((",I0,",",I0,"),(",I0,",",I0,"),(",I0,",",I0,")) = (",I0,",",I0,")")') &
            'M', result%canonical_representative(1,1,k), result%canonical_representative(1,2,k), &
            result%canonical_representative(2,1,k), result%canonical_representative(2,2,k), &
            result%canonical_representative(3,1,k), result%canonical_representative(3,2,k), &
            result%iq_pairs(1,k), result%iq_pairs(2,k)
      end do
   end subroutine triplet_result

end program V3_permutations
