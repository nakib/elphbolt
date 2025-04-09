program V3_permutations
   !! generation and classification of irreducible triplets using permutation symmetry
   !!
   !! Returns triplets consisting of (iband, ik) pairs
   !! lambda is defined as (iband, i * k), where iband is the band index and ik wave-vector index
   !! calculates all triplet permutations using Johnson and Trotter algorithm
   !! and groups them under lexicographically smallest representatives

   use precision, only: i64, r64
   use misc, only: print_message, subtitle, timer, exit_with_message, mux_vector, demux_state, &
      mux_state, twonorm, permutations, demux_vector
   use numerics_module, only: numerics
   use crystal_module, only: crystal
   use symmetry_module, only: symmetry
   use phonon_module, only: phonon

   implicit none

   type TripletSet
      !! reps is canonical representatives of equivalent permutation triplets (pick 1 unique triplet representation to represent the whole triplets)
      !! counts tells how many permutations belong to reps
      !! stores all permutations that map to reps

      integer(i64), allocatable :: reps(:, :, :)
      integer(i64), allocatable :: counts(:)
      integer(i64), allocatable :: perm_list(:, :, :, :)
   end type TripletSet

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
      implicit none

      integer(i64) :: nbands, test
      integer(i64), parameter :: mesh_size(3) = [2, 2, 2]
      integer(i64), allocatable :: lambda_1(:), lambda_2(:)
      type(TripletSet) :: result
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
      !! demux_state is used to extract (iband, ik) from the 1D lambda index
      !! demux_vector is used to convert ik to a 3D q-vector q(i) on the mesh
      !! ik3 is computed from momentum conservation q1 - q2 + q3 = 0 mod mesh
      !! mux_vector is used to convert q3 back into the 1D index ik3
      !! for each triplet, all 3! = 6 permutations are generated using
      !! the Johnson–Trotter algorithm (via the permutations() function)

      integer(i64), intent(in) :: nbands, mesh_size(3)
      integer(i64), intent(in) :: lambda_1(:), lambda_2(:)
      type(TripletSet), intent(out) :: triplet_data

      integer(i64) :: iq1, iq2, iband1, ik1, iband2, ik2, iband3, ik3
      integer(i64) :: q1(3), q2(3), q3(3), m1, m2, base
      integer(i64) :: triplet(3, 2), permuted(3, 2), sorted(3, 2)
      integer(i64), allocatable :: all_perms(:, :)
      integer(i64) :: i, j, k, perm_index, n_stored, max_triplets
      integer(i64) :: perm(3)
      logical :: new_triplet

      ! the maximum number of irreducible triplets could be stored
      max_triplets = int(1.2d0 * size(lambda_1) * size(lambda_2) * nbands)

      allocate(triplet_data%reps(3, 2, max_triplets))
      allocate(triplet_data%counts(max_triplets))
      allocate(triplet_data%perm_list(3, 2, 6, max_triplets))

      triplet_data%counts = 0_i64
      n_stored = 0

      all_perms = permutations(3_i64)

      do iq1 = 1, size(lambda_1)
         m1 = lambda_1(iq1)
         ! Demux state for lambda_1: for each lambda_1 value, convert its index to (iband1, ik1)
         call demux_state(m1, nbands, iband1, ik1)
         ! Demux vector to get q1
         call demux_vector(ik1, q1, mesh_size, base = 0_i64)

         do iq2 = 1, size(lambda_2)
            m2 = lambda_2(iq2)
            ! Demux state for lambda_2: for each lambda_2 value, convert its index to (iband2, ik2)
            call demux_state(m2, nbands, iband2, ik2)
            ! Demux vector to get q1
            call demux_vector(ik2, q2, mesh_size, base = 0_i64)

            ! Compute q3 using modular arithmetic: such that momentum is conserved: q₁ - q₂ + q₃ ≡ 0
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

               ! Check and Store Unique Triplets to avoid redundant
               new_triplet = .true.
               do k = 1, n_stored
                  if(all(sorted == triplet_data%reps(:, :, k))) then
                     new_triplet = .false.
                     exit
                  end if
               end do

               if(new_triplet) then
                  n_stored = n_stored + 1
                  triplet_data%reps(:, :, n_stored) = sorted
                  triplet_data%counts(n_stored) = 0
               end if

               ! Mapping current triplet permutation to its canonical representative group
               do k = 1, n_stored
                  if(all(sorted == triplet_data%reps(:, :, k))) then
                     perm_index = triplet_data%counts(k) + 1
                     triplet_data%perm_list(:, :, perm_index, k) = triplet
                     triplet_data%counts(k) = perm_index
                     exit
                  end if
               end do
            end do
         end do
      end do
   end subroutine generate_triplets

   pure logical function lex_less(a, b)
      !! lex_less is essential for enforcing permutation symmetry and finding irreducible triplets
      !!
      !! is a comparison function -- boolean comparator
      !! this function compares two triplets a and b lexicographically.
      !! checks if a comes before b when sorted in lexicographic order.
      !! Each a and b is a 3×2 array
      !! lex_less(a, b) returns .true. if a is lexicographically smaller than b.

      integer(i64), intent(in) :: a(3, 2), b(3, 2)

      !Local
      integer :: i

      do i = 1, 3
         if(a(i, 1) < b(i, 1)) then
            lex_less = .true.
            return
         else if(a(i, 1) > b(i, 1)) then
            lex_less = .false.
            return
         else if(a(i, 2) < b(i, 2)) then
            lex_less = .true.
            return
         else if(a(i, 2) > b(i, 2)) then
            lex_less = .false.
            return
         end if
      end do

      lex_less = .false.
   end function lex_less

   subroutine triplet_result(result)
      type(TripletSet), intent(in) :: result
      integer(i64) :: i, k

      do k = 1, size(result%counts)
         if (result%counts(k) == 0) exit
         write(*, '(A, 3("(",I0,",",I0,")",:), A, I0)') 'Irreducible Triplet: ', &
            result%reps(1, 1, k), result%reps(1, 2, k), &
            result%reps(2, 1, k), result%reps(2, 2, k), &
            result%reps(3, 1, k), result%reps(3, 2, k), &
            ' -> Count: ', result%counts(k)
         write(*,*) 'All Permutations:'
         do i = 1, result%counts(k)
            write(*, '(3("(",I0,",",I0,")",:))') result%perm_list(1, 1, i, k), result%perm_list(1, 2, i, k), &
               result%perm_list(2, 1, i, k), result%perm_list(2, 2, i, k), &
               result%perm_list(3, 1, i, k), result%perm_list(3, 2, i, k)
         end do
      end do
   end subroutine triplet_result

end program V3_permutations
