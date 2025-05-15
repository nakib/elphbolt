program test_V3_permutation

#ifdef _OPENACC
   use openacc
#endif

   use precision, only: i64, r64
   use misc, only: print_message, subtitle, timer, exit_with_message, mux_vector, demux_state, &
      mux_state, twonorm, demux_vector, permutations, lex_less_1d, map_triplet_full_to_reduced
   use V3offload, only: calculate_3ph_interaction, calculate_3ph_interaction_minimalset, Vm2_3ph_reference
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
   integer :: count_full, count_minimal
   integer(i64) :: lambda1, lambda2, lambda3
   integer(i64) :: istate1, iq2, iq3_minus, s2, s3, nstates_irred, s1, iq1_ibz, &
      iq1
   real(r64) :: val !value

   if(this_image() == 1) then
      write(*, '(A)')  'V3offload playground'
      write(*, '(A, I5)') 'Number of coarray images = ', num_images()
   end if

   !Set up crystal
   !call crys%initialize

   !Set up numerics data
   !call num%initialize(crys)

   !Calculate crystal and BZ symmetries
   !call sym%calculate_symmetries(crys, num%qmesh)

   !Calculate phonons
   !call ph%initialize(crys, sym, num)

   !Calculate ph-ph vertex (cpu, original)
   call  Vm2_3ph_reference(ev1_s1, ev2_s2, ev3_s3, Index_i, Index_j, Index_k, ifc3, phases_q2q3, ntrip, nb)
   Vm2_calculator => Vm2_3ph_reference
   call t_event%start_timer('reference V- on cpu')
   call calculate_3ph_interaction(ph, crys, num, V2, Vm2_calculator)
   call t_event%end_timer('reference V- on cpu')
   print*, 'value = ', twonorm(pack(V2, .true.))
   ! print*, V2(3, 1, 2, 3, 3)
   ! print*, V2(3, 2, 4, 3, 2)
   ! print*, V2(4, 1, 4, 1, 4)
   ! print*, V2(5, 1, 5, 1, 5)
   ! print*, V2(6, 1, 6, 1, 6)
   count_full = 0
   print *, '---------------------------------------------------------------'
   print *, '   lambda1    lambda2    lambda3      Value'
   print *, '---------------------------------------------------------------'

   do lambda1 = 1, 6
      call demux_state(lambda1, ph%numbands, s1, iq1)

      do lambda2 = 1, 10
         call demux_state(lambda2, ph%numbands, s2, iq2)

         do lambda3 = 1, 10
            call demux_state(lambda3, ph%numbands, s3, iq3_minus)

            if(s3 <= size(V2, 1) .and. iq3_minus <= size(V2, 2) .and. &
               s2 <= size(V2, 3) .and. iq2 <= size(V2, 4) .and. &
               lambda1 <= size(V2, 5)) then

               val = V2(s3, iq3_minus, s2, iq2, lambda1)

               !if(abs(val) > 1.0e-7_r64) then
               write(*,'(3I10, 2X, F16.12)') lambda1, lambda2, lambda3, val
               count_full = count_full + 1
               !end if
            end if
         end do
      end do
   end do
   print *, 'Number of V2 elements:', count_full

   !Calculate V2_minimal_set (reference)
   call  Vm2_3ph_reference(ev1_s1, ev2_s2, ev3_s3, Index_i, Index_j, Index_k, ifc3, phases_q2q3, ntrip, nb)
   Vm2_calculator => Vm2_3ph_reference
   call t_event%start_timer('reference V2 minimal set')
   call calculate_3ph_interaction_minimalset(ph, crys, num, V2_minimal_set, M, Vm2_calculator)
   call t_event%end_timer('reference V2 minimal set')
   !print*, 'value = ', twonorm(pack(V2_minimal_set, .true.))
   ! print*, V2(3, 1, 2, 3, 3)
   ! print*, V2(3, 2, 4, 3, 2)
   ! print*, V2(4, 1, 4, 1, 4)
   ! print*, V2(5, 1, 5, 1, 5)
   ! print*, V2(6, 1, 6, 1, 6)
   count_minimal = 0
   print *, '---------------------------------------------------------------'
   print *, '   lambda1    lambda2    lambda3      Value'
   print *, '---------------------------------------------------------------'

   do lambda1 = 1, 6
      call demux_state(lambda1, ph%numbands, s1, iq1)

      do lambda2 = 1, 10
         call demux_state(lambda2, ph%numbands, s2, iq2)

         do lambda3 = 1, 10
            call demux_state(lambda3, ph%numbands, s3, iq3_minus)

            !canonical value from V2_minimal_set
            if(s3 <= size(V2_minimal_set, 1) .and. iq3_minus <= size(V2_minimal_set, 2) .and. &
               s2 <= size(V2_minimal_set, 3) .and. iq2 <= size(V2_minimal_set, 4) .and. &
               lambda1 <= size(V2_minimal_set, 5)) then

               val = V2_minimal_set(s3, iq3_minus, s2, iq2, lambda1)

               if(abs(val) /= 1.0_r64) then
                  write(*,'(3I10, 2X, F16.12)') lambda1, lambda2, lambda3, val
                  count_minimal = count_minimal + 1
               end if
            end if
         end do
      end do
   end do
   print *, 'Number of V2 minimal set elements:', count_minimal

   ! Reduction factor between V2 and V2 minimal set
   print *, 'Reduction factor (minimal/full):', real(count_minimal)/real(count_full)
   print *, 'Symmetry saving (in %):', (1.0 - real(count_minimal)/real(count_full)) * 100.0

end program test_V3_permutation
