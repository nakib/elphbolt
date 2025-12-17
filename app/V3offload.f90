program V3offload

#ifdef _OPENACC
  use openacc
#endif

  use precision, only: i64, r64
  use misc, only: print_message, subtitle, timer, exit_with_message, mux_vector, demux_state, &
       mux_state, twonorm, demux_vector, int_div
  use delta, only: delta_fn_triang, delta_fn_tetra
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use phonon_module, only: phonon
  use iso_c_binding, only: c_int
  use gpu_interface, only: compute_V2_on_gpu

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
  real(r64), allocatable :: V2_minus(:, :, :, :, :), V2_plus(:, :, :, :, :)
  integer(c_int) :: istat
  integer(i64) :: count_minus, count_plus

  !gpu side
  real(r64), allocatable :: V2gpu(:)
  integer(i64), allocatable :: S_list(:, :)
  integer(i64) :: S_count
  ! cpu side
  real(r64), allocatable :: V2_cpu(:)
  integer(i64), allocatable :: S_list_cpu(:,:)
  integer(i64) :: S_count_cpu

  real(r64) :: estimated_size_gb, gpu_memory_gb

  ! ESR / delta-function data
  integer(i64), allocatable :: simplex_map(:, :, :)
  integer(i64), allocatable :: simplex_count(:)
  real(r64), allocatable :: simplex_evals(:, :, :)
  real(r64), allocatable :: ens(:, :)
  logical :: use_tetra
  integer(i64) :: nTot

  if(this_image() == 1) then
     write(*, '(A)')  'V3offload playground'
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

  write(*, '(A)') 'Phonon data dimensions:'
  write(*, '(A, I0)') 'ph%numbands = ', ph%numbands
  write(*, '(A, I0)') 'ph%nwv_irred = ', ph%nwv_irred
  write(*, '(A, I0)') 'ph%nwv = ', ph%nwv
  write(*, '(A, I0)') 'ph%numtriplets = ', ph%numtriplets

  if (ph%numbands <= 0 .or. ph%nwv_irred <= 0 .or. ph%nwv <= 0) then
     write(*, '(A)') 'Invalid phonon dimensions!'
     stop
  end if
  estimated_size_gb = real(ph%numbands*ph%nwv*ph%numbands*ph%nwv* &
       ph%nwv_irred*ph%numbands*8, r64) / 1.0e9_r64  !1Gb 10**9 bytes

  print*, 'Estimated full V2 size: ', estimated_size_gb, ' GB'

  gpu_memory_gb = 20.0_r64  

  if(allocated(ph%simplex_map)) then
     ! Use existing ESR data from phonon module
     simplex_map = ph%simplex_map
     simplex_count = ph%simplex_count
     simplex_evals = ph%simplex_evals
     ens = ph%ens
     use_tetra = num%tetrahedra
  else
     allocate(simplex_map(2, ph%nwv, 1))
     allocate(simplex_count(ph%nwv))
     allocate(simplex_evals(1, ph%numbands, 3))  ! triangles (3 vertices)
     allocate(ens(ph%nwv, ph%numbands))

     simplex_map = 1_i64
     simplex_count = 0_i64
     simplex_evals = 0.0_r64

     ! Copy phonon energies if available
     if(allocated(ph%ens)) then
        ens = ph%ens
     else
        ens = 0.0_r64
     end if

     use_tetra = .true.

     write(*, '(A)') 'ESR data not found in phonon module.'
  end if

  !  !Calculate ph-ph vertex (cpu, original)
  !  Vm2_calculator => Vm2_3ph_reference
  !  call t_event%start_timer('reference V- on cpu')
  !  call calculate_3ph_interaction(ph, crys, num, V2, Vm2_calculator)
  !  call t_event%end_timer('reference V- on cpu')
  !  print*, 'value = ', twonorm(pack(V2, .true.))

  ! !  !Calculate ph-ph vertex (cpu, refactor)
  !  Vm2_calculator => Vm2_3ph_refactor
  !  call t_event%start_timer('refactored V- on cpu')
  !  call calculate_3ph_interaction(ph, crys, num, V2, Vm2_calculator)
  !  call t_event%end_timer('refactored V- on cpu')
  !  print*, 'value = ', twonorm(pack(V2, .true.))
  !
  !  !Calculate ph-ph vertex (gpu, algo 1)
  !  call t_event%start_timer('V- on gpu, algo 1')
  !  call calculate_3ph_interaction_gpu(ph, crys, num, V2)
  !  call t_event%end_timer('V- on gpu, algo 1')
  !  print*, 'value = ', twonorm(pack(V2, .true.))
  !
  !  call t_event%start_timer('V- on gpu, low transfer')
  !  call calculate_3ph_interaction_lowtransfer(ph, crys, num, V2)
  !  call t_event%end_timer('V- on gpu, low transfer')
  !  print*, 'value = ', twonorm(pack(V2, .true.))
  !
  !  !Calculate ph-ph vertex (gpu, algo 2)
  !  call t_event%start_timer('V- on gpu, algo 2')
  !  call calculate_3ph_interaction_gpu_algo2(ph, crys, num, V2)
  !  call t_event%end_timer('V- on gpu, algo 2')
  !  print*, 'value = ', twonorm(pack(V2, .true.))

  !    !Calculate ph-ph vertex (V on gpu)- 2 kernel MSR
  !    call t_event%start_timer('V- on gpu, CUDA Fortran (cuf)')
  !    call compute_V2_on_gpu( ph%numbands, ph%nwv_irred, ph%nwv, ph%numtriplets, &
  !                         ph%evecs, ph%Index_i, ph%Index_j, ph%Index_k, ph%ifc3, &
  !                         ph%indexlist_irred, ph%wavevecs, ph%wvmesh, &
  !                         crys%reclattvecs, ph%R_j, ph%R_k, V2, istat)
  !    call t_event%end_timer('V- on gpu, CUDA Fortran (cuf)')
  !    print*, 'value = ', twonorm(pack(V2, .true.))

  !  !Calculate ph-ph vertex (V on cpu to compare with gpu)
  !  call t_event%start_timer('V on cpu')
  !  call compute_V2_cpu(ph%numbands, ph%nwv_irred, ph%nwv, ph%numtriplets, &
  !       ph%evecs, ph%Index_i, ph%Index_j, ph%Index_k, ph%ifc3, &
  !       ph%indexlist_irred, ph%wavevecs, ph%wvmesh, &
  !       crys%reclattvecs,  ph%R_j, ph%R_k, &
  !       simplex_map, simplex_count, simplex_evals, &
  !       ens, use_tetra, V2_cpu, S_list_cpu, S_count_cpu)
  !  call t_event%end_timer('V on cpu')
  !  print*, 'V2 norm = ', twonorm(V2_cpu)

  !Calculate ph-ph vertex (V on gpu)- tripletlist on cpu and Vm2 on gpu
  call t_event%start_timer('Triplet list on cpu and V on gpu, CUDA Fortran (cuf)')
  call compute_V2_on_gpu( ph%numbands, ph%nwv_irred, ph%nwv, ph%numtriplets, &
       ph%evecs, ph%Index_i, ph%Index_j, ph%Index_k, ph%ifc3, &
       ph%indexlist_irred, ph%wavevecs, ph%wvmesh, &
       crys%reclattvecs, ph%R_j, ph%R_k, &
       simplex_map, simplex_count, simplex_evals, ens, use_tetra, &
       V2gpu, S_list, S_count, istat)
  call t_event%end_timer('Triplet list on cpu and V on gpu, CUDA Fortran (cuf)')

  if(istat /= 0) then
     write(*, '(A, I0)') 'GPU computation failed with status = ', istat
  else
     print*, 'Number of valid transitions (S_count) = ', S_count
     print*, 'V2 norm = ', twonorm(V2gpu)
  end if

  ! Cleanup
  if(allocated(V2gpu)) deallocate(V2gpu)
  if(allocated(S_list)) deallocate(S_list)
  deallocate( simplex_map, simplex_count, simplex_evals, ens)

  !  !Calculate ph-ph vertex (V on gpu kernel)-MSR-ESR
  !  call t_event%start_timer('V- on gpu, CUDA Fortran (cuf)')
  !  call compute_V2_on_gpu( &
  !       ph%numbands, ph%nwv_irred, ph%nwv, ph%numtriplets, &
  !       ph%evecs, ph%Index_i, ph%Index_j, ph%Index_k, ph%ifc3, &
  !       ph%indexlist_irred, ph%wavevecs, ph%wvmesh, &
  !       crys%reclattvecs, ph%R_j, ph%R_k, &
  !       ph%ens, V2_minus, V2_plus, istat)
  !  call t_event%end_timer('V- on gpu, CUDA Fortran (cuf)')
  !  print *, 'V2_minus = ', twonorm(pack(V2_minus, .true.))
  !  print *, 'V2_plus = ', twonorm(pack(V2_plus, .true.))
  !
!!$
!!$  !Calculate ph-ph vertex (gpu, algo 3)
!!$  call t_event%start_timer('V- on gpu, algo 3')
!!$  call calculate_3ph_interaction_gpu_algo3(ph, crys, num, V2)
!!$  call t_event%end_timer('V- on gpu, algo 3')
!!$  print*, 'value = ', twonorm(pack(V2, .true.))

!!$  !Calculate ph-ph vertex (gpu, algo 3)
!!$  call t_event%start_timer('V- on gpu, algo 3')
!!$  call calculate_3ph_interaction_gpu_algo3(ph, crys, num, V2)
!!$  call t_event%end_timer('V- on gpu, algo 3')
!!$  print*, 'value = ', twonorm(pack(V2, .true.))

contains

  subroutine compute_V2_cpu(nb, nwv_irred, nwv, ntrip, evecs, Index_i, Index_j, Index_k, ifc3, &
       indexlist_irred, wavevecs, wvmesh, reclatt, Rj, Rk, &
       simplex_map, simplex_count, simplex_evals, ens, use_tetra, &
       V2_cpu, S_list_cpu, S_count_cpu)
    !! CPU reference that builds Slist and compute V2

    integer(i64), intent(in) :: nb, nwv_irred, nwv, ntrip
    complex(r64), intent(in) :: evecs(nwv, nb, nb)
    integer(i64), intent(in) :: Index_i(ntrip), Index_j(ntrip), Index_k(ntrip)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)
    integer(i64), intent(in) :: indexlist_irred(nwv_irred)
    real(r64), intent(in) :: wavevecs(nwv, 3)
    integer(i64), intent(in) :: wvmesh(3)
    real(r64), intent(in) :: reclatt(3, 3)
    real(r64), intent(in) :: Rj(3, ntrip), Rk(3, ntrip)
    integer(i64), intent(in) :: simplex_map(:, :, :)
    integer(i64), intent(in) :: simplex_count(:)
    real(r64), intent(in) :: simplex_evals(:, :, :)
    real(r64), intent(in) :: ens(nwv, nb)
    logical, intent(in) :: use_tetra

    ! for the outputs
    real(r64), allocatable, intent(out) :: V2_cpu(:)
    integer(i64), allocatable, intent(out) :: S_list_cpu(:, :)
    integer(i64), intent(out) :: S_count_cpu

    ! Local variables
    integer(i64) :: nTot, idx, i0, q, r, current_capacity, initial_capacity
    integer(i64) :: iq1_ibz, iq1, iq2, iq3_minus, iq3_plus, s1, s2, s3, i
    integer(i64) :: ilambda1, ilambda2, ilambda3_minus, ilambda3_plus
    real(r64) :: q1f(3), q2f(3), q3f(3), en1, en2, en3_minus, en3_plus
    integer(i64) :: q1i(3), q2i(3), q3i(3), neg_q2i(3), neg_iq2
    real(r64) :: delta_minus, delta_plus
    logical :: has_minus, has_plus
    integer(i64) :: stat_minus, stat_plus, stat_both, stat_neither
    real(r64) :: q2c(3), q3c(3), Vm2_val
    complex(r64) :: ev1(nb), ev2(nb), ev3(nb)
    real(r64) :: time_slist_total, time_v2_total
    real(r64) :: tstart_v2, tend_v2, tstart_all, tend_all
    real(r64) :: tstart, tend, time_s_full, time_all

    S_count_cpu = 0_i64
    stat_minus = 0_i64
    stat_plus = 0_i64
    stat_both = 0_i64
    stat_neither = 0_i64 
    time_slist_total = 0.0_r64
    time_v2_total = 0.0_r64

    nTot = nwv_irred*nwv*nb*nb*nb
    print *, "Total possible transitions:", nTot

    ! Start with 5% of total capacity
    initial_capacity = max(int(0.05_r64*real(nTot, r64), i64), 1000_i64)
    allocate(S_list_cpu(3, initial_capacity))
    allocate(V2_cpu(initial_capacity))

    !S_count_cpu = 0_i64
    current_capacity = initial_capacity

    call cpu_time(tstart_all)
    ! Compute all possible transitions
    do idx = 1_i64, nTot
       i0 = idx - 1_i64

       ! Decompose index into (iq1_ibz, iq2, s1, s2, s3)
       call int_div(i0, nwv_irred, q, r)
       iq1_ibz = r + 1_i64 !because int_div return 0 so +1 for fortran index.

       call int_div(q, nwv, q, r)
       iq2 = r + 1_i64

       call int_div(q, nb, q, r)
       s1 = r + 1_i64

       call int_div(q, nb, q, r)
       s3 = r + 1_i64

       call int_div(q, nb, q, r)
       s2 = r + 1_i64

       ! Map IBZ to FBZ
       iq1 = indexlist_irred(iq1_ibz)

       ! Momentum selection rule: calculate q3
       q1f = wavevecs(iq1, :)
       q2f = wavevecs(iq2, :)

       q1i = int(nint(q1f*real(wvmesh, r64)), i64)
       q2i = int(nint(q2f*real(wvmesh, r64)), i64)

       ! q3- = fold(q1 - q2)
       q3i = modulo(q1i - q2i, wvmesh)
       iq3_minus = mux_vector(q3i, wvmesh, 0_i64)

       ! q3+ = fold(q1 + q2)
       q3i = modulo(q1i + q2i, wvmesh)
       iq3_plus = mux_vector(q3i, wvmesh, 0_i64)

       ! Index of -q2 for plus process
       neg_q2i = modulo(-q2i, wvmesh)
       neg_iq2 = mux_vector(neg_q2i, wvmesh, 0_i64)

       ! Energy selection rule
       en1 = ens(iq1, s1)
       en2 = ens(iq2, s2)
       en3_minus = ens(iq3_minus, s3)
       en3_plus = ens(iq3_plus, s3)

       ! Skip if all energies are zero
       if(en1*en2*en3_minus == 0.0_r64 .and. en1*en2*en3_plus == 0.0_r64) then
          stat_neither = stat_neither + 1_i64
          cycle
       end if

       ! Calculate delta functions for energy conservation
       if(use_tetra) then
          delta_minus = delta_fn_tetra(en1 - en3_minus, iq2, s2, wvmesh, &
               simplex_map, simplex_count, simplex_evals)
          delta_plus = delta_fn_tetra(en3_plus - en1, neg_iq2, s2, wvmesh, &
               simplex_map, simplex_count, simplex_evals)
       else
          delta_minus = delta_fn_triang(en1 - en3_minus, iq2, s2, wvmesh, &
               simplex_map, simplex_count, simplex_evals)
          delta_plus = delta_fn_triang(en3_plus - en1, neg_iq2, s2, wvmesh, &
               simplex_map, simplex_count, simplex_evals)
       end if

       has_minus = (delta_minus > 0.0_r64)
       has_plus = (delta_plus > 0.0_r64)

       ! Count statistics
       if(has_minus .and. has_plus) then
          stat_both = stat_both + 1_i64
          stat_minus = stat_minus + 1_i64
          stat_plus = stat_plus + 1_i64
       else if(has_minus) then
          stat_minus = stat_minus + 1_i64
       else if(has_plus) then
          stat_plus = stat_plus + 1_i64
       else
          stat_neither = stat_neither + 1_i64
       end if

       ! Convert to state indices (ilambda)
       ilambda1 = mux_state(nb, s1, iq1)
       ilambda2 = mux_state(nb, s2, iq2)
       ilambda3_minus = mux_state(nb, s3, iq3_minus)
       ilambda3_plus = mux_state(nb, s3, iq3_plus)

       ! Minus process
       if(has_minus) then
          ! call cpu_time(tstart_slist)
          !      ilambda3_minus = mux_state(nb, s3, iq3_minus)

          S_count_cpu = S_count_cpu + 1_i64

          ! Expand arrays if needed
          if(S_count_cpu > current_capacity) then
             call expand_2d(S_list_cpu, 1.3_r64)
             call expand_1d(V2_cpu, 1.3_r64)
             ! update the capacity.
             current_capacity = size(S_list_cpu, 2)
          end if

          ! Store indices in Slist
          S_list_cpu(1, S_count_cpu) = ilambda1
          S_list_cpu(2, S_count_cpu) = ilambda2
          S_list_cpu(3, S_count_cpu) = ilambda3_minus

          !call cpu_time(tend_slist)
          !time_slist_total = time_slist_total + (tend_slist - tstart_slist)

          !timing V2 computation
          call cpu_time(tstart_v2)

          ! Calculate V2 for this transition
          q3f = real(modulo(q1i - q2i, wvmesh), r64) / real(wvmesh, r64)
          q2c = matmul(reclatt, q2f)
          q3c = matmul(reclatt, q3f)

          ! Gather eigenvectors
          do i = 1, nb
             ev1(i) = evecs(iq1, s1, i)
             ev2(i) = evecs(iq2, s2, i)
             ev3(i) = evecs(iq3_minus, s3, i)
          end do

          ! Compute norm V2
          Vm2_val = Vm2_3ph_cpu(ev1, ev2, ev3, Index_i, Index_j, Index_k, &
               ifc3, q2c, q3c, Rj, Rk, ntrip, nb)

          V2_cpu(S_count_cpu) = Vm2_val
          call cpu_time(tend_v2)
          time_v2_total = time_v2_total + (tend_v2 - tstart_v2)
       end if

       ! Plus process
       if(has_plus) then

          S_count_cpu = S_count_cpu + 1_i64

          ! Expand arrays if needed
          if(S_count_cpu > current_capacity) then
             call expand_2d(S_list_cpu, 1.3_r64)
             call expand_1d(V2_cpu, 1.3_r64)
             current_capacity = size(S_list_cpu, 2)
          end if

          ! Store indices
          S_list_cpu(1, S_count_cpu) = ilambda1
          S_list_cpu(2, S_count_cpu) = ilambda2
          S_list_cpu(3, S_count_cpu) = ilambda3_plus

          !timing V2 computation
          call cpu_time(tstart_v2)

          ! Calculate V2 for this transition
          q3f = real(modulo(q1i + q2i, wvmesh), r64) / real(wvmesh, r64)
          q2c = matmul(reclatt, q2f)
          q3c = matmul(reclatt, q3f)

          ! eigenvectors (nwv, nb, nb)
          do i = 1, nb
             ev1(i) = evecs(iq1, s1, i)
             ev2(i) = evecs(iq2, s2, i)
             ev3(i) = evecs(iq3_plus, s3, i)
          end do

          ! Compute norm V2
          Vm2_val = Vm2_3ph_cpu(ev1, ev2, ev3, Index_i, Index_j, Index_k, &
               ifc3, q2c, q3c, Rj, Rk, ntrip, nb)

          V2_cpu(S_count_cpu) = Vm2_val
          call cpu_time(tend_v2)
          time_v2_total = time_v2_total + (tend_v2 - tstart_v2)
       end if
    end do
    call cpu_time(tend_all)
    time_all = tend_all - tstart_all
    time_s_full = time_all - time_v2_total

    print *, ""
    print *, "=== Timing Results ==="
    print *, "Time for S_list construction (s):", time_s_full
    print *, "Time for V2 computation (s):", time_v2_total
    print *, "Total time (S_list + V2) (s):", time_all
    print *, "Percentage in S_list:", 100.0_r64*time_s_full / (time_s_full + time_v2_total), "%"
    print *, "Percentage in V2:", 100.0_r64*time_v2_total / (time_s_full + time_v2_total), "%"
    print *, "======================"

    print *, ""
    print *, "=== S_list cpu debugg ==="
    print *, "S_list_cpu(1:3, 1):", S_list_cpu(:, 1)
    if(S_count_cpu >= 100) then
       print *, "S_list_cpu(1:3, 100):", S_list_cpu(:, 100)
    end if
    print *, "S_list_cpu(1:3, S_count_cpu):", S_list_cpu(:, S_count_cpu)
    print *, "Min ilambda1:", minval(S_list_cpu(1, 1:S_count_cpu))
    print *, "Max ilambda1:", maxval(S_list_cpu(1, 1:S_count_cpu))
    print *, "Min ilambda2:", minval(S_list_cpu(2, 1:S_count_cpu))
    print *, "Max ilambda2:", maxval(S_list_cpu(2, 1:S_count_cpu))
    print *, "Min ilambda3:", minval(S_list_cpu(3, 1:S_count_cpu))
    print *, "Max ilambda3:", maxval(S_list_cpu(3, 1:S_count_cpu))
    print *, "Expected max ilambda:", nwv*nb
    print *, "=============================="

    print *, ""
    print *, "Statistics"
    print *, "Total transitions checked:", nTot
    print *, "Neither satisfied:", stat_neither
    print *, "Only minus satisfied:", stat_minus - stat_both
    print *, "Only plus satisfied:", stat_plus - stat_both
    print *, "Both satisfied:", stat_both
    print *, "Total minus processes:", stat_minus
    print *, "Total plus processes:", stat_plus
    print *, "Expected S_count:", stat_minus + stat_plus
    print *, "Actual S_count:", S_count_cpu
    print *, "Discrepancy:", S_count_cpu - (stat_minus + stat_plus)
    print *, ""

    ! Shrink to actual size
    call shrink_2d(S_list_cpu, S_count_cpu)
    call shrink_1d(V2_cpu, S_count_cpu)
  end subroutine compute_V2_cpu

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
         idim, jdim, s2s3
    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         aux
    complex(r64) :: phases(ph%numtriplets)

    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))

    V2 = 0.0

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
             do s2s3 = 1, ph%numbands**2
                s2 = int((s2s3 - 1)/ph%numbands) + 1 !changes slow
                s3 = modulo(s2s3 - 1, ph%numbands) + 1 !changes fast

                aux = Vm2_calculator(ph%evecs(iq1, s1, :), &
                     ph%evecs(iq2, s2, :), ph%evecs(iq3_minus, s3, :), &
                     ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:, :, :, :), &
                     phases(:), ph%numtriplets, ph%numbands)

                V2(s3, iq3_minus, s2, iq2, istate1) = aux
             end do
          end do
       end do
    end do
  end subroutine calculate_3ph_interaction

  subroutine calculate_3ph_interaction_gpu(ph, crys, num, V2)
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    real(r64), allocatable, intent(out) :: V2(:, :, :, :, :)

    !Local variables
    integer(i64) :: istate1, nstates_irred, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
         idim, jdim, s2s3
    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         aux
    complex(r64) :: phases(ph%numtriplets)

    complex(r64) :: R1(3, 3, ph%numtriplets), R3(ph%numtriplets)

    !$acc data copyin(ph%ifc3, ph%Index_i, ph%Index_j, ph%Index_k) &
    !$acc      create(R1, R3, ev1, ev2, ev3)

    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))

    V2 = 0.0

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

             ev1 = reshape(ph%evecs(iq1, s1, :), shape = [3, int(crys%numatoms, 4)])
             !$acc update device(ev1)

             do s2 = 1, ph%numbands
                ev2 = reshape(ph%evecs(iq2, s2, :), shape = [3, int(crys%numatoms, 4)])
                !$acc update device(ev2)

                do s3 = 1, ph%numbands
                   ev3 = reshape(ph%evecs(iq3_minus, s3, :), shape = [3, int(crys%numatoms, 4)])
                   !$acc update device(ev3)

                   aux = Vm2_3ph_gpu(ev1, ev2, ev3, &
                        ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:,:,:,:), &
                        phases(:), ph%numtriplets, ph%numbands, R1, R3)

                   V2(s3, iq3_minus, s2, iq2, istate1) = aux
                end do
             end do
          end do
       end do
    end do
    !$acc end data
  end subroutine calculate_3ph_interaction_gpu

  subroutine calculate_3ph_interaction_lowtransfer(ph, crys, num, V2)
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    real(r64), allocatable, intent(out) :: V2(:, :, :, :, :)

    !Local variables
    integer(i64) :: istate1, nstates_irred, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
         idim, jdim, s2s3
    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         aux
    complex(r64) :: phases(ph%numtriplets), &
         evecs(ph%nwv, ph%numbands, 3, crys%numatoms)

    complex(r64) :: R1(3, 3, ph%numtriplets), R3(ph%numtriplets)


    !TODO Reshape full eigenvectors list
    print*, 'reshaping evecs'
    do iq1 = 1, ph%nwv
       do s1 = 1, ph%numbands
          evecs(iq1, s1, :, :) = reshape(ph%evecs(iq1, s1, :), &
               shape = [3, int(crys%numatoms, 4)])
       end do
    end do
    print*, 'reshaping done'

    print*, 'copying data to gpu'

    !$acc data copyin(ph%ifc3, ph%Index_i, ph%Index_j, ph%Index_k, evecs) &
    !$acc      create(iq1, s1, iq2, s2, iq3_minus, s3, R1, R3)

    print*, 'copying data done'

    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))

    V2 = 0.0

    !Run over first phonon IBZ states
    !do istate1 = 1, nstates_irred
    !Demux state index into branch (s) and wave vector (iq) indices
    !call demux_state(istate1, ph%numbands, s1, iq1_ibz)
    do iq1_ibz = 1, ph%nwv_irred

       !Muxed index of wave vector from the IBZ index list.
       !This will be used to access IBZ information from the FBZ quantities.
       iq1 = ph%indexlist_irred(iq1_ibz)

       !$acc update device(iq1)

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

          !$acc update device(iq2, iq3_minus)

          q2_cart = matmul(crys%reclattvecs, q2)
          q3_minus_cart = matmul(crys%reclattvecs, q3_minus)

          phases = exp((0.0_r64, -1.0_r64)* &
               (matmul(q2_cart, ph%R_j) + matmul(q3_minus_cart, ph%R_k)))

          do s1 = 1, ph%numbands
             !$acc update device(s1)

             istate1 = mux_state(ph%numbands, s1, iq1_ibz)

             do s2 = 1, ph%numbands
                !$acc update device(s2)

                do s3 = 1, ph%numbands
                   !$acc update device(s3)

!!$acc update device(iq1, iq2, iq3_minus, s1, s2, s3)

                   aux = Vm2_3ph_gpu_lowtransfer(evecs, &
                        iq1, s1, iq2, s2, iq3_minus, s3, &
                        ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:,:,:,:), &
                        phases(:), ph%numtriplets, ph%numbands, R1, R3)

                   V2(s3, iq3_minus, s2, iq2, istate1) = aux
                end do
             end do
          end do
       end do
    end do
    !$acc end data
  end subroutine calculate_3ph_interaction_lowtransfer

  real(r64) function Vm2_3ph_gpu_lowtransfer(evecs, &
       iq1, s1, iq2, s2, iq3_minus, s3, &
       Index_i, Index_j, Index_k, ifc3, phases, ntrip, nb, R1, R3)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.

    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), &
         Index_k(ntrip), nb, iq1, s1, iq2, s2, iq3_minus, s3
    complex(r64), intent(in) :: phases(ntrip), evecs(:, :, :, :)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)
    complex(r64), intent(inout) :: R1(3, 3, ntrip), R3(ntrip)

    !Local variables
    integer(i64) :: it, a, b, c, aind, bind, cind
    integer(i64) :: ijk, nijk, ndim, idim
    complex(r64) :: aux1, aux2, ev1(3), ev2(3), ev3(3)

    nijk = ntrip
    ndim = 3

    !$acc data create(ev1, ev2, ev3)

    !$acc parallel loop
    do ijk = 1, nijk       
       ev1(:) = evecs(iq1, s1, :, Index_i(ijk))
       ev2(:) = evecs(iq2, s2, :, Index_j(ijk))
       ev3(:) = evecs(iq3_minus, s3, :, Index_k(ijk))

       do c = 1, ndim
          do b = 1, ndim
             aux1 = (0.0_r64, 0.0_r64)
             do idim = 1, ndim
                !aux1 = aux1 + ifc3(idim, b, c, ijk)*evecs(iq1, s1, idim, Index_i(ijk))
                aux1 = aux1 + ifc3(idim, b, c, ijk)*ev1(idim)
             end do
             R1(b, c, ijk) = aux1
          end do
       end do

       do c = 1, ndim
          aux1 = (0.0_r64, 0.0_r64)
          do idim = 1, ndim
             !aux1 = aux1 + R1(idim, c, ijk)*conjg(evecs(iq2, s2, idim, Index_j(ijk)))
             aux1 = aux1 + R1(idim, c, ijk)*conjg(ev2(idim))
          end do
          R1(c, 1, ijk) = aux1
       end do

       R3(ijk) = (0.0_r64, 0.0_r64)
       do idim = 1, ndim
          R3(ijk) = R3(ijk) + R1(idim, 1, ijk)*&
               conjg(ev3(idim))
          !conjg(evecs(iq3_minus, s3, idim, Index_k(ijk)))
       end do
    end do
    !$acc update host(R3)

    !$acc end data

    !And the Fourier transform
    Vm2_3ph_gpu_lowtransfer = abs(dot_product(conjg(R3), phases))**2
  end function Vm2_3ph_gpu_lowtransfer

!!$  subroutine calculate_3ph_interaction_gpu_flat(ph, crys, num, V2)
!!$    type(phonon), intent(in) :: ph
!!$    type(crystal), intent(in) :: crys
!!$    type(numerics), intent(in) :: num
!!$    real(r64), allocatable, intent(out) :: V2(:, :, :, :, :)
!!$
!!$    !Local variables
!!$    integer(i64) :: istate1, nstates_irred, &
!!$         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
!!$         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
!!$         idim, jdim, s2s3
!!$    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), &
!!$         q3_minus_cart(3), aux, ifc3_flat(27*ph%numtriplets)
!!$    
!!$    !complex(r64) :: ev1(3, crys%numatoms), ev2(3, crys%numatoms), ev3(3, crys%numatoms)
!!$    complex(r64) :: phases(ph%numtriplets)
!!$
!!$    !complex(r64) :: R3(ph%numtriplets)
!!$    !complex(r64) :: ev12(3, 3, ph%numtriplets)
!!$    !complex(r64) :: ev123(3, 3, 3, ph%numtriplets)
!!$    complex(r64) :: ev123_flat(27*ph%numtriplets)
!!$
!!$    ifc3_flat = pack(ph%ifc3, mask = .true.)
!!$    
!!$    !$acc data copyin(ifc3_flat, ph%Index_i, ph%Index_j, ph%Index_k, ph%evecs) &
!!$    !$acc      create(ev123_flat)
!!$
!!$    !Total number of IBZ blocks states
!!$    nstates_irred = ph%nwv_irred*ph%numbands
!!$
!!$    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))
!!$
!!$    V2 = 0.0
!!$
!!$    !Run over first phonon IBZ states
!!$    !do istate1 = 1, nstates_irred
!!$    !Demux state index into branch (s) and wave vector (iq) indices
!!$    !call demux_state(istate1, ph%numbands, s1, iq1_ibz)
!!$    do iq1_ibz = 1, ph%nwv_irred
!!$
!!$       !Muxed index of wave vector from the IBZ index list.
!!$       !This will be used to access IBZ information from the FBZ quantities.
!!$       iq1 = ph%indexlist_irred(iq1_ibz)
!!$
!!$       !Initial (IBZ blocks) wave vector (crystal coords.)
!!$       q1 = ph%wavevecs(iq1, :)
!!$
!!$       !Convert from crystal to 0-based index vector
!!$       q1_indvec = nint(q1*ph%wvmesh)
!!$
!!$       do iq2 = 1, ph%nwv
!!$          !Initial (IBZ blocks) wave vector (crystal coords.)
!!$          q2 = ph%wavevecs(iq2, :)
!!$
!!$          !Convert from crystal to 0-based index vector
!!$          q2_indvec = nint(q2*ph%wvmesh)
!!$
!!$          !Folded final phonon wave vector
!!$          q3_minus_indvec = modulo(q1_indvec - q2_indvec, ph%wvmesh) !0-based index vector
!!$          q3_minus = q3_minus_indvec/dble(ph%wvmesh) !crystal coords.
!!$
!!$          !Muxed index of q3_minus
!!$          iq3_minus = mux_vector(q3_minus_indvec, ph%wvmesh, 0_i64)
!!$
!!$          q2_cart = matmul(crys%reclattvecs, q2)
!!$          q3_minus_cart = matmul(crys%reclattvecs, q3_minus)
!!$
!!$          phases = exp((0.0_r64, -1.0_r64)* &
!!$               (matmul(q2_cart, ph%R_j) + matmul(q3_minus_cart, ph%R_k)))
!!$
!!$          do s1 = 1, ph%numbands
!!$
!!$             istate1 = mux_state(ph%numbands, s1, iq1_ibz)
!!$
!!$             !ev1 = reshape(ph%evecs(iq1, s1, :), shape = [3, crys%numatoms])
!!$             !!$acc update device(ev1)
!!$
!!$             do s2 = 1, ph%numbands
!!$                !ev2 = reshape(ph%evecs(iq2, s2, :), shape = [3, crys%numatoms])                     
!!$                !!$acc update device(ev2)
!!$
!!$                do s3 = 1, ph%numbands
!!$                   !ev3 = reshape(ph%evecs(iq3_minus, s3, :), shape = [3, crys%numatoms])
!!$                   !!$acc update device(ev3)
!!$
!!$                   aux = Vm2_3ph_gpu_flat(s1, iq1, iq2, s2, iq3_minus, s3, &
!!$                        ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), &
!!$                        phases(:), ph%numtriplets, ph%numbands)
!!$
!!$                   V2(s3, iq3_minus, s2, iq2, istate1) = aux
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$    !$acc end data
!!$  end subroutine calculate_3ph_interaction_gpu_flat
!!$
!!$  real(r64) function Vm2_3ph_gpu_flat(s1, iq1, iq2, s2, iq3_minus, s3, &
!!$       Index_i, Index_j, Index_k, phases, ntrip, nb)
!!$    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.
!!$
!!$    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), &
!!$         nb, s1, iq1, iq2, s2, iq3_minus, s3
!!$    complex(r64), intent(in) :: phases(ntrip)
!!$
!!$    !Local variables
!!$    integer(i64) :: it, a, b, c, aind, bind, cind, numatoms, idx
!!$    integer(i64) :: ijk, nijk, ndim, idim, j, atom1_ind, atom2_ind, atom3_ind
!!$    complex(r64) :: aux1, aux2, R3(ntrip)
!!$
!!$    nijk = ntrip
!!$    ndim = 3
!!$    numatoms = nb/3
!!$
!!$    R3 = (0.0_r64, 0.0_r64)
!!$
  !acc data copyin(iq1, iq2, iq3_minus, s1, s2, s3) copy(R3)
!!$    !$acc data copy(R3)
!!$
!!$    !$acc parallel loop
!!$    do ijk = 1, nijk
!!$       !0-based
!!$       atom1_ind = 3*(Index_i(ijk) - 1)
!!$       atom2_ind = 3*(Index_j(ijk) - 1)
!!$       atom3_ind = 3*(Index_k(ijk) - 1)
!!$       do c = 1, 3
!!$          do b = 1, 3
!!$             do a = 1, 3
!!$                idx = a + (b - 1)*3 + (c - 1)*9 + (ijk - 1)*27
!!$
!!$                ev123_flat(idx) = ph%evecs(iq1, s1, atom1_ind + a)*&
!!$                     conjg(ph%evecs(iq2, s2, atom2_ind + b))*&
!!$                     conjg(ph%evecs(iq3_minus, s3, atom3_ind + c))
!!$             end do
!!$          end do
!!$       end do
!!$
!!$       idx = (ijk - 1)*27
!!$       R3(ijk) = R3(ijk) + &
!!$            dot_product(ifc3_flat(idx + 1:idx + 27), ev123_flat(idx + 1:idx + 27))
!!$    end do
!!$    !$acc end data
!!$
!!$    !And the Fourier transform
!!$    Vm2_3ph_gpu_flat = abs(dot_product(conjg(R3), phases))**2
!!$  end function Vm2_3ph_gpu_flat

  subroutine calculate_3ph_interaction_gpu_algo2(ph, crys, num, V2)
    type(phonon), intent(in) :: ph
    type(crystal), intent(in) :: crys
    type(numerics), intent(in) :: num
    real(r64), allocatable, intent(out) :: V2(:, :, :, :, :)

    !Local variables
    integer(i64) :: istate1, nstates_irred, &
         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
         idim, jdim, s2s3
    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
         aux
    complex(r64) :: phases(ph%numtriplets)

    complex(r64) :: R1(3, 3, ph%numtriplets), R3(ph%numtriplets)

    !$acc data copyin(ph%ifc3, ph%Index_i, ph%Index_j, ph%Index_k, ph%R_j, ph%R_k) &
    !$acc      create(R1, R3, ev1, ev2, ev3, phases, q2_cart, q3_minus_cart)

    !Total number of IBZ blocks states
    nstates_irred = ph%nwv_irred*ph%numbands

    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))

    V2 = 0.0

    !Run over first phonon IBZ states
    do istate1 = 1, nstates_irred
       !Demux state index into branch (s) and wave vector (iq) indices
       call demux_state(istate1, ph%numbands, s1, iq1_ibz)

       !Muxed index of wave vector from the IBZ index list.
       !This will be used to access IBZ information from the FBZ quantities.
       iq1 = ph%indexlist_irred(iq1_ibz)

       !Initial (IBZ blocks) wave vector (crystal coords.)
       q1 = ph%wavevecs(iq1, :)

       !Convert from crystal to 0-based index vector
       q1_indvec = nint(q1*ph%wvmesh)

       ev1 = reshape(ph%evecs(iq1, s1, :), shape = [3, int(crys%numatoms, 4)])
       !$acc update device(ev1)

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
          !$acc update device(q2_cart, q3_minus_cart)

          call calculate_phases_on_gpu(q2_cart, q3_minus_cart, ph%R_j, ph%R_k, phases, ph%numtriplets)

          do s2 = 1, ph%numbands
             ev2 = reshape(ph%evecs(iq2, s2, :), shape = [3, int(crys%numatoms, 4)])
             !$acc update device(ev2)

             do s3 = 1, ph%numbands
                ev3 = reshape(ph%evecs(iq3_minus, s3, :), shape = [3, int(crys%numatoms, 4)])
                !$acc update device(ev3)

                aux = Vm2_3ph_gpu_algo2(ev1, ev2, ev3, &
                     ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:, :, :, :), &
                     phases, ph%numtriplets, ph%numbands, R1, R3)

                V2(s3, iq3_minus, s2, iq2, istate1) = aux
             end do
          end do
       end do
    end do
    !$acc end data

    !end associate
  end subroutine calculate_3ph_interaction_gpu_algo2
!!$
!!$  subroutine calculate_3ph_interaction_gpu_algo3(ph, crys, num, V2)
!!$    type(phonon), intent(in) :: ph
!!$    type(crystal), intent(in) :: crys
!!$    type(numerics), intent(in) :: num
!!$    real(r64), allocatable, intent(out) :: V2(:, :, :, :, :)
!!$
!!$    !Local variables
!!$    integer(i64) :: istate1, nstates_irred, &
!!$         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
!!$         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
!!$         idim, jdim, s2s3
!!$    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
!!$         aux
!!$    complex(r64) :: phases(ph%numtriplets)
!!$
!!$    complex(r64) :: R1(3, 3, ph%numtriplets), R3(ph%numtriplets)
!!$
!!$    !$acc data copyin(ph%ifc3, ph%Index_i, ph%Index_j, ph%Index_k) &
!!$    !$acc      create(R1, R3, ev1, ev2, ev3)
!!$
!!$    !Total number of IBZ blocks states
!!$    nstates_irred = ph%nwv_irred*ph%numbands
!!$
!!$    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))
!!$
!!$    V2 = 0.0
!!$
!!$    !Run over first phonon IBZ states
!!$    !do istate1 = 1, nstates_irred
!!$    !Demux state index into branch (s) and wave vector (iq) indices
!!$    !call demux_state(istate1, ph%numbands, s1, iq1_ibz)
!!$    do iq1_ibz = 1, ph%nwv_irred
!!$
!!$       !Muxed index of wave vector from the IBZ index list.
!!$       !This will be used to access IBZ information from the FBZ quantities.
!!$       iq1 = ph%indexlist_irred(iq1_ibz)
!!$
!!$       !Initial (IBZ blocks) wave vector (crystal coords.)
!!$       q1 = ph%wavevecs(iq1, :)
!!$
!!$       !Convert from crystal to 0-based index vector
!!$       q1_indvec = nint(q1*ph%wvmesh)
!!$
!!$       do iq2 = 1, ph%nwv
!!$          !Initial (IBZ blocks) wave vector (crystal coords.)
!!$          q2 = ph%wavevecs(iq2, :)
!!$
!!$          !Convert from crystal to 0-based index vector
!!$          q2_indvec = nint(q2*ph%wvmesh)
!!$
!!$          !Folded final phonon wave vector
!!$          q3_minus_indvec = modulo(q1_indvec - q2_indvec, ph%wvmesh) !0-based index vector
!!$          q3_minus = q3_minus_indvec/dble(ph%wvmesh) !crystal coords.
!!$
!!$          !Muxed index of q3_minus
!!$          iq3_minus = mux_vector(q3_minus_indvec, ph%wvmesh, 0_i64)
!!$
!!$          q2_cart = matmul(crys%reclattvecs, q2)
!!$          q3_minus_cart = matmul(crys%reclattvecs, q3_minus)
!!$
!!$          phases = exp((0.0_r64, -1.0_r64)* &
!!$               (matmul(q2_cart, ph%R_j) + matmul(q3_minus_cart, ph%R_k)))
!!$
!!$          do s1 = 1, ph%numbands
!!$
!!$             istate1 = mux_state(ph%numbands, s1, iq1_ibz)
!!$
!!$             ev1 = reshape(ph%evecs(iq1, s1, :), shape = [3, int(crys%numatoms, 4)])
!!$             !$acc update device(ev1)
!!$
!!$             do s2 = 1, ph%numbands
!!$                ev2 = reshape(ph%evecs(iq2, s2, :), shape = [3, int(crys%numatoms, 4)])
!!$                !$acc update device(ev2)
!!$
!!$                do s3 = 1, ph%numbands
!!$                   ev3 = reshape(ph%evecs(iq3_minus, s3, :), shape = [3, int(crys%numatoms, 4)])
!!$                   !$acc update device(ev3)
!!$
!!$                   aux = Vm2_3ph_gpu_algo3(ev1, ev2, ev3, &
!!$                        ph%Index_i(:), ph%Index_j(:), ph%Index_k(:), ph%ifc3(:,:,:,:), &
!!$                        phases(:), ph%numtriplets, ph%numbands, R1, R3)
!!$
!!$                   V2(s3, iq3_minus, s2, iq2, istate1) = aux
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$    !$acc end data
!!$  end subroutine calculate_3ph_interaction_gpu_algo3

!!$  subroutine calculate_3ph_interaction_gpu_algo3(ph, crys, num, V2)
!!$    type(phonon), intent(in) :: ph
!!$    type(crystal), intent(in) :: crys
!!$    type(numerics), intent(in) :: num
!!$    real(r64), allocatable, intent(out) :: V2(:, :, :, :, :)
!!$
!!$    !Local variables
!!$    integer(i64) :: istate1, nstates_irred, &
!!$         nprocs, s1, s2, s3, iq1_ibz, iq1, iq2, iq3_minus, it, &
!!$         q1_indvec(3), q2_indvec(3), q3_minus_indvec(3), &
!!$         idim, jdim, s2s3, q3_list(ph%nwv, ph%nwv_irred), numtriplets_gpu, numbands_gpu, nwv_gpu, &
!!$         aind, bind, cind, a, b, c
!!$    real(r64) :: en1, en2, en3, q1(3), q2(3), q3_minus(3), q2_cart(3), q3_minus_cart(3), &
!!$         aux
!!$    real(r64), allocatable :: V2_device(:, :, :, :)
!!$    complex(r64) :: phases(ph%numtriplets), ev1(ph%numbands), ev2(ph%numbands), ev3(ph%numbands), &
!!$         V0, aux1, aux2, aux3
!!$
!!$    !complex(r64) :: R1(3, 3, ph%numtriplets), R3(ph%numtriplets)
!!$
!!$    !Total number of IBZ blocks states
!!$    nstates_irred = ph%nwv_irred*ph%numbands
!!$
!!$    numtriplets_gpu = ph%numtriplets
!!$    numbands_gpu = ph%numbands
!!$    nwv_gpu = ph%nwv
!!$
!!$    !Precompute list of q3s
!!$    do iq1_ibz = 1, ph%nwv_irred
!!$       !Initial (IBZ blocks) wave vector (crystal coords.)
!!$       q1 = ph%wavevecs(ph%indexlist_irred(iq1_ibz), :)
!!$
!!$       !Convert from crystal to 0-based index vector
!!$       q1_indvec = nint(q1*ph%wvmesh)
!!$
!!$       do iq2 = 1, ph%nwv
!!$          !Initial (IBZ blocks) wave vector (crystal coords.)
!!$          q2 = ph%wavevecs(iq2, :)
!!$
!!$          !Convert from crystal to 0-based index vector
!!$          q2_indvec = nint(q2*ph%wvmesh)
!!$
!!$          !Folded final phonon wave vector
!!$          q3_minus_indvec = modulo(q1_indvec - q2_indvec, ph%wvmesh) !0-based index vector
!!$          q3_minus = q3_minus_indvec/dble(ph%wvmesh) !crystal coords.
!!$
!!$          !Muxed index of q3_minus
!!$          iq3_minus = mux_vector(q3_minus_indvec, ph%wvmesh, 0_i64)
!!$
!!$          q3_list(iq2, iq1_ibz) = iq3_minus
!!$       end do
!!$    end do
!!$    !
!!$
!!$    allocate(V2(ph%numbands, ph%nwv, ph%numbands, ph%nwv, nstates_irred))
!!$    allocate(V2_device(ph%numbands, ph%nwv, ph%numbands, ph%nwv))
!!$
!!$    associate(wavevecs=>ph%wavevecs, evecs=>ph%evecs, R_j=>ph%R_j, R_k=>ph%R_k, &
!!$         Index_i=>ph%Index_i, Index_j=>ph%Index_j, Index_k=>ph%Index_k, ifc3=>ph%ifc3)
!!$
!!$      !$acc data copyin(ifc3, Index_i, Index_j, Index_k, R_j, R_k, &
!!$      !$acc             wavevecs, q3_list, evecs, numtriplets_gpu, numbands_gpu, nwv_gpu) &
!!$      !$acc      create(ev1, ev2, ev3, phases, q2_cart, q3_minus_cart, V2_device)
!!$
!!$      !Run over first phonon IBZ states
!!$      do istate1 = 1, nstates_irred
!!$         !Demux state index into branch (s) and wave vector (iq) indices
!!$         call demux_state(istate1, ph%numbands, s1, iq1_ibz)
!!$
!!$         !Muxed index of wave vector from the IBZ index list.
!!$         !This will be used to access IBZ information from the FBZ quantities.
!!$         iq1 = ph%indexlist_irred(iq1_ibz)
!!$
!!$         !$acc data copyin(iq1, s1)
!!$         !$acc parallel loop
!!$         do iq2 = 1, nwv_gpu
!!$            !1st phonon eigenvector
!!$            ev1 = evecs(iq1, s1, :)
!!$
!!$            !Muxed index of q3_minus
!!$            iq3_minus = q3_list(iq2, iq1_ibz)
!!$
!!$            q2_cart = matmul(crys%reclattvecs, wavevecs(iq2, :))
!!$            q3_minus_cart = matmul(crys%reclattvecs, wavevecs(iq3_minus, :))
!!$
!!$            do it = 1, numtriplets_gpu
!!$               phases(it) = exp((0.0_r64, -1.0_r64)* &
!!$                    (dot_product(q2_cart, R_j(:, it)) + &
!!$                    dot_product(q3_minus_cart, R_k(:, it))))
!!$            end do
!!$
!!$            do s2 = 1, numbands_gpu
!!$               ev2 = evecs(iq2, s2, :)
!!$
!!$               do s3 = 1, numbands_gpu
!!$                  ev3 = evecs(iq3_minus, s3, :)
!!$
!!$                  V2_device(s3, iq3_minus, s2, iq2) = &
!!$                       Vm2_3ph_reference(ev1, ev2, ev3, &
!!$                       Index_i(:), Index_j(:), Index_k(:), ifc3(:,:,:,:), &
!!$                       phases, numtriplets_gpu, numbands_gpu)
!!$               end do
!!$            end do
!!$         end do
!!$         !$acc update host(V2_device)
!!$         !$acc end data
!!$
!!$         V2(:, :, :, :, istate1) = V2_device(:, :, :, :)
!!$      end do
!!$      !$acc end data
!!$
!!$    end associate
!!$  end subroutine calculate_3ph_interaction_gpu_algo3

  real(r64) function Vm2_3ph_reference(ev1_s1, ev2_s2, ev3_s3, &
       Index_i, Index_j, Index_k, ifc3, phases_q2q3, ntrip, nb)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.

    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
    complex(r64), intent(in) :: phases_q2q3(ntrip), ev1_s1(nb), ev2_s2(nb), ev3_s3(nb)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)

    !Local variables
    integer(i64) :: it, a, b, c, aind, bind, cind
    complex(r64) :: aux1, aux2, aux3, V0

    !$acc routine seq

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

  real(r64) function Vm2_3ph_refactor(ev1_s1, ev2_s2, ev3_s3, &
       Index_i, Index_j, Index_k, ifc3, phases_q2q3, ntrip, nb)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.

    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
    complex(r64), intent(in) :: phases_q2q3(ntrip), ev1_s1(nb), ev2_s2(nb), ev3_s3(nb)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)

    !Local variables
    integer(i64) :: it, a, b, c, aind, bind, cind
    integer(i64) :: i, j, k, ijk, nijk, ndim, nat
    complex(r64) :: aux1, aux2, aux3, V0
    complex(r64) :: ev1(3, nb/3), ev2(3, nb/3), ev3(3, nb/3)
    complex(r64) :: R1(3, 3, ntrip), R2(3, ntrip), R3(ntrip)

    nijk = ntrip
    ndim = 3
    nat = nb/3
    ev1 = reshape(ev1_s1, shape = [ndim, nat])
    ev2 = reshape(ev2_s2, shape = [ndim, nat])
    ev3 = reshape(ev3_s3, shape = [ndim, nat])

    !The tensor contraction section
    do ijk = 1, nijk
       i = Index_i(ijk)
       j = Index_j(ijk)
       k = Index_k(ijk)

       do concurrent (c = 1:ndim, b = 1:ndim)
          R1(b, c, ijk) = dot_product(ifc3(:, b, c, ijk), ev1(:, i))
       end do

       do c = 1, ndim
          R2(c, ijk) = dot_product(ev2(:, j), R1(:, c, ijk))
       end do

       R3(ijk) = dot_product(ev3(:, k), R2(:, ijk))
    end do

    !And the Fourier transform
    Vm2_3ph_refactor = abs(dot_product(conjg(R3), phases_q2q3))**2
  end function Vm2_3ph_refactor

  real(r64) function Vm2_3ph_gpu(ev1, ev2, ev3, &
       Index_i, Index_j, Index_k, ifc3, phases, ntrip, nb, R1, R3)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.

    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
    complex(r64), intent(in) :: phases(ntrip), ev1(3, nb/3), ev2(3, nb/3), ev3(3, nb/3)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)
    complex(r64), intent(inout) :: R1(3, 3, ntrip), R3(ntrip)

    !Local variables
    integer(i64) :: it, a, b, c, aind, bind, cind
    integer(i64) :: ijk, nijk, ndim, idim
    complex(r64) :: aux1, aux2

    nijk = ntrip
    ndim = 3

    !$acc parallel loop
    do ijk = 1, nijk
       do c = 1, ndim
          do b = 1, ndim
             aux1 = (0.0_r64, 0.0_r64)
             do idim = 1, ndim
                aux1 = aux1 + ifc3(idim, b, c, ijk)*ev1(idim, Index_i(ijk))
             end do
             R1(b, c, ijk) = aux1
          end do
       end do

       do c = 1, ndim
          aux1 = (0.0_r64, 0.0_r64)
          do idim = 1, ndim
             aux1 = aux1 + R1(idim, c, ijk)*conjg(ev2(idim, Index_j(ijk)))
          end do
          R1(c, 1, ijk) = aux1
       end do

       R3(ijk) = (0.0_r64, 0.0_r64)
       do idim = 1, ndim
          R3(ijk) = R3(ijk) + R1(idim, 1, ijk)*conjg(ev3(idim, Index_K(ijk)))
       end do
    end do
    !$acc update host(R3)

    !And the Fourier transform
    Vm2_3ph_gpu = abs(dot_product(conjg(R3), phases))**2
  end function Vm2_3ph_gpu

  real(r64) function Vm2_3ph_gpu_algo2(ev1, ev2, ev3, &
       Index_i, Index_j, Index_k, ifc3, phases, ntrip, nb, R1, R3)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.

    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
    complex(r64), intent(in) :: phases(ntrip), ev1(3, nb/3), ev2(3, nb/3), ev3(3, nb/3)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)
    complex(r64), intent(inout) :: R1(3, 3, ntrip), R3(ntrip)

    !Local variables
    integer(i64) :: it, a, b, c, aind, bind, cind
    integer(i64) :: ijk, nijk, ndim, idim
    complex(r64) :: aux1, aux2

    nijk = ntrip
    ndim = 3

    !$acc parallel loop
    do ijk = 1, nijk
       do concurrent (c = 1:ndim, b = 1:ndim)
          aux1 = (0.0_r64, 0.0_r64)
          do idim = 1, ndim
             aux1 = aux1 + ifc3(idim, b, c, ijk)*ev1(idim, Index_i(ijk))
          end do
          R1(b, c, ijk) = aux1
       end do

       do c = 1, ndim
          aux1 = (0.0_r64, 0.0_r64)
          do idim = 1, ndim
             aux1 = aux1 + R1(idim, c, ijk)*conjg(ev2(idim, Index_j(ijk)))
          end do
          R1(c, 1, ijk) = aux1
       end do

       R3(ijk) = (0.0_r64, 0.0_r64)
       do idim = 1, ndim
          R3(ijk) = R3(ijk) + R1(idim, 1, ijk)*conjg(ev3(idim, Index_k(ijk)))
       end do

       R3(ijk) = R3(ijk)*phases(ijk)
    end do
    !$acc update host(R3)

    Vm2_3ph_gpu_algo2 = abs(sum(R3))**2
  end function Vm2_3ph_gpu_algo2
!!$
!!$  real(r64) function Vm2_3ph_gpu_algo3(ev1, ev2, ev3, &
!!$       Index_i, Index_j, Index_k, ifc3, phases, ntrip, nb, R1, R3)
!!$    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.
!!$
!!$    integer(i64), intent(in) :: ntrip, Index_i(ntrip), Index_j(ntrip), Index_k(ntrip), nb
!!$    complex(r64), intent(in) :: phases(ntrip), ev1(3, nb/3), ev2(3, nb/3), ev3(3, nb/3)
!!$    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip)
!!$    complex(r64), intent(inout) :: R1(3, 3, ntrip), R3(ntrip)
!!$
!!$    !Local variables
!!$    integer(i64) :: it, a, b, c, iat, jnd, knd
!!$    integer(i64) :: ijk, nijk, ndim, idim
!!$    complex(r64) :: aux1, aux2
!!$
!!$    nijk = ntrip
!!$    ndim = 3
!!$
!!$    !$acc parallel loop gang vector private(aux1, idim, iat, jnd, knd)
!!$    do ijk = 1, nijk
!!$       !$acc cache(ev1, ev2, ev3)
!!$       iat = Index_i(ijk)
!!$       jnd = Index_j(ijk)
!!$       knd = Index_k(ijk)
!!$
!!$       !R1(b, c, ijk)
!!$       do c = 1, ndim
!!$          do b = 1, ndim
!!$             aux1 = (0.0_r64, 0.0_r64)
!!$             do idim = 1, ndim
!!$                aux1 = aux1 + ifc3(idim, b, c, ijk)*ev1(idim, iat)
!!$             end do
!!$             R1(b, c, ijk) = aux1
!!$          end do
!!$       end do
!!$
!!$       !Contrat with ev2
!!$       do c = 1, ndim
!!$          aux1 = (0.0_r64, 0.0_r64)
!!$          do idim = 1, ndim
!!$             aux1 = aux1 + R1(idim, c, ijk)*conjg(ev2(idim, jnd))
!!$          end do
!!$          R1(c, 1, ijk) = aux1
!!$       end do
!!$
!!$       !Contract with ev3
!!$       R3(ijk) = (0.0_r64, 0.0_r64)
!!$       do idim = 1, ndim
!!$          R3(ijk) = R3(ijk) + R1(idim, 1, ijk)*conjg(ev3(idim, Index_K(ijk)))
!!$       end do
!!$
!!$       !Fourier phase
!!$       R3(ijk) = R3(ijk)*phases(ijk)
!!$    end do
!!$    !$acc update host(R3)
!!$
!!$    !And the Fourier transform
!!$    Vm2_3ph_gpu_algo3 = abs(sum(R3))**2
!!$  end function Vm2_3ph_gpu_algo3
!!$
  subroutine calculate_phases_on_gpu(q2_cart, q3_minus_cart, R_j, R_k, phases, ntrip)
    real(r64), intent(in) :: q2_cart(3), q3_minus_cart(3)
    real(r64), intent(in) :: R_j(3, ntrip), R_k(3, ntrip)
    complex(r64), intent(out) :: phases(ntrip)
    integer(i64), intent(in) :: ntrip

    integer(i64) :: it

    !$acc parallel loop
    do it = 1, ntrip
       phases(it) = exp((0.0_r64, -1.0_r64)* &
            (dot_product(q2_cart, R_j(:, it)) + &
            dot_product(q3_minus_cart, R_k(:, it))))
    end do
  end subroutine calculate_phases_on_gpu
!!$
!!$  subroutine outer_complex(A, B, C)
!!$    !! Outer product of A and B
!!$    !!
!!$    !! C_ij = A_i.B_j
!!$
!!$    complex(r64), intent(in) :: A(:), B(:)
!!$    complex(r64), intent(out) :: C(:, :) !2D-array
!!$    integer :: j, len_a, len_b
!!$
!!$    !$acc routine seq
!!$    len_a = size(A)
!!$    len_b = size(B)
!!$    if(len_a /= size(C, 1) &
!!$         .or. len_b /= size(C, 2)) then
!!$       print *, 'Dimension mismatch. Exiting.'
!!$       call exit
!!$    end if
!!$
!!$    do j = 1, len_b
!!$       C(:, j) = A(:)*B(j)
!!$    end do
!!$  end subroutine outer_complex

  subroutine shrink_2d(arr, actual_count)
    !! Shrink 2D integer array to actual used size
    !!
    !! Keeps first dimension, shrinks second dimension to actual_count

    integer(i64), allocatable, intent(inout) :: arr(:, :) !2D: (6 indices per transition, S_count transitions)
    integer(i64), intent(in) :: actual_count
    integer(i64) :: dim1
    integer(i64), allocatable :: tmparr(:,:)

    dim1 = size(arr, 1)

    ! Allocate smaller array
    ! what does in my case? allocate(tmp(3, S_count))
    !tmp(:, 1:S_count) = S_list(:, 1:S_count)
    !call move_alloc(tmp, S_list), that means during the copy, I temporarily have: old S_list and new tmp.
    allocate(tmparr(dim1, actual_count))

    ! Copy only the used portion
    tmparr = arr(:, 1:actual_count)
    call move_alloc(tmparr, arr)

    print *, "Shrunk S_list to actual size:", actual_count
  end subroutine shrink_2d

  subroutine expand_2d(arr, factor)
    !! Expand 2D integer array along second dimension by given factor
    !!
    !! First dimension size is preserved

    integer(i64), allocatable, intent(inout) :: arr(:, :) !2D: (6 indices per transition, S_count transitions)
    real(r64), intent(in) :: factor
    integer(i64) :: dim1, dim2, new_dim2
    integer(i64), allocatable :: tmparr(:, :)
    real(r64) :: expansion_factor

    ! Calculate new size
    dim1 = size(arr, 1)
    dim2 = size(arr, 2)
    new_dim2 = int(real(dim2, r64)*factor, i64)

    ! Allocate new larger array
    allocate(tmparr(dim1, new_dim2))
    tmparr = 0_i64
    tmparr(:, 1:dim2) = arr

    ! Move allocation
    call move_alloc(tmparr, arr)

    print *, "Expanded S_list from", dim2, "to", new_dim2
  end subroutine expand_2d

  subroutine expand_1d(arr, factor)
    !! Expand 1D real array by given factor

    real(r64), allocatable, intent(inout) :: arr(:)
    real(r64), intent(in) :: factor
    integer(i64) :: dim1, new_dim1
    real(r64), allocatable :: tmparr(:)

    dim1 = size(arr)
    new_dim1 = int(real(dim1, r64)*factor, i64)

    allocate(tmparr(new_dim1))
    tmparr = 0.0_r64
    tmparr(1:dim1) = arr

    call move_alloc(tmparr, arr)

    print *, "Expanded V2_cpu from", dim1, "to", new_dim1
  end subroutine expand_1d

  subroutine shrink_1d(arr, actual_count)
    !! Shrink 1D real array to actual used size

    real(r64), allocatable, intent(inout) :: arr(:)
    integer(i64), intent(in) :: actual_count
    real(r64), allocatable :: tmparr(:)

    allocate(tmparr(actual_count))
    tmparr = arr(1:actual_count)
    call move_alloc(tmparr, arr)

    print *, "Shrunk V2_cpu to actual size:", actual_count
  end subroutine shrink_1d

  pure real(r64) function Vm2_3ph_cpu(ev1_s1, ev2_s2, ev3_s3, &
       Index_i, Index_j, Index_k, ifc3, q2c, q3c, Rj, Rk, ntrip, nb)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.

    integer(i64), intent(in) :: ntrip, nb
    integer(i64), intent(in) :: Index_i(ntrip), Index_j(ntrip), Index_k(ntrip)
    complex(r64), intent(in) :: ev1_s1(nb), ev2_s2(nb), ev3_s3(nb)
    real(r64), intent(in) :: ifc3(3, 3, 3, ntrip), q2c(3), q3c(3)
    real(r64), intent(in) :: Rj(3, ntrip), Rk(3, ntrip)

    integer(i64) :: it, a, b, c, aind, bind, cind
    complex(r64) :: aux1, aux2, aux3, V0, phase
    real(r64) :: arg

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
       arg = dot_product(q2c, Rj(:, it)) + dot_product(q3c, Rk(:, it))
       phase = exp((0.0_r64, -1.0_r64)*arg)
       aux1 = aux1 + V0*phase
    end do
    Vm2_3ph_cpu = abs(aux1)**2
  end function Vm2_3ph_cpu
end program V3offload
