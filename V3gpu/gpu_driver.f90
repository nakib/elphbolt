module gpu_3ph_driver
  use iso_fortran_env, only: int64, real64
  use cudafor
  use gpu_3ph_kernels !contains: build_Slist_host (host), calculate_Vm2_kernel, mux_state ...

  implicit none

  private
  public :: compute_V2_on_gpu

contains

  subroutine expand_int64_2d(arr, factor)
    !! Expand 2D integer array along second dimension by given factor
    !!
    !! First dimension size is preserved

    integer(int64), allocatable, intent(inout) :: arr(:, :) !2D: (3 indices per transition, triplet_count transitions)
    real(real64), intent(in) :: factor
    integer(int64) :: dim1, dim2, new_dim2
    integer(int64), allocatable :: tmparr(:, :)
    real(real64) :: expansion_factor

    ! Calculate new size
    dim1 = size(arr, 1)
    dim2 = size(arr, 2)
    new_dim2 = int(real(dim2, real64)*factor, int64)

    ! Allocate new larger array
    allocate(tmparr(dim1, new_dim2))
    tmparr = 0_int64
    tmparr(:, 1:dim2) = arr

    ! Move allocation
    call move_alloc(tmparr, arr)

    print *, "Expanded S_list from", dim2, "to", new_dim2
  end subroutine expand_int64_2d

  subroutine shrink_int64_2d(arr, actual_count)
    !! Shrink 2D integer array to actual used size
    !!
    !! Keeps first dimension, shrinks second dimension to actual_count

    integer(int64), allocatable, intent(inout) :: arr(:, :) !2D: (3 indices per transition, triplet_count transitions)
    integer(int64), intent(in) :: actual_count
    integer(int64) :: dim1
    integer(int64), allocatable :: tmparr(:,:)

    dim1 = size(arr, 1)

    ! Allocate smaller array
    ! what does in my case? allocate(tmp(3, triplet_count))
    !tmp(:, 1:triplet_count) = triplet_list(:, 1:triplet_count)
    !call move_alloc(tmp, triplet_list), that means during the copy, I temporarily have: old triplet_list and new tmp.
    allocate(tmparr(dim1, actual_count))

    ! Copy only the used portion
    tmparr = arr(:, 1:actual_count)
    call move_alloc(tmparr, arr)

    print *, "Shrunk S_list to actual size:", actual_count
  end subroutine shrink_int64_2d

  subroutine compute_V2_on_gpu(nb, nwv_irred, nwv, ntrip, &
       evecs, Index_i, Index_j, Index_k, ifc3, &
       indexlist_irred, wavevecs, wvmesh, reclatt, Rj, Rk, &
       simplex_map, simplex_count, simplex_evals, ens, use_Tetrahedra, &
       V2, triplet_list_out, triplet_count_out, istat)
    !! compute_V2_on_gpu
    !!
    !! Build triplet list on host side, compact the triplet list of valid tuples (with delta functions)
    !! then compute Vm2 per tuple and store in 1D V2 array.

    ! Arrays (host)
    integer(int64), intent(in) :: nb, nwv_irred, nwv, ntrip
    complex(real64), intent(in) :: evecs(nwv, nb, nb)
    integer(int64), intent(in) :: Index_i(ntrip), Index_j(ntrip), Index_k(ntrip)
    real(real64), intent(in) :: ifc3(3, 3, 3, ntrip)
    integer(int64), intent(in) :: indexlist_irred(nwv_irred)
    real(real64), intent(in) :: wavevecs(nwv, 3)
    integer(int64), intent(in) :: wvmesh(3)
    real(real64), intent(in) :: reclatt(3, 3)
    real(real64), intent(in) :: Rj(3, ntrip), Rk(3, ntrip)

    ! Parameters for energy selection (host side)
    integer(int64), intent(in) :: simplex_map(:, :, :)  ! triangle or tetra map
    integer(int64), intent(in) :: simplex_count(:)      ! counts per k-point
    real(real64), intent(in) :: simplex_evals(:, :, :)  ! simplex energies
    real(real64), intent(in) :: ens(nwv, nb)            ! phonon energies
    logical, intent(in) :: use_Tetrahedra               ! true=tetra, false=triangle

    !final result(host) - V2 1D sized to triplet_count
    real(real64), allocatable, intent(out) :: V2(:)

    !Each column of triplet_list_out describes one valid three-phonon process, using state indices (ilambda).
    integer(int64), allocatable, intent(out) :: triplet_list_out(:, :)
    integer(int64), intent(out) :: triplet_count_out

    integer, intent(out) :: istat

    ! Derived sizes
    integer(int64) :: nstates_irred
    integer(int64) :: maxtriplet, Thread_per_block, nBlocks, nTot, initial_size
    integer :: ierr, ist
    real(real64) :: tstart_gpu, tend_gpu, time_gpu_kernel, tstart_cpu, tend_cpu, time_cpu, &
         t_all0, t_all1, t_total, memory_GB
    integer(kind = cuda_count_kind) :: heapsize

    !Host side that contains all valid processes, but physically has extra unused capacity.
    integer(int64), allocatable :: triplet_list_h(:, :)

    ! Device arrays
    complex(real64), device, allocatable :: evecs_d(:, :, :)
    real(real64), device, allocatable :: ifc3_d(:, :, :, :)
    integer(int64), device, allocatable :: Index_i_d(:), Index_j_d(:), Index_k_d(:)
    integer(int64), device, allocatable :: indexlist_irred_d(:)
    real(real64), device, allocatable :: wavevecs_d(:, :)
    integer(int64), device, allocatable :: wvmesh_d(:)
    real(real64), device, allocatable :: reclatt_d(:, :)
    real(real64), device, allocatable :: Rj_d(:, :), Rk_d(:, :)
    real(real64), device, allocatable :: V2_d(:)

    !Device (GPU) copy of triplet_list_h. Used by CUDA kernel to map (ilambda1, ilambda2, ilambda3) into indices into V2_d...
    integer(int64), device, allocatable :: triplet_list_d(:, :) !(3, maxS)
    real(real64), device, allocatable :: Vm2_list_d(:)

    ! Host side
    integer(int64) :: triplet_count_h, neither_count
    integer(int64) :: both_count, only_minus, only_plus, expected

    real(real64), parameter :: frac_valid_ref = 0.032_real64 ! from Si 30^3, valid fraction = 3.17%
    real(real64), parameter :: safety_factor = 2.0_real64
    real(real64), parameter :: max_tripletlist_GB = 6.0_real64 ! cap, you can modify it. 

    istat = 0

    ! Increase heap size for dynamic allocation
    heapsize = 512_cuda_count_kind*1024_cuda_count_kind*1024_cuda_count_kind
    ist = cudaDeviceSetLimit(cudaLimitMallocHeapSize, heapsize)

    print*, "Setting heap size to:", real(heapsize, real64) / 1.0e9_real64, "GB"

    if(ist /= 0) then
       write(*, *) 'cudaDeviceSetLimit(heap) failed, ist=', ist
    end if

    ! Total number of IBZ block-states
    nstates_irred = nwv_irred*nb
    ! Total possible processes
    nTot = nwv_irred*nwv*nb*nb*nb

    ! based on the observed 3% valid process
    initial_size = int(frac_valid_ref*safety_factor*real(nTot, real64), int64)   
    memory_GB = real(3_int64*initial_size*8_int64, real64) / 1.0e9_real64
    print *, "Initial triplet_list processes:", initial_size
    print *, "Initial triplet_list memory estimate: ~", memory_GB, "GB"

    if(memory_GB > max_tripletlist_GB) then
       print *, "triplet_list would use", memory_GB, "GB, capping to", max_tripletlist_GB, "GB"
       initial_size = int(max_tripletlist_GB*1.0e9_real64 / real(3_int64*8_int64, real64), int64)
       memory_GB = real(3_int64*initial_size*8_int64, real64) / 1.0e9_real64
       print *, "Capped initial_size:", initial_size, " to ~", memory_GB, "GB"
    end if

    allocate(triplet_list_h(3, initial_size))

    print*, 'triplet_list allocation done.'

    ! Allocate & copy to device
    allocate(evecs_d(nwv, nb, nb))
    allocate(ifc3_d(3, 3, 3, ntrip))
    allocate(Index_i_d(ntrip), Index_j_d(ntrip), Index_k_d(ntrip))
    allocate(indexlist_irred_d(nwv_irred))
    allocate(wavevecs_d(nwv, 3))
    allocate(wvmesh_d(3))
    allocate(reclatt_d(3, 3))
    allocate(Rj_d(3, ntrip), Rk_d(3, ntrip))

    print*, 'Other stuff allocation done.'

    evecs_d = evecs
    ifc3_d = ifc3
    Index_i_d = Index_i
    Index_j_d = Index_j
    Index_k_d = Index_k
    indexlist_irred_d = indexlist_irred
    wavevecs_d = wavevecs
    wvmesh_d = wvmesh
    reclatt_d = reclatt
    Rj_d = Rj
    Rk_d = Rk

    call cpu_time(t_all0)
    !build triplet_list on cpu with MSR and ESR    
    call cpu_time(tstart_cpu)
    call build_tripletlist_host_dynamic_without_mask( triplet_list_h, triplet_count_h, &
         indexlist_irred, wavevecs, wvmesh, reclatt, &
         simplex_map, simplex_count, simplex_evals, ens, use_Tetrahedra, &
         nwv_irred, nwv, nb )

    print*, 'build_triplet list call finished'

    call cpu_time(tend_cpu)
    time_cpu = tend_cpu - tstart_cpu
    print *, 'CPU build Slist on host took (s): ', time_cpu

    print *, 'Total possible transitions: ', nTot
    print *, 'Valid tuples (triplet_count_h): ', triplet_count_h
    print*, "Invalid transitions filtered:", nTot - triplet_count_h
    print*, "Percentage kept (valid):", &
         100.0_real64*real(triplet_count_h, real64) / real(nTot, real64), "%"
    print*, "Percentage filtered out:", &
         100.0_real64*real(nTot - triplet_count_h, real64) / real(nTot, real64), "%"
    print*, "Reduction factor:", &
         real(nTot, real64) / real(triplet_count_h, real64), "x"

    if(triplet_count_h > 2_int64*nTot) then
       print *, "gpu_3ph_driver: triplet_list (triplet_count_h=", triplet_count_h, ")"
       istat = -1 ! signals that "Triplet list overflow"
       deallocate(triplet_list_h)
       return !Must exit on error, prevent GPU corruption.
    end if

    if(triplet_count_h <= 0_int64) then
       allocate(V2(0))
       allocate(triplet_list_out(3, 0))
       triplet_count_out = 0
       deallocate(triplet_list_h)
       return ! Skip GPU: no valid transitions, nothing to compute.
    end if

    !device side
    allocate(triplet_list_d(3, triplet_count_h))
    triplet_list_d = triplet_list_h(:, 1:triplet_count_h)

    print*, "=== S_list gpu debugg ==="
    print*, "triplet_list_h(1:3, 1):", triplet_list_h(:, 1)
    print*, "triplet_list_h(1:3, 100):", triplet_list_h(:, 100)
    print*, "triplet_list_h(1:3, triplet_count_h):", triplet_list_h(:, triplet_count_h)
    print*, "Min ilambda1:", minval(triplet_list_h(1, 1:triplet_count_h))
    print*, "Max ilambda1:", maxval(triplet_list_h(1, 1:triplet_count_h))
    print*, "Min ilambda2:", minval(triplet_list_h(2, 1:triplet_count_h))
    print*, "Max ilambda2:", maxval(triplet_list_h(2, 1:triplet_count_h))
    print*, "Min ilambda3:", minval(triplet_list_h(3, 1:triplet_count_h))
    print*, "Max ilambda3:", maxval(triplet_list_h(3, 1:triplet_count_h))
    print*, "Expected max ilambda:", nwv*nb
    print*, "=========================="

    ! Allocate 1D V2 arrays (device and host)
    allocate(V2_d(triplet_count_h))
    allocate(Vm2_list_d(triplet_count_h))
    V2_d = 0.0_real64

    Thread_per_block = 256_int64
    nBlocks = int((triplet_count_h + Thread_per_block - 1_int64) / Thread_per_block)

    print*, "Launching kernel with:"
    print*, "Blocks:", nBlocks
    print*, "Thread_per_block", Thread_per_block
    print*, "Total threads:", nBlocks*Thread_per_block
    print*, "triplet_count_h:", triplet_count_h

    ! check kernel lauch eror
    !source: https://forums.developer.nvidia.com/t/newbie-question-about-maximum-number-of-blocks/42015
    if(nBlocks > 2147483647_int64) then
       print *, "nBlocks exceeds cuda maximum!" !2^31-1 blocks (the first dimension maximum) 
       istat = -2 ! is an error code that signals "too many blocks > 2^31-1" to the calling function.
       call dealloc_all()
       return
    end if

    ! Kernel: compute Vm2 over triplet_list and store directly in V2_d
    call cpu_time(tstart_gpu)
    call calculate_Vm2_kernel<<<nBlocks, Thread_per_block>>>( &
         V2_d, evecs_d, triplet_list_d, &
         Index_i_d, Index_j_d, Index_k_d, ifc3_d, indexlist_irred_d, wavevecs_d, &
         wvmesh_d, reclatt_d, Rj_d, Rk_d, triplet_count_h, nb, nwv, ntrip)

    ierr = cudaGetLastError()
    if(ierr /= cudaSuccess) then
       print *, "===================================="
       print *, "Error: Kernel launch failed!"
       print *, "Error code:", ierr
       print *, "This means the kernel did not start"
       print *, "===================================="
       istat = ierr
       call dealloc_all()
       return
    end if

    print*, 'V2 kernel launched (waiting for completion...)'

    ierr = cudaDeviceSynchronize()
    if(ierr /= 0) then
       print *, "===================================="
       print *, "Error: Kernel execution failed!"
       print *, "Error code:", ierr
       print *, "The kernel started but crashed"
       print *, "===================================="
       istat = ierr
       call dealloc_all()
       return
    end if
    !    print*, 'V2 kernel computed.'
    !
    !    ierr = cudaDeviceSynchronize()
    !    if(ierr /= 0) then
    !       istat = ierr
    !       call dealloc_all()
    !       return !GPU error occurred, cleanup and propagate error code. Must exit, results are invalid.
    !    end if

    call cpu_time(tend_gpu)
    time_gpu_kernel = tend_gpu - tstart_gpu
    print*, "Elapsed time gpu kernel for calculate Vm2:", time_gpu_kernel

    ! Copy back to host   ! I think here another problem
    allocate(V2(triplet_count_h))
    V2 = V2_d
    ierr = cudaDeviceSynchronize()
    if(ierr /= 0) then
       print *, "Copy-back synchronization failed:", ierr
    end if

    ! debug evecs and wavecs
    print*, "V2(1) =", V2(1)
    print*, "V2(100) =", V2(100)
    print*, "V2(triplet_count) =", V2(triplet_count_h)
    print*, "V2 array size:", size(V2)
    print*, "Memory used by V2 (GB):", real(size(V2)*8_int64, real64) / 1.0e9_real64
    print*, "Number of zeros:", count(V2 == 0.0_real64)
    print*, "Number of positive values:", count(V2 > 0.0_real64)

    if(allocated(triplet_list_out)) deallocate(triplet_list_out)
    ! move_alloc does not allocate or copy data, just transfers the allocation descriptor from triplet_list_h to triplet_list_out.
    call move_alloc(triplet_list_h, triplet_list_out)
    triplet_count_out = triplet_count_h

    print *, "triplet_list_out shape after move_alloc:", size(triplet_list_out, 1), size(triplet_list_out, 2)
    print *, "Logical triplet_count_out:", triplet_count_out !can see it’s still (3, 250000000) physically but I know the meaningful part is (3, 189240633).

    call cpu_time(t_all1)
    t_total = t_all1 - t_all0
    print *, "Total compute_V2_on_gpu (s):", t_total

    ! Clean up & return
    call dealloc_all()
    return ! results computed and copied to host, all temp memory freed.but this can be removed, execution would exit naturally here anyway.

  contains
    subroutine dealloc_all()
      deallocate(evecs_d, V2_d, ifc3_d, Index_i_d, Index_j_d, Index_k_d, indexlist_irred_d, wavevecs_d, &
           wvmesh_d, reclatt_d, Rj_d, Rk_d, triplet_list_d, Vm2_list_d)
    end subroutine dealloc_all

    subroutine build_tripletlist_host_dynamic_without_mask( triplet_list, triplet_count, indexlist_irred, wavevecs, wvmesh, reclatt, &
         simplex_map, simplex_count, simplex_evals, ens, use_Tetrahedra, &
         nwv_irred, nwv, nb )
      !! This subroutine build triplet list with momentum and energy selection rules 
      !!
      !! expand and shrink when needed. 

      integer(int64), value :: nwv_irred, nwv, nb
      logical, value :: use_Tetrahedra

      ! inputs
      integer(int64), intent(in) :: indexlist_irred(nwv_irred)
      real(real64), intent(in) :: wavevecs(nwv, 3), reclatt(3, 3), ens(nwv, nb)
      integer(int64), intent(in) :: wvmesh(3)
      integer(int64), intent(in) :: simplex_map(:, :, :), simplex_count(:)
      real(real64), intent(in) :: simplex_evals(:, :, :)

      ! outputs
      integer(int64), allocatable, intent(inout) :: triplet_list(:, :)
      integer(int64), intent(out) :: triplet_count

      ! locals
      integer(int64) :: nTot, idx, i0, q, r, current_capacity
      integer(int64) :: iq1_ibz, iq1, iq2, iq3_minus, iq3_plus, s1, s2, s3
      integer(int64) :: ilambda1, ilambda2, ilambda3_minus, ilambda3_plus
      real(real64) :: q1f(3), q2f(3), en1, en2, en3_minus, en3_plus
      integer(int64) :: q1i(3), q2i(3), q3i(3), neg_q2_indvec(3), neg_iq2
      real(real64) :: delta_minus, delta_plus
      logical :: has_minus, has_plus

      ! Statistics counters (instead of masks)
      integer(int64) :: stat_minus, stat_plus, stat_both, stat_neither

      ! Shrink_int64_2d will want to reduce the array to (3, triplet_count).
      logical, parameter :: use_shrink = .true.

      triplet_count = 0_int64
      stat_minus = 0_int64
      stat_plus = 0_int64  
      stat_both = 0_int64
      stat_neither = 0_int64

      ! Total possible processes.
      nTot = nwv_irred*nwv*nb*nb*nb
      current_capacity = size(triplet_list, 2)

      do idx = 1_int64, nTot
         i0 = idx - 1_int64

         ! iq1_ibz in [1, nwv_irred]
         call int_div(i0, nwv_irred, q, r)
         iq1_ibz = r + 1_int64

         ! iq2 in [1, nwv]
         call int_div(q, nwv, q, r)
         iq2 = r + 1_int64

         ! s1 in [1, nb]
         call int_div(q, nb, q, r)
         s1 = r + 1_int64

         ! s3 in [1, nb]
         call int_div(q, nb, q, r)
         s3 = r + 1_int64

         ! s2 in [1, nb]
         call int_div(q, nb, q, r)
         s2 = r + 1_int64

         iq1 = indexlist_irred(iq1_ibz)

         ! Calculate q3 for momentum selection rule
         q1f = wavevecs(iq1, :)
         q2f = wavevecs(iq2, :)

         q1i = int(nint(q1f*real(wvmesh, real64)), int64)
         q2i = int(nint(q2f*real(wvmesh, real64)), int64)

         ! q3- = fold(q1 - q2)
         q3i = imod(q1i - q2i, wvmesh)
         iq3_minus = mux_vector_dev(q3i, wvmesh, 0_int64)

         ! q3+ = fold(q1 + q2)
         q3i = imod(q1i + q2i, wvmesh)
         iq3_plus = mux_vector_dev(q3i, wvmesh, 0_int64)

         ! index of -q2
         neg_q2_indvec = imod(-q2i, wvmesh)
         neg_iq2 = mux_vector_dev(neg_q2_indvec, wvmesh, 0_int64)

         ! Energies selection rules
         en1 = ens(iq1, s1)
         en2 = ens(iq2, s2)
         en3_minus = ens(iq3_minus, s3)
         en3_plus = ens(iq3_plus, s3)

         if(en1*en2*en3_minus == 0.0_real64 .and. en1*en2*en3_plus == 0.0_real64) then
            stat_neither = stat_neither + 1_int64
            cycle
         end if

         ! Energy-conservation deltas
         if(use_Tetrahedra) then
            delta_minus = delta_fn_tetra(en1 - en3_minus, iq2, s2, wvmesh, simplex_map, simplex_count, simplex_evals)
            delta_plus  = delta_fn_tetra(en3_plus - en1, neg_iq2, s2, wvmesh, simplex_map, simplex_count, simplex_evals)
         else
            delta_minus = delta_fn_triang(en1 - en3_minus, iq2, s2, wvmesh, simplex_map, simplex_count, simplex_evals)
            delta_plus  = delta_fn_triang(en3_plus - en1, neg_iq2, s2, wvmesh, simplex_map, simplex_count, simplex_evals)
         end if

         has_minus = (delta_minus > 0.0_real64)
         has_plus = (delta_plus > 0.0_real64)

         ! Count
         if(has_minus .and. has_plus) then
            stat_both = stat_both + 1_int64
            stat_minus = stat_minus + 1_int64
            stat_plus = stat_plus + 1_int64
         else if(has_minus) then
            stat_minus = stat_minus + 1_int64
         else if(has_plus) then
            stat_plus = stat_plus + 1_int64
         else
            stat_neither = stat_neither + 1_int64
         end if

         ! Convert to state indices
         ilambda1 = mux_state(nb, s1, iq1)
         ilambda2 = mux_state(nb, s2, iq2)
         ilambda3_minus = mux_state(nb, s3, iq3_minus)
         ilambda3_plus = mux_state(nb, s3, iq3_plus)

         ! Build triplet_list - minus process
         if(has_minus) then
            triplet_count = triplet_count + 1_int64
            if(triplet_count > current_capacity) then
               call expand_int64_2d(triplet_list, 1.3_real64)
               current_capacity = size(triplet_list, 2)
            end if
            triplet_list(1, triplet_count) = ilambda1
            triplet_list(2, triplet_count) = ilambda2
            triplet_list(3, triplet_count) = ilambda3_minus
         end if

         ! Build triplet_list - Plus process
         if(has_plus) then
            triplet_count = triplet_count + 1_int64
            if(triplet_count > current_capacity) then
               call expand_int64_2d(triplet_list, 1.3_real64)
               current_capacity = size(triplet_list, 2)
            end if
            triplet_list(1, triplet_count) = ilambda1
            triplet_list(2, triplet_count) = ilambda2
            triplet_list(3, triplet_count) = ilambda3_plus
         end if

         ! For debug suff: Progress indicator every 100M transitions
         if(mod(idx, 100000000_int64) == 0) then
            print *, "Progress:", idx, "/", nTot, &
                 "( triplet_count=", triplet_count, ", capacity=", current_capacity, ")"
         end if

      end do

      if(use_shrink) then
         call shrink_int64_2d(triplet_list, triplet_count)
      end if
      print *, "Memory required for S_list (3 int64 per process): approx", &
           real(3_int64*triplet_count*8_int64, real64) / 1.0e9_real64, "GB"
      print *, "Physical capacity still allocated (processes):", size(triplet_list, 2)

      !statistics
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
      print *, "Actual S_count:", triplet_count
      print *, "Discrepancy:", triplet_count - (stat_minus + stat_plus)
      print *, ""

    end subroutine build_tripletlist_host_dynamic_without_mask
  end subroutine compute_V2_on_gpu
end module gpu_3ph_driver
