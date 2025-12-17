module gpu_interface
  !! This module provides a clean Fortran API that wraps the C-compatible interface
  !!
  !! handling all the C interoperability details internally.

  use iso_c_binding
  use iso_fortran_env, only: int64, real64

  implicit none

  interface
     subroutine compute_V2_on_gpu_c(nb, nwv_irred, nwv, ntrip, &
          evecs, Index_i, Index_j, Index_k, ifc3, &
          indexlist_irred, wavevecs, wvmesh, reclatt, Rj, Rk, &
          simplex_map, simplex_count, simplex_evals, ens, use_tetra, &
          V2_p, triplet_list_p, triplet_count, istat, &
          nincmax, nsimplex, nverts) bind(C, name = "compute_V2_on_gpu_c")

       import :: c_int64_t, c_double_complex, c_double, c_int, c_ptr

       integer(c_int64_t), value :: nb, nwv_irred, nwv, ntrip
       integer(c_int64_t), value :: nincmax, nsimplex, nverts
       complex(c_double_complex) :: evecs(nwv, nb, nb)
       integer(c_int64_t) :: Index_i(*), Index_j(*), Index_k(*)
       real(c_double) :: ifc3(3, 3, 3, *)
       integer(c_int64_t) :: indexlist_irred(*)
       real(c_double) :: wavevecs(nwv, 3)
       integer(c_int64_t) :: wvmesh(3)
       real(c_double) :: reclatt(3, 3)
       real(c_double) :: Rj(3, *), Rk(3, *)

       ! ESR / delta-function data
       integer(c_int64_t) :: simplex_map(2, nwv, *)
       integer(c_int64_t) :: simplex_count(*)
       real(c_double) :: simplex_evals(nsimplex, nb, nverts)
       real(c_double) :: ens(nwv, nb)
       integer(c_int), value :: use_tetra

       !outputs
       type(c_ptr) :: V2_p
       type(c_ptr) :: triplet_list_p
       integer(c_int64_t) :: triplet_count
       integer(c_int) :: istat !error status.
     end subroutine compute_V2_on_gpu_c
  end interface

contains

  subroutine compute_V2_on_gpu(nb, nwv_irred, nwv, ntrip, &
       evecs, Index_i, Index_j, Index_k, ifc3, &
       indexlist_irred, wavevecs, wvmesh, reclatt, Rj, Rk, &
       simplex_map, simplex_count, simplex_evals, ens, use_tetra, &
       V2, triplet_list_out, triplet_count_out, istat)
    !! Fortran interface to GPU three-phonon matrix element calculation
    !!
    !! Computes V^2 matrix elements for three-phonon processes that satisfy
    !! both momentum selection rules (MSR) and energy selection rules (ESR).
    !! The calculation proceeds in two stages:
    !! First, build list of valid triplets satisfying selection rules (on CPU)
    !! then ompute matrix elements for valid triplets (on GPU).
    !! Output Parameters: 
    !! triplet_count_out: number of valid processes found, typically 1-5% of total possible transitions.
    !! triplet_list_out: state index triplets [3, triplet_count_out]. ecach column (lambda1, lambda2, lambda3)
    !! where lambda = (k, band)

    integer(int64), intent(in) :: nb, nwv_irred, nwv, ntrip
    complex(real64), intent(in) :: evecs(nwv, nb, nb)
    integer(int64), intent(in) :: Index_i(ntrip), Index_j(ntrip), Index_k(ntrip)
    real(real64), intent(in) :: ifc3(3, 3, 3, ntrip)
    integer(int64), intent(in) :: indexlist_irred(nwv_irred)
    real(real64), intent(in) :: wavevecs(nwv, 3)
    integer(int64), intent(in) :: wvmesh(3)
    real(real64), intent(in) :: reclatt(3, 3)
    real(real64), intent(in) :: Rj(3, ntrip), Rk(3, ntrip)

    ! ESR / delta-function data
    integer(int64), intent(in) :: simplex_map(:, :, :)
    integer(int64), intent(in) :: simplex_count(:)
    real(real64), intent(in) :: simplex_evals(:, :, :)
    real(real64), intent(in) :: ens(nwv, nb)
    logical, intent(in) :: use_tetra

    ! outputs
    real(real64), allocatable, intent(out) :: V2(:)
    integer(int64), allocatable, intent(out) :: triplet_list_out(:, :)
    integer(int64), intent(out) :: triplet_count_out
    integer, intent(out) :: istat

    ! Local variables
    integer(c_int) :: use_tetra_c, istat_c
    integer(c_int64_t) :: nincmax, nsimplex, nverts, triplet_count_c
    type(c_ptr) :: V2_p, triplet_list_p
    real(c_double), pointer :: V2_temp(:)
    integer(c_int64_t), pointer :: triplet_list_temp(:, :)
    integer(int64) :: i, j

    ! Convert logical to C int
    use_tetra_c = merge(1_c_int, 0_c_int, use_tetra)

    ! Get dimensions from assumed-shape arrays
    nincmax = int(size(simplex_map, 3), c_int64_t)
    nsimplex = int(size(simplex_evals, 1), c_int64_t)
    nverts = int(size(simplex_evals, 3), c_int64_t)

    ! Call C wrapper
    call compute_V2_on_gpu_c(nb, nwv_irred, nwv, ntrip, &
         evecs, Index_i, Index_j, Index_k, ifc3, &
         indexlist_irred, wavevecs, wvmesh, reclatt, Rj, Rk, &
         simplex_map, simplex_count, simplex_evals, ens, use_tetra_c, &
         V2_p, triplet_list_p, triplet_count_c, istat_c, &
         nincmax, nsimplex, nverts)

    istat = int(istat_c, kind = kind(istat))
    triplet_count_out = triplet_count_c

    ! copy data from C pointers to fortran allocatable arrays
    if(c_associated(V2_p) .and. triplet_count_c > 0) then
       call c_f_pointer(V2_p, V2_temp, [triplet_count_c])
       allocate(V2(triplet_count_c))
       ! Copy element by element to ensure data is captured
       do i = 1, triplet_count_c
          V2(i) = V2_temp(i)
       end do
       ! Disassociates a pointer from a target.
       nullify(V2_temp)  ! I don't deallocate it cz memory owned by wrapper
    else
       allocate(V2(0))
    end if

    if(c_associated(triplet_list_p) .and. triplet_count_c > 0) then
       call c_f_pointer(triplet_list_p, triplet_list_temp, [3_int64, triplet_count_c])
       allocate(triplet_list_out(3, triplet_count_c))
       ! Copy element by element
       do j = 1, triplet_count_c
          do i = 1, 3
             triplet_list_out(i, j) = triplet_list_temp(i, j)
          end do
       end do
       nullify(triplet_list_temp)  ! I don't deallocate cz memory owned by wrapper
    else
       ! No valid triplets found, then return empty array
       allocate(triplet_list_out(3, 0))
    end if

  end subroutine compute_V2_on_gpu
end module gpu_interface

