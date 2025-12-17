module gpu_wrapper_c
  !! C-interface wrapper for GPU three-phonon interaction calculations
  !!
  !! This module provides a C-compatible interface to the Fortran GPU driver,
  !! managing memory storage for data that must outlive individual function calls.

  use iso_c_binding
  use iso_fortran_env, only: int64, real64
  use gpu_3ph_driver, only: compute_V2_on_gpu

  implicit none

  ! Where we keep the computed V2 values.
  real(real64), allocatable, save :: V2_storage(:)
  ! Where we keep the list of valid processes.
  integer(int64), allocatable, save :: triplet_list_storage(:, :)

contains

  subroutine compute_V2_on_gpu_c(nb, nwv_irred, nwv, ntrip, &
       evecs, Index_i, Index_j, Index_k, ifc3, &
       indexlist_irred, wavevecs, wvmesh, reclatt, Rj, Rk, &
       simplex_map, simplex_count, simplex_evals, ens, use_tetra, &
       V2_p, triplet_list_p, triplet_count, istat, &
       nincmax, nsimplex, nverts) bind(C, name = "compute_V2_on_gpu_c")
    !! This subroutine computes V^2 matrix elements for valid three-phonon processes
    !!
    !! that satisfy both momentum and energy conservation rules.
    !! Inputs dimensions: nb: number of phonon bands, nwv_irred: number of irreducible k-points,
    !! nwv: total number of k-points in full grid, ntrip: number of triplets, 
    !! nincmax: max simplicies any k-points belongs to, nsimplex: total number of simplicies (triangles/tetrahedra)
    !! nverts: number of vertices/simplex; 3 for triangles and 4 for tetrahedra.
    !! input arrays passed from C:
    !! evecs: phonon eigenvectors [nwv, nb, nb].
    !! Index_i, Index_j, Index_k: atom indices for force constants [ntrip].
    !! ifc3: third-order force constant [3, 3, 3, ntrip].
    !! indexlist_irred: mapping to apping to irreducible k-points [wv_irred].
    !! wavevecs: K-point coordinates [nwv, 3], nvmesh: K-point grid dimensions [3].
    !! reclatt: Reciprocal lattice vectors [3, 3], Rj, Rk: Lattice vectors for phases [3, ntrip].
    !! simplex_map: K-point to simplex mapping [2, nwv, nincmax], simplex_count: number of simplices per k-point [nwv].
    !! simplex_evals: energies at simplex vertices [nsimplex, nb, nverts].
    !! ens: phonon energies [nwv, nb], use_tetra: 1=tetrahedron method, 0=triangle method.
    !! Output (returned via pointers):
    !! V2_p: pointer to the matrix element V2 values [ntrip], triplet_list_p: pointer to state index triplets[3, ntrip].
    !! triplet_count: number of valid triplet satisfying selection rules. istat: 0=success, negative=error code.

    integer(c_int64_t), value :: nb, nwv_irred, nwv, ntrip
    integer(c_int64_t), value :: nincmax, nsimplex, nverts

    complex(c_double_complex) :: evecs(nwv, nb, nb)
    integer(c_int64_t) :: Index_i(ntrip), Index_j(ntrip), Index_k(ntrip)
    real(c_double) :: ifc3(3, 3, 3, ntrip)
    integer(c_int64_t) :: indexlist_irred(nwv_irred)
    real(c_double) :: wavevecs(nwv, 3)
    integer(c_int64_t) :: wvmesh(3)
    real(c_double) :: reclatt(3, 3)
    real(c_double) :: Rj(3, ntrip), Rk(3, ntrip)

    ! ESR / delta-function data
    integer(c_int64_t) :: simplex_map(2, nwv, nincmax)
    integer(c_int64_t) :: simplex_count(nwv)
    real(c_double) :: simplex_evals(nsimplex, nb, nverts)
    real(c_double) :: ens(nwv, nb)
    integer(c_int), value :: use_tetra

    ! Output pointers to 1D V2 and 2D S_list
    type(c_ptr) :: V2_p  ! Will point to V2_storage(valid_triplet_count)
    type(c_ptr) :: triplet_list_p ! Will point to valid triplet(3, valid_triplet_count)
    integer(c_int64_t) :: triplet_count ! Number of valid transitions

    integer(c_int) :: istat

    ! Local variables
    logical :: use_tetra_fortran
    integer :: istat_fortran
    integer(int64) :: valid_triplet_count

    use_tetra_fortran = (use_tetra /= 0)

    ! Deallocate previous persistent arrays if they exist
    if(allocated(V2_storage)) deallocate(V2_storage)
    if(allocated(triplet_list_storage)) deallocate(triplet_list_storage)

    call compute_V2_on_gpu(nb, nwv_irred, nwv, ntrip, &
         evecs, Index_i, Index_j, Index_k, ifc3, &
         indexlist_irred, wavevecs, wvmesh, reclatt, Rj, Rk, &
         simplex_map, simplex_count, simplex_evals, ens, use_tetra_fortran, &
         V2_storage, triplet_list_storage, valid_triplet_count, istat_fortran)

    istat = int(istat_fortran, kind = c_int)
    triplet_count = valid_triplet_count

    ! These pointers remain valid until cleanup_gpu_arrays() is called
    if(allocated(V2_storage) .and. size(V2_storage) > 0) then
       V2_p = c_loc(V2_storage(1))
    else
       V2_p = c_null_ptr
    end if

    ! triplet_list_storage is 2D: (3, triplet_count)
    if(allocated(triplet_list_storage) .and. size(triplet_list_storage, 2) > 0) then
       ! Return pointer to first element of 2D array
       triplet_list_p = c_loc(triplet_list_storage(1, 1))
    else
       triplet_list_p = c_null_ptr
    end if

  end subroutine compute_V2_on_gpu_c

  subroutine cleanup_gpu_arrays() bind(C, name = "cleanup_gpu_arrays")
    !! Free memory allocated for GPU computation results.
    !!
    !! After calling this, any pointers returned by compute_V2_on_gpu_c become invalid.

    if(allocated(V2_storage)) deallocate(V2_storage)
    if(allocated(triplet_list_storage)) deallocate(triplet_list_storage)
  end subroutine cleanup_gpu_arrays

end module gpu_wrapper_c
