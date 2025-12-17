module gpu_3ph_kernels
  use cudafor
  use iso_fortran_env, only: int64, real64

  implicit none

contains

  attributes(device) real(real64) function Vm2_3ph_reference_dev( &
       ev1_s1, ev2_s2, ev3_s3, Index_i, Index_j, Index_k, ifc3, &
       q2c, q3c, Rj, Rk, ntrip, nb)
    !! Function to calculate the squared 3-ph interaction vertex |V-|^2.

    integer(int64), value :: ntrip, nb
    integer(int64) :: Index_i(ntrip), Index_j(ntrip), Index_k(ntrip)
    complex(real64) :: ev1_s1(nb), ev2_s2(nb), ev3_s3(nb)
    real(real64) :: ifc3(3, 3, 3, ntrip), q2c(3), q3c(3), Rj(3,ntrip), Rk(3,ntrip)

    integer(int64) :: it, a, b, c, aind, bind, cind
    complex(real64) :: aux1, aux2, aux3, V0, phase
    real(real64) :: arg

    aux1 = (0.0_real64, 0.0_real64)
    do it = 1, ntrip
       aind = 3*(Index_k(it) - 1)
       bind = 3*(Index_j(it) - 1)
       cind = 3*(Index_i(it) - 1)
       V0 = (0.0_real64, 0.0_real64)
       do a = 1, 3
          aux2 = conjg(ev3_s3(a + aind))
          do b = 1, 3
             aux3 = aux2*conjg(ev2_s2(b + bind))
             do c = 1, 3
                if(ifc3(c, b, a, it) /= 0.0_real64) then
                   V0 = V0 + ifc3(c, b, a, it)*ev1_s1(c + cind)*aux3
                end if
             end do
          end do
       end do
       arg = dot_product(q2c, Rj(:, it)) + dot_product(q3c, Rk(:, it))
       phase = exp((0.0_real64, -1.0_real64)*arg)
       aux1 = aux1 + V0*phase
    end do
    Vm2_3ph_reference_dev = abs(aux1)**2
  end function Vm2_3ph_reference_dev

  pure real(real64) function twonorm_real_rank1(v)
    !! 2-norm of a rank-1 real vector

    real(real64), intent(in) :: v(:)
    integer(int64) :: i, s

    s = size(v)
    twonorm_real_rank1 = 0.0_real64
    do i = 1, s
       twonorm_real_rank1 = v(i)**2 + twonorm_real_rank1
    end do
    twonorm_real_rank1 = sqrt(twonorm_real_rank1)
  end function twonorm_real_rank1

  attributes(host, device) pure integer(int64) function mux_state(nbands, iband, ik)
    !! Multiplex a (band index, wave vector index) pair into a state index
    !!
    !! nbands is the number of bands
    !! iband is the band index
    !! ik is the wave vector index

    integer(int64), value :: nbands, ik, iband

    mux_state = (ik - 1_int64)*nbands + iband
  end function mux_state

  attributes(host, device) pure subroutine demux_state(nbands, istate, iband, ik)
    !! Demultiplex a state index into (band index, wave vector index) pair
    !!
    !! nbands is the number of bands
    !! istate is the multiplexed state index
    !! iband is the band index
    !! ik is the wave vector index

    integer(int64), value :: nbands, istate
    integer(int64), intent(out) :: iband, ik

    iband = imod(istate - 1_int64, nbands) + 1_int64
    ik = (istate - 1_int64) / nbands + 1_int64
  end subroutine demux_state

  attributes(host, device) subroutine int_div(num, denom, q, r)
    !! Quotient(q) and remainder(r) of the integer division num/denom.

    integer(int64), value :: num, denom
    integer(int64) :: q, r   !can make outputs explicit: ,intent(out)

    q = num / denom
    r = mod(num, denom)
  end subroutine int_div

  attributes(host, device) elemental integer(int64) function imod(a, m)
    !! return only the positive remainder for 64-bit integers.
    !!
    !! Added the elemental to make it work on arrays as mod or modulo.

    integer(int64), value :: a, m

    imod = mod(a, m)
    if(imod < 0_int64) imod = imod + m
  end function imod

  attributes(host, device) integer(int64) function mux_vector_dev(v, mesh, base)
    !! Multiplex index of a single wave vector.
    !! Output is always 1-based.
    !! v is the demultiplexed triplet of a wave vector.
    !! mesh is the number of wave vectors along the three reciprocal lattice vectors.
    !! base states whether v has 0- or 1-based indexing.

    integer(int64), value :: base
    integer(int64), intent(in) :: v(3), mesh(3)

    if(base == 0_int64) then
       ! v is 0-based; return 1-based linear index
       mux_vector_dev = (v(3)*mesh(2) + v(2))*mesh(1) + v(1) + 1_int64
    else
       ! v is 1-based; return 1-based linear index
       mux_vector_dev = ((v(3) - 1_int64)*mesh(2) + (v(2) - 1_int64))*mesh(1) + v(1)
    end if
  end function mux_vector_dev

  attributes(host) pure real(real64) function delta_fn_dispatch( &
       e, ik, ib, mesh, map3, count, evals, use_tetra) result(val)
    !! Return pointer to either tetrahedra or tringular delta
    !! function evaulator.

    real(real64), intent(in) :: e
    integer(int64), intent(in) :: ik, ib
    integer(int64), intent(in) :: mesh(3)
    integer(int64), intent(in) :: map3(:, :, :)
    integer(int64), intent(in) :: count(:)
    real(real64), intent(in) :: evals(:, :, :)
    logical, intent(in) :: use_tetra

    if(use_tetra) then
       ! map3/count/evals are the tetra variants here
       val = delta_fn_tetra(e, ik, ib, mesh, map3, count, evals)
    else
       ! map3/count/evals are the triangle variants here
       val = delta_fn_triang(e, ik, ib, mesh, map3, count, evals)
    end if
  end function delta_fn_dispatch

  attributes(host) pure real(real64) function safe_div(a, b) result(q)
    !! Returns a/b when abs(b) is sufficiently large; otherwise returns 0.0.

    real(real64), intent(in) :: a, b
    real(real64), parameter :: eps = 1.0e-14_real64
    if (abs(b) > eps) then
       q = a / b
    else
       q = 0.0_real64
    end if
  end function safe_div

  attributes(host) pure real(real64) function delta_fn_triang(e, ik, ib, mesh, triangmap, triangcount, triang_evals)
    !! Calculate delta function using the triangle method a la
    !! Kurganskii et al. Phys. Stat. Sol.(b) 129, 293 (1985)
    !!
    !! e Sample energy
    !! ik Wave vector index
    !! ib Band index
    !! mesh Wave vector grid
    !! triangmap Wave vector to (triangle, vertex) mapping
    !! triangcount Number of triangles in which a wave vector belongs
    !! triang_evals Triangles populated with the eigenvalues

    real(real64), intent(in) :: e
    integer(int64), intent(in) :: ik, ib
    integer(int64), intent(in) :: mesh(3), triangmap(:,:,:), triangcount(:)
    real(real64), intent(in) :: triang_evals(:,:,:)

    !Local variables
    integer(int64) :: iv, it, itk, num, numtriangs
    logical :: c1, c2, c3, c4
    real(real64) :: e1, e2, e3, E12, E21, E13, E31, E23, E32, tmp

    tmp = 0.0_real64
    delta_fn_triang = 0.0_real64

    !Total number of triangles in the system
    numtriangs = product(mesh)*2

    !Grab number of triangles in which wave vector belongs
    num = triangcount(ik)

    do itk = 1, num !Run over triangles
       it = triangmap(1, ik, itk) !Grab triangle
       iv = triangmap(2, ik, itk) !Grab vertex

       !Grab vertex energies
       e1 = triang_evals(it, ib, 1)
       e2 = triang_evals(it, ib, 2)
       e3 = triang_evals(it, ib, 3)

       !Evaluate the four possible cases
       c1 = e <= e1
       c2 = e1 < e .and. e <= e2
       c3 = e2 < e .and. e <= e3
       c4 = e3 < e

       tmp = 0.0_real64

       if(c1 .or. c4) cycle

       !Define Eij
       ! Note that at this stage the quantities below might
       ! be ill defined due to degeneracies. But the conditionals
       ! that will follow will take this into account.
       E12 = (e - e2)/(e1 - e2)
       E21 = (e - e1)/(e2 - e1)
       E13 = (e - e3)/(e1 - e3)
       E31 = (e - e1)/(e3 - e1)
       E23 = (e - e3)/(e2 - e3)
       E32 = (e - e2)/(e3 - e2)

       select case(iv)
       case(1)
          if(c2) then
             tmp = E21*(E12 + E13)/(e3 - e1)
          else if(c3) then
             tmp = E23*E13/(e3 - e1)
          end if
       case(2)
          if(c2) then
             tmp = E21*E21/(e3 - e1)
          else if(c3) then
             tmp = E23*E23/(e3 - e1)
          end if
       case(3)
          if(c2) then
             tmp = E21*E31/(e3 - e1)
          else if(c3) then
             tmp = E23*(E31 + E32)/(e3 - e1)
          end if
       end select

       delta_fn_triang = delta_fn_triang + tmp
    end do !itk

    if(delta_fn_triang < 1.0e-12_real64) delta_fn_triang = 0.0_real64

    !Normalize with the total number of triangles
    delta_fn_triang = delta_fn_triang/numtriangs
  end function delta_fn_triang

  attributes(host) pure real(real64) function delta_fn_tetra(e, ik, ib, mesh, tetramap, tetracount, tetra_evals)
    !! Calculate delta function using the tetraheron method.
    !!
    !! e Sample energy
    !! ik Wave vector index
    !! ib Band index
    !! mesh Wave vector grid
    !! tetramap Wave vector to (tetrahedron, vertex) mapping
    !! tetracount Number of tetrahedra in which a wave vector belongs
    !! tetra_evals Tetrahedra populated with the eigenvalues

    real(real64), intent(in) :: e
    integer(int64), intent(in) :: ik, ib
    integer(int64), intent(in) :: mesh(3), tetramap(:,:,:), tetracount(:)
    real(real64), intent(in) :: tetra_evals(:,:,:)

    !Local variables
    integer(int64) :: iv, it, itk, num, numtetra
    logical :: c1, c2, c3
    real(real64) :: e1, e2, e3, e4, e1e, e2e, e3e, e4e, &
         e21, e31, e41, e32, e42, e43, tmp ! eji \equiv ej - ei

    tmp = 0.0_real64
    delta_fn_tetra = 0.0_real64

    !Total number of tetrahedra in the system
    numtetra = product(mesh)*6

    !Grab number of tetrahedra in which wave vector belongs
    num = tetracount(ik)

    do itk = 1, num !Run over tetrahedra
       it = tetramap(1, ik, itk) !Grab tetrahedron
       iv = tetramap(2, ik, itk) !Grab vertex

       !Grab vertex energies
       e1 = tetra_evals(it, ib, 1)
       e2 = tetra_evals(it, ib, 2)
       e3 = tetra_evals(it, ib, 3)
       e4 = tetra_evals(it, ib, 4)

       !Define the energy differences
       e1e = e1 - e
       e2e = e2 - e
       e3e = e3 - e
       e4e = e4 - e
       e21 = e2 - e1
       e31 = e3 - e1
       e41 = e4 - e1
       e32 = e3 - e2
       e42 = e4 - e2
       e43 = e4 - e3

       !Evaluate the three cases
       c1 = e1 <= e .and. e <= e2
       c2 = e2 <= e .and. e <= e3
       c3 = e3 <= e .and. e <= e4

       if(.not. (e < e1 .or. e > e4)) then
          !Evaluate the expressions for the three cases
          select case(iv)
          case(1)
             if(c1) then
                tmp = (e2e/e21 + e3e/e31 + e4e/e41)*(e1e**2)/e41/e31/e21

                if(e1 == e2) then
                   tmp = 0.0_real64
                end if
             else if(c2) then
                tmp = -0.5_real64*(e3e/(e31**2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41) &
                     + e4e/(e41**2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32))

                if(e2 == e3) then
                   tmp = -0.5_real64*(e4e*e1e/e41/e42 + e1e/e41 &
                        + e4e/(e41**2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31))
                end if
             else if(c3) then
                tmp = (e4e**3)/(e41**2)/e42/e43

                if(e3 == e4) then
                   tmp = (e4e**2)/(e41**2)/e42
                end if
             end if
          case(2)
             if(c1) then
                tmp = -(e1e**3)/(e21**2)/e31/e41

                if(e1 == e2) then
                   tmp = 0.0_real64
                end if
             else if(c2) then
                tmp = -0.5_real64*(e3e/(e32**2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41) &
                     + e4e/(e42**2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41))

                if(e2 == e3) then
                   tmp = -0.5_real64*(0.0 + e4e/e42/e41 + 0.0 &
                        + e4e/(e42**2)*(0.0 + e4e*e1e/e41/e31 + 1.0))
                end if
             else if(c3) then
                tmp = (e4e**3)/e41/(e42**2)/e43

                if(e3 == e4) then
                   tmp = 0.0_real64
                end if
             end if
          case(3)
             if(c1) then
                tmp = -(e1e**3)/e21/(e31**2)/e41

                if(e1 == e2) then
                   tmp = 0.0_real64
                end if
             else if(c2) then
                tmp = 0.5_real64*(e2e/(e32**2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41) &
                     + e1e/(e31**2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41))

                if(e2 == e3) then
                   tmp = 0.5_real64*(0.0 + e4e/e42/e41 + e1e/e31/e41 &
                        + e1e/(e31**2)*(0.0 + e4e*e1e/e41/e42 + e1e/e41))
                end if
             else if(c3) then
                tmp = (e4e**3)/e41/e42/(e43**2)

                if(e3 == e4) then
                   tmp = 0.0_real64
                end if
             end if
          case(4)
             if(c1) then
                tmp = -(e1e**3)/e21/e31/(e41**2)
                if(e1 == e2) then
                   tmp = 0.0_real64
                end if
             else if(c2) then
                tmp = 0.5_real64*(e2e/(e42**2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41) &
                     + e1e/(e41**2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32))

                if(e2 == e3) then
                   tmp = 0.5_real64*(0.0 &
                        + e1e/(e41**2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31))
                end if
             else if(c3) then
                tmp = -(e3e/e43 + e2e/e42 + e1e/e41)*(e4e**2)/e41/e42/e43

                if(e3 == e4) then
                   tmp = 0.0_real64
                end if
             end if
          end select

          if ((e1 == e2) .and. (e1 == e3) .and. (e1 == e4) .and. (e == e1)) then
             tmp = 0.25_real64
          end if

          delta_fn_tetra = delta_fn_tetra + tmp
       end if ! .not. (e <= e1 .or. e >= e4)
    end do !itk

    if(delta_fn_tetra < 1.0e-12_real64) delta_fn_tetra = 0.0_real64

    !Normalize with the total number of tetrahedra
    delta_fn_tetra = delta_fn_tetra/numtetra
  end function delta_fn_tetra

  attributes(global) subroutine calculate_Vm2_kernel(V2, evecs, triplet_list, &
       Index_i, Index_j, Index_k, ifc3, indexlist_irred, wavevecs, &
       wvmesh, reclatt, Rj, Rk, triplet_count, nb, nwv, ntrip)
    !! Kernel: Calculate Vm2 for each entry in triplet_list and store directly in 1D V2 array.
    !!
    !! V2 is now a 1D array with size triplet_count

    integer(int64), value :: triplet_count, nb, nwv, ntrip

    !1D output array
    real(real64) :: V2(triplet_count)

    ! Input arrays
    complex(real64) :: evecs(nwv, nb, nb)
    integer(int64) :: triplet_list(3, triplet_count) ! (ilambda1, ilambda2, ilambda3)
    integer(int64) :: Index_i(ntrip), Index_j(ntrip), Index_k(ntrip)
    real(real64) :: ifc3(3, 3, 3, ntrip)
    integer(int64) :: indexlist_irred(*)
    real(real64) :: wavevecs(nwv, 3)
    integer(int64) :: wvmesh(3)
    real(real64) :: reclatt(3, 3)
    real(real64) :: Rj(3, ntrip), Rk(3, ntrip)

    integer(int64) :: idx, iq1_ibz, iq1, iq2, iq3, s1, s2, s3, i
    integer(int64) :: ilambda1, ilambda2, ilambda3
    complex(real64) :: ev1(nb), ev2(nb), ev3(nb)
    real(real64) :: q1f(3), q2f(3), q3f(3), q2c(3), q3c(3)
    integer(int64) :: q1i(3), q2i(3), q3i(3)

    ! This line is CUDA Fortran syntax, used to assign a unique thread ID in a GPU kernel.
    idx = (int(blockIdx%x, int64) - 1_int64)*int(blockDim%x, int64) + int(threadIdx%x, int64)
    if(idx < 1_int64 .or. idx > triplet_count) return

    ! Read compressed indices
    ilambda1 = triplet_list(1, idx)
    ilambda2 = triplet_list(2, idx)
    ilambda3 = triplet_list(3, idx)

    ! Demultiplex back to (wavevector, band) pairs
    call demux_state(nb, ilambda1, s1, iq1)
    call demux_state(nb, ilambda2, s2, iq2)
    call demux_state(nb, ilambda3, s3, iq3)

    ! Calculate q vectors
    ! Initial (IBZ blocks) wave vector (crystal coords.)
    q1f = wavevecs(iq1, :)
    q2f = wavevecs(iq2, :)
    q3f = wavevecs(iq3, :)

    ! Convert to Cartesian coordinates
    q2c = matmul(reclatt, q2f)
    q3c = matmul(reclatt, q3f) 

    ! Gather eigenvectors
    do i = 1_int64, nb
       ev1(i) = evecs(iq1, s1, i)
       ev2(i) = evecs(iq2, s2, i)
       ev3(i) = evecs(iq3, s3, i)
    end do

    ! Calculate matrix element and store directly in 1D V2 array
    V2(idx) = Vm2_3ph_reference_dev(ev1, ev2, ev3, Index_i, Index_j, Index_k, &
         ifc3, q2c, q3c, Rj, Rk, ntrip, nb)
  end subroutine calculate_Vm2_kernel
end module gpu_3ph_kernels
