program V3offload

#ifdef _OPENACC
  use openacc
#endif

  use precision, only: i64, r64
  use misc, only: print_message, subtitle, timer, exit_with_message, mux_vector, demux_state, &
       mux_state, twonorm, demux_vector
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
  real(r64) :: val !value

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

  !Calculate ph-ph vertex (cpu, original)
  Vm2_calculator => Vm2_3ph_reference
  call t_event%start_timer('reference V- on cpu')
  call calculate_3ph_interaction(ph, crys, num, V2, Vm2_calculator)
  call t_event%end_timer('reference V- on cpu')
  print*, 'value = ', twonorm(pack(V2, .true.))

  !Calculate ph-ph vertex (cpu, refactor)
  Vm2_calculator => Vm2_3ph_refactor
  call t_event%start_timer('refactored V- on cpu')
  call calculate_3ph_interaction(ph, crys, num, V2, Vm2_calculator)
  call t_event%end_timer('refactored V- on cpu')
  print*, 'value = ', twonorm(pack(V2, .true.))

  !Calculate ph-ph vertex (gpu, algo 1)
  call t_event%start_timer('V- on gpu, algo 1')
  call calculate_3ph_interaction_gpu(ph, crys, num, V2)
  call t_event%end_timer('V- on gpu, algo 1')
  print*, 'value = ', twonorm(pack(V2, .true.))

  !Calculate ph-ph vertex (gpu, algo 2)
  call t_event%start_timer('V- on gpu, algo 2')
  call calculate_3ph_interaction_gpu_algo2(ph, crys, num, V2)
  call t_event%end_timer('V- on gpu, algo 2')
  print*, 'value = ', twonorm(pack(V2, .true.))

!!$  !Calculate ph-ph vertex (gpu, algo 3)
!!$  call t_event%start_timer('V- on gpu, algo 3')
!!$  call calculate_3ph_interaction_gpu_algo3(ph, crys, num, V2)
!!$  call t_event%end_timer('V- on gpu, algo 3')
!!$  print*, 'value = ', twonorm(pack(V2, .true.))

contains

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
end program V3offload
