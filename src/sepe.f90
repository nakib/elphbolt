! Copyright 2024 elphbolt contributors.
! This file is part of elphbolt <https://github.com/nakib/elphbolt>.
!
! elphbolt is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! elphbolt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with elphbolt. If not, see <http://www.gnu.org/licenses/>.

module SEPE_module
  !! Module containing type and procedures related to the solution of the
  !! Semiconductor Electron-Phonon Equations (SEPE) a la Stefanucci & Perfetto.

  use precision, only: r64, i64
  use params, only: qe, kB, hbar_eVps, oneI, complex_zero
!!$  use misc, only: Bose, mux_vector, binsearch, timer, subtitle, &
!!$       print_message, distribute_points, demux_state, trace
  use misc, only: print_message, exit_with_message, write2file_rank2_real, &
       distribute_points, demux_state, binsearch, interpolate, demux_vector, mux_vector, &
       trace, subtitle, append2file_transport_tensor, write2file_response, &
       linspace, readfile_response, write2file_spectral_tensor, subtitle, timer, &
       twonorm, write2file_rank1_real, precompute_interpolation_corners_and_weights, &
       interpolate_using_precomputed, Jacobian, cross_product, qdist, Bose
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use phonon_module, only: phonon
  use electron_module, only: electron
  use interactions, only: calculate_ph_rta_rates, read_transition_probs_e, &
       calculate_el_rta_rates, calculate_bound_scatt_rates, calculate_thinfilm_scatt_rates, &
       calculate_4ph_rta_rates, calculate_W3ph_OTF, calculate_Y_OTF, calculate_Xee_OTF, &
       calculate_Xee_13_OTF, calculate_ph_rta_coherence_rates
  use bz_sums, only: calculate_transport_coeff, calculate_spectral_transport_coeff, &
       calculate_cumulative_transport_coeff

  implicit none

  private
  public sepe

  type sepe
     !! Data and procedures related to the BTE.

     real(r64), allocatable :: ph_rta_rates_iso_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to isotope scattering.
     real(r64), allocatable :: ph_rta_rates_subs_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to substitution scattering.
     real(r64), allocatable :: ph_rta_rates_bound_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to boundary scattering.
     real(r64), allocatable :: ph_rta_rates_thinfilm_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to thin-film scattering.
     real(r64), allocatable :: ph_rta_rates_3ph_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to 3-ph interactions.
     real(r64), allocatable :: ph_rta_rates_4ph_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to 4-ph interactions.
     real(r64), allocatable :: ph_rta_rates_phe_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to ph-e interactions.
     real(r64), allocatable :: ph_rta_rates_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ.
     real(r64), allocatable :: ph_rta_coherence_rates_ibz(:,:)
     !! Phonon RTA coherence rates on the IBZ.
     real(r64), allocatable :: ph_field_term_T(:,:,:)
     !! Phonon field coupling term for gradT field on the FBZ.
     real(r64), allocatable :: ph_response_T(:,:,:)
     !! Phonon response function for gradT field on the FBZ.
     real(r64), allocatable :: ph_field_term_E(:,:,:)
     !! Phonon field coupling term for E field on the FBZ.
     real(r64), allocatable :: ph_response_E(:,:,:)
     !! Phonon response function for E field on the FBZ.
     complex(r64), allocatable :: ph_coherence_T(:, :, :)
     !! Phonon coherence term for gradT field on the FBZ
     complex(r64), allocatable :: ph_coherence_E(:, :, :)
     !! Phonon coherence term for E field on the FBZ

     real(r64), allocatable :: el_rta_rates_echimp_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to charged impurity scattering.
     real(r64), allocatable :: el_rta_rates_bound_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to boundary scattering.
     real(r64), allocatable :: el_rta_rates_eph_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to e-ph interactions.
     real(r64), allocatable :: el_rta_rates_ee_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to e-e interactions.
     real(r64), allocatable :: el_rta_rates_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ.
     real(r64), allocatable :: el_field_term_T(:,:,:)
     !! Electron field coupling term for gradT field on the FBZ.
     real(r64), allocatable :: el_response_T(:,:,:)
     !! Electron response function for gradT field on the FBZ.
     real(r64), allocatable :: el_field_term_E(:,:,:)
     !! Electron field coupling term for E field on the FBZ.
     real(r64), allocatable :: el_response_E(:,:,:)
     !! Electron response function for E field on the FBZ.
   contains

     procedure :: solve_sepe=>sepe_driver

  end type sepe

  type transport_coeffs
     !! Module level private data pack for all the transport coefficients.

     real(r64), allocatable :: ph_kappa(:,:,:), ph_alphabyT(:,:,:), &
          dummy(:,:,:), I_diff(:,:,:), I_drag(:,:,:), el_kappa0(:,:,:), el_alphabyT(:,:,:), &
          el_sigma(:, :,:), el_sigmaS(:, :, :), ph_drag_term_T(:,:,:), ph_drag_term_E(:,:,:)

     real(r64) :: ph_kappa_scalar, ph_kappa_scalar_old, ph_alphabyT_scalar, ph_alphabyT_scalar_old, &
          el_kappa0_scalar, el_kappa0_scalar_old, el_alphabyT_scalar, el_alphabyT_scalar_old, &
          el_sigma_scalar, el_sigma_scalar_old, el_sigmaS_scalar, el_sigmaS_scalar_old, KO_dev, lambda, &
          tot_alphabyT_scalar

   contains

     procedure :: initialize_ph=>allocate_ph_transport_coeffs, &
          initialize_el=>allocate_el_transport_coeffs

  end type transport_coeffs

contains

  subroutine allocate_ph_transport_coeffs(self, ph_numbands)
    !! Allocator of the phonon transport coefficients in the pack.

    class(transport_coeffs), intent(inout) :: self
    integer(i64), intent(in) :: ph_numbands

    allocate(self%ph_kappa(ph_numbands, 3, 3), self%ph_alphabyT(ph_numbands, 3, 3), &
         self%dummy(ph_numbands, 3, 3))
  end subroutine allocate_ph_transport_coeffs

  subroutine allocate_el_transport_coeffs(self, el_numbands)
    !! Allocator of the electron transport coefficients in the pack.

    class(transport_coeffs), intent(inout) :: self
    integer(i64), intent(in) :: el_numbands

    allocate(self%el_sigma(el_numbands, 3, 3), self%el_sigmaS(el_numbands, 3, 3), &
         self%el_alphabyT(el_numbands, 3, 3), self%el_kappa0(el_numbands, 3, 3))
  end subroutine allocate_el_transport_coeffs

  subroutine sepe_driver(self, num, crys, sym, ph, el)
    !! Subroutine to orchestrate the SEPE calculations.
    !!
    !! self SEPE object
    !! num Numerics object
    !! crys Crystal object
    !! sym Symmertry object
    !! ph Phonon object
    !! el Electron object

    class(sepe), intent(inout) :: self
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(phonon), intent(in) :: ph
    type(electron), intent(in), optional :: el

    call subtitle("Calculating SEPEs...")

    !Phonon RTA
    if(.not. num%onlyebte) &
         call dragless_phbte_RTA(num%cwd_T, self, num, crys, sym, ph, el)

    !Electron RTA
    if(.not. num%onlyphbte) &
         call dragless_ebte_RTA(num%cwd_T, self, num, crys, sym, el, ph)

    !Dragful electron-phonon BTEs
    if(num%drag) &
         call dragfull_ephbtes(num%cwd_T, self, num, crys, sym, ph, el)
  end subroutine sepe_driver

  subroutine calculate_field_term(species, field, nequiv, ibz2fbz_map, &
       T, chempot, ens, vels, rta_rates_ibz, field_term, el_indexlist)
    !! Subroutine to calculate the field coupling term of the BTE.
    !!
    !! species Type of particle
    !! field Type of field
    !! nequiv List of the number of equivalent points for each IBZ wave vector
    !! ibz2fbz_map Map from an FBZ wave vectors to its IBZ wedge image
    !! T Temperature in K
    !! ens FBZ energies
    !! vels FBZ velocities
    !! chempot Chemical potential (should be 0 for phonons)
    !! rta_rates_ibz IBZ RTA scattering rates
    !! field_term FBZ field-coupling term of the BTE
    !! el_indexlist [Optional] 

    character(len = 2), intent(in) :: species
    character(len = 1), intent(in) :: field
    integer(i64), intent(in) :: nequiv(:), ibz2fbz_map(:,:,:)
    real(r64), intent(in) :: T, chempot, ens(:,:), vels(:,:,:), rta_rates_ibz(:,:)
    real(r64), allocatable, intent(out) :: field_term(:,:,:)
    integer(i64), intent(in), optional :: el_indexlist(:)

    !Local variables
    integer(i64) :: ik_ibz, ik_fbz, ieq, ib, nk_ibz, nk, nbands, pow, &
         chunk, num_active_images, start, end
    real(r64) :: A
    logical :: trivial_case

    !Set constant and power of energy depending on species and field type
    if(species == 'ph') then
       A = 1.0_r64/T
       pow = 1
       if(chempot /= 0.0_r64) then
          call exit_with_message("Phonon chemical potential non-zero in calculate_field_term. Exiting.")
       end if
    else if(species == 'el') then
       if(field == 'T') then
          A = 1.0_r64/T
          pow = 1
       else if(field == 'E') then
          A = qe
          pow = 0
       else
          call exit_with_message("Unknown field type in calculate_field_term. Exiting.")
       end if
    else
       call exit_with_message("Unknown particle species in calculate_field_term. Exiting.")
    end if

    !Number of IBZ wave vectors
    nk_ibz = size(rta_rates_ibz(:,1))

    !Number of FBZ wave vectors
    nk = size(ens(:,1))

    !Number of bands
    nbands = size(ens(1,:))

    !Allocate and initialize field term
    allocate(field_term(nk, nbands, 3))
    field_term(:,:,:) = 0.0_r64

    !No field-coupling case
    trivial_case = species == 'ph' .and. field == 'E'

    if(.not. trivial_case) then
       !Divide IBZ states among images
       call distribute_points(nk_ibz, chunk, start, end, num_active_images)

       !Only work with the active images
       if(this_image() <= num_active_images) then
          do ik_ibz = start, end
             do ieq = 1, nequiv(ik_ibz)
                if(species == 'ph') then
                   ik_fbz = ibz2fbz_map(ieq, ik_ibz, 2)
                else
                   !Find index of electron in indexlist
                   call binsearch(el_indexlist, ibz2fbz_map(ieq, ik_ibz, 2), ik_fbz)
                end if
                do ib = 1, nbands
                   if(rta_rates_ibz(ik_ibz, ib) /= 0.0_r64) then
                      field_term(ik_fbz, ib, :) = A*vels(ik_fbz, ib, :)*&
                           (ens(ik_fbz, ib) - chempot)**pow/rta_rates_ibz(ik_ibz, ib)
                   end if
                end do
             end do
          end do
       end if

       !Reduce field term
       !Units:
       ! nm.eV/K for phonons, gradT-field
       ! nm.eV/K for electrons, gradT-field
       ! nm.C for electrons, E-field
       call co_sum(field_term)
    end if
  end subroutine calculate_field_term

  subroutine dragless_phbte_RTA(Tdir, self, num, crys, sym, ph, el)
    !! Dragless phonon BTE calculator in the relaxation time approximation.
    !! It is impure as it mutates the phonon sector of the bte data type and
    !! writes to disk. It should be kept private to this data type unless made safer.

    class(sepe), intent(inout) :: self !Mutation alert!
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(phonon), intent(in) :: ph
    type(electron), intent(in) :: el
    character(*), intent(in) :: Tdir

    !Locals
    real(r64) :: ph_kappa_scalar, ph_kappa_scalar_old, ph_alphabyT_scalar, ph_alphabyT_scalar_old
    type(timer) :: t
    integer(i64) :: iq
    type(transport_coeffs) :: trans

    call trans%initialize_ph(ph%numbands)

    call t%start_timer('RTA ph BTE')

    !Allocate total RTA scattering rates
    allocate(self%ph_rta_rates_ibz(ph%nwv_irred, ph%numbands))

    !Calculate RTA scattering rates

    ! 3-phonon and, optionally, phonon-electron 
    if(num%phe) then
       call calculate_ph_rta_rates(self%ph_rta_rates_3ph_ibz, self%ph_rta_rates_phe_ibz, num, crys, ph, el)
    else
       call calculate_ph_rta_rates(self%ph_rta_rates_3ph_ibz, self%ph_rta_rates_phe_ibz, num, crys, ph)
    end if

    ! 4-ph scattering rates
    call calculate_4ph_rta_rates(self%ph_rta_rates_4ph_ibz, num, crys, ph)

    ! phonon-boundary
    call calculate_bound_scatt_rates(ph%prefix, num%phbound, crys%bound_length, &
         ph%vels, ph%indexlist_irred, self%ph_rta_rates_bound_ibz)

    !Matthiessen's rule without thin-film
    self%ph_rta_rates_ibz = self%ph_rta_rates_3ph_ibz + self%ph_rta_rates_phe_ibz + &
         self%ph_rta_rates_iso_ibz + self%ph_rta_rates_subs_ibz + &
         self%ph_rta_rates_bound_ibz + self%ph_rta_rates_4ph_ibz

    ! phonon-thin-film
    call calculate_thinfilm_scatt_rates(ph%prefix, num%phthinfilm, num%phthinfilm_ballistic, &
         crys%specfac, crys%thinfilm_height, crys%thinfilm_normal, &
         ph%vels, ph%indexlist_irred, self%ph_rta_rates_ibz, self%ph_rta_rates_thinfilm_ibz)

    !Matthiessen's rule with thin-film
    self%ph_rta_rates_ibz = self%ph_rta_rates_ibz + self%ph_rta_rates_thinfilm_ibz

    !Allocate RTA phonon coherence rates
    allocate(self%ph_rta_coherence_rates_ibz(ph%nwv_irred, ph%numbands))

    !Calculate phonon RTA coherence rates
    call calculate_ph_rta_coherence_rates(self%ph_rta_coherence_rates_ibz, num, crys, ph, el)

    !gradT field:
    ! Calculate field term (gradT=>F0)
    call calculate_field_term('ph', 'T', ph%nequiv, ph%ibz2fbz_map, &
         crys%T, 0.0_r64, ph%ens, ph%vels, self%ph_rta_rates_ibz, self%ph_field_term_T)

    ! Symmetrize field term
    do iq = 1, ph%nwv
       self%ph_field_term_T(iq,:,:)=transpose(&
            matmul(ph%symmetrizers(:,:,iq),transpose(self%ph_field_term_T(iq,:,:))))
    end do

    ! RTA solution of BTE
    allocate(self%ph_response_T(ph%nwv, ph%numbands, 3))
    self%ph_response_T = self%ph_field_term_T

    ! Calculate transport coefficient
    call calculate_transport_coeff('ph', 'T', crys%T, 1_i64, 0.0_r64, ph%ens, ph%vels, &
         crys%volume, ph%wvmesh, self%ph_response_T, sym, trans%ph_kappa, trans%dummy)
    !---------------------------------------------------------------------------------!

    !E field:
    ! Calculate field term (E=>G0)
    call calculate_field_term('ph', 'E', ph%nequiv, ph%ibz2fbz_map, &
         crys%T, 0.0_r64, ph%ens, ph%vels, self%ph_rta_rates_ibz, self%ph_field_term_E)

    ! RTA solution of BTE
    allocate(self%ph_response_E(ph%nwv, ph%numbands, 3))
    self%ph_response_E = self%ph_field_term_E

    ! Calculate transport coefficient
    call calculate_transport_coeff('ph', 'E', crys%T, 1_i64, 0.0_r64, ph%ens, ph%vels, &
         crys%volume, ph%wvmesh, self%ph_response_E, sym, trans%ph_alphabyT, trans%dummy)
    trans%ph_alphabyT = trans%ph_alphabyT/crys%T
    !---------------------------------------------------------------------------------!

    !Change to data output directory
    call chdir(trim(adjustl(Tdir)))

    !Write T-dependent RTA scattering rates to file
    call write2file_rank2_real('ph.W_coherence_rta_phe', self%ph_rta_coherence_rates_ibz)
    call write2file_rank2_real('ph.W_rta_3ph', self%ph_rta_rates_3ph_ibz)
    call write2file_rank2_real('ph.W_rta_4ph', self%ph_rta_rates_4ph_ibz)
    call write2file_rank2_real('ph.W_rta_phe', self%ph_rta_rates_phe_ibz)
    call write2file_rank2_real('ph.W_rta', self%ph_rta_rates_ibz)

    !Change back to cwd
    call chdir(trim(adjustl(num%cwd)))

    !Calculate and print transport scalars
    !gradT:
    ph_kappa_scalar = trace(sum(trans%ph_kappa, dim = 1))/crys%dim
    !E:
    ph_alphabyT_scalar = trace(sum(trans%ph_alphabyT, dim = 1))/crys%dim

    !if(.not. num%drag .and. this_image() == 1) then
    if(this_image() == 1) then
       call print_message("RTA solution:")
       call print_message("-------------")
       write(*,*) "iter    k_ph[W/m/K]"
       write(*,"(I3, A, 1E16.8)") 0, "    ", ph_kappa_scalar
    end if

    ph_kappa_scalar_old = ph_kappa_scalar
    ph_alphabyT_scalar_old = ph_alphabyT_scalar

    ! Append RTA coefficients in no-drag files
    ! Change to data output directory
    call chdir(trim(adjustl(Tdir)))
    call append2file_transport_tensor('nodrag_ph_kappa_', 0, trans%ph_kappa)

    ! Print RTA band/branch resolved response functions
    call write2file_response('RTA_F0_', self%ph_response_T) !gradT, ph

    ! Change back to cwd
    call chdir(trim(adjustl(num%cwd)))

    call t%end_timer('RTA ph BTE')

    sync all
  end subroutine dragless_phbte_RTA

  subroutine dragless_ebte_RTA(Tdir, self, num, crys, sym, el, ph)
    !! Dragless electron BTE calculator in the relaxation time approximation.
    !! It is impure as it mutates the electron sector of the bte data type and
    !! writes to disk. It should be kept private to this data type unless made safer.

    class(sepe), intent(inout) :: self !Mutation alert!
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(phonon), intent(in) :: ph
    type(electron), intent(in) :: el
    character(*), intent(in) :: Tdir

    !Locals
    real(r64) :: el_kappa0_scalar, el_kappa0_scalar_old, el_alphabyT_scalar, el_alphabyT_scalar_old, &
         el_sigma_scalar, el_sigma_scalar_old, el_sigmaS_scalar, el_sigmaS_scalar_old
    type(timer) :: t
    integer(i64) :: ik
    type(transport_coeffs) :: trans

    call trans%initialize_el(el%numbands)

    call t%start_timer('RTA e BTE')

    !Calculate RTA scattering rates
    ! e-ph and e-impurity
    call calculate_el_rta_rates(self%el_rta_rates_eph_ibz, self%el_rta_rates_echimp_ibz, &
         self%el_rta_rates_ee_ibz, num, crys, el)

    ! e-boundary
    call calculate_bound_scatt_rates(el%prefix, num%elbound, crys%bound_length, &
         el%vels, el%indexlist_irred, self%el_rta_rates_bound_ibz)

    !Allocate total RTA scattering rates
    allocate(self%el_rta_rates_ibz(el%nwv_irred, el%numbands))

    !Matthiessen's rule
    self%el_rta_rates_ibz = self%el_rta_rates_eph_ibz + self%el_rta_rates_echimp_ibz + &
         self%el_rta_rates_ee_ibz + self%el_rta_rates_bound_ibz

    !gradT field:
    ! Calculate field term (gradT=>I0)
    call calculate_field_term('el', 'T', el%nequiv, el%ibz2fbz_map, &
         crys%T, el%chempot, el%ens, el%vels, self%el_rta_rates_ibz, &
         self%el_field_term_T, el%indexlist)

    ! Symmetrize field term
    do ik = 1, el%nwv
       self%el_field_term_T(ik,:,:) = &
            transpose(matmul(el%symmetrizers(:, :, ik), transpose(self%el_field_term_T(ik, :, :))))
    end do

    ! RTA solution of BTE
    allocate(self%el_response_T(el%nwv, el%numbands, 3))
    self%el_response_T = self%el_field_term_T

    ! Calculate transport coefficient
    call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, el%ens, &
         el%vels, crys%volume, el%wvmesh, self%el_response_T, sym, trans%el_kappa0, trans%el_sigmaS)

    !E field:
    ! Calculate field term (E=>J0)
    call calculate_field_term('el', 'E', el%nequiv, el%ibz2fbz_map, &
         crys%T, el%chempot, el%ens, el%vels, self%el_rta_rates_ibz, &
         self%el_field_term_E, el%indexlist)

    ! Symmetrize field term
    do ik = 1, el%nwv
       self%el_field_term_E(ik, :, :) = &
            transpose(matmul(el%symmetrizers(:, :, ik), transpose(self%el_field_term_E(ik, :, :))))
    end do

    ! RTA solution of BTE
    allocate(self%el_response_E(el%nwv, el%numbands, 3))
    self%el_response_E = self%el_field_term_E

    ! Calculate transport coefficient
    call calculate_transport_coeff('el', 'E', crys%T, el%spindeg, el%chempot, el%ens, el%vels, &
         crys%volume, el%wvmesh, self%el_response_E, sym, trans%el_alphabyT, trans%el_sigma)
    trans%el_alphabyT = trans%el_alphabyT/crys%T
    !--!

    !Change to data output directory
    call chdir(trim(adjustl(Tdir)))

    !Write e-ph RTA scattering rates to file
    call write2file_rank2_real('el.W_rta_eph', self%el_rta_rates_eph_ibz)

    !Write e-e RTA scattering rates to file
    call write2file_rank2_real('el.W_rta_ee', self%el_rta_rates_ee_ibz)

    !Write e-chimp RTA scattering rates to file
    call write2file_rank2_real('el.W_rta_echimp', self%el_rta_rates_echimp_ibz)

    !Change back to cwd
    call chdir(trim(adjustl(num%cwd)))

    !Calculate and print transport scalars
    !gradT:
    el_kappa0_scalar = trace(sum(trans%el_kappa0, dim = 1))/crys%dim
    el_sigmaS_scalar = trace(sum(trans%el_sigmaS, dim = 1))/crys%dim

    !E:
    el_sigma_scalar = trace(sum(trans%el_sigma, dim = 1))/crys%dim
    el_alphabyT_scalar = trace(sum(trans%el_alphabyT, dim = 1))/crys%dim

    !if(.not. num%drag .and. this_image() == 1) then
    if(this_image() == 1) then
       call print_message("RTA solution:")
       call print_message("-------------")
       write(*,*) "iter    k0_el[W/m/K]        sigmaS[A/m/K]", &
            "         sigma[1/Ohm/m]      alpha_el/T[A/m/K]"
       write(*,"(I3, A, 1E16.8, A, 1E16.8, A, 1E16.8, A, 1E16.8)") 0, &
            "    ", el_kappa0_scalar, "     ", el_sigmaS_scalar, &
            "     ", el_sigma_scalar, "     ", el_alphabyT_scalar
    end if

    el_kappa0_scalar_old = el_kappa0_scalar
    el_sigmaS_scalar_old = el_sigmaS_scalar
    el_sigma_scalar_old = el_sigma_scalar 
    el_alphabyT_scalar_old = el_alphabyT_scalar

    ! Append RTA coefficients in no-drag files
    ! Change to data output directory
    call chdir(trim(adjustl(Tdir)))
    call append2file_transport_tensor('nodrag_el_sigmaS_', 0, trans%el_sigmaS, el%bandlist)
    call append2file_transport_tensor('nodrag_el_sigma_', 0, trans%el_sigma, el%bandlist)
    call append2file_transport_tensor('nodrag_el_alphabyT_', 0, trans%el_alphabyT, el%bandlist)
    call append2file_transport_tensor('nodrag_el_kappa0_', 0, trans%el_kappa0, el%bandlist)

    ! Print RTA band/branch resolved response functions
    call write2file_response('RTA_I0_', self%el_response_T, el%bandlist) !gradT, el
    call write2file_response('RTA_J0_', self%el_response_E, el%bandlist) !E, el

    ! Change back to cwd
    call chdir(trim(adjustl(num%cwd)))

    call t%end_timer('RTA e BTE')

    sync all
  end subroutine dragless_ebte_RTA

  subroutine dragfull_ephbtes(Tdir, self, num, crys, sym, ph, el)
    !! Dragful electron-phonon BTEs calculator.
    !! It is impure as it mutates the the bte data type and
    !! writes to disk. It should be kept private to this data type unless made safer.

    class(sepe), intent(inout) :: self !Mutation alert!
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(phonon), intent(in) :: ph
    type(electron), intent(in) :: el
    character(*), intent(in) :: Tdir

    !Locals
    real(r64) :: ph_kappa_scalar, ph_kappa_scalar_old, ph_alphabyT_scalar, ph_alphabyT_scalar_old, &
         el_kappa0_scalar, el_kappa0_scalar_old, el_alphabyT_scalar, el_alphabyT_scalar_old, &
         el_sigma_scalar, el_sigma_scalar_old, el_sigmaS_scalar, el_sigmaS_scalar_old, KO_dev, &
         lambda, lambda_diag(3), tot_alphabyT_scalar
    real(r64), allocatable :: I_diff(:,:,:), I_drag(:,:,:), &
         ph_drag_term_T(:,:,:), ph_drag_term_E(:,:,:), &
         ph_coherence_term_T(:,:,:), ph_coherence_term_E(:,:,:), widc(:,:)
    integer(i64), allocatable :: idc(:,:) , ksint(:,:)
    integer :: it_ph_occ, it_ph_coh, it_el, icart
    integer(i64) :: ik
    character(:), allocatable :: tableheader
    type(timer) :: t
    type(transport_coeffs) :: trans

    call trans%initialize_el(el%numbands)
    call trans%initialize_ph(ph%numbands)

    !Allocate phonon coherences
    allocate(self%ph_coherence_E(ph%nwv, ph%numbands, 3), self%ph_coherence_T(ph%nwv, ph%numbands, 3))
    self%ph_coherence_E = complex_zero
    self%ph_coherence_T = complex_zero

    call t%start_timer('Coupled e-ph BTEs')

    allocate(widc(product(el%wvmesh),6), idc(product(el%wvmesh),9), &
         ksint(product(el%wvmesh),3))
    do ik = 1, size(ksint,1)
       call demux_vector(ik, ksint(ik,:), el%wvmesh, 0_i64)
    end do
    call precompute_interpolation_corners_and_weights(ph%wvmesh,  &
         el%mesh_ref_array, ksint, idc, widc)

    tot_alphabyT_scalar = el_alphabyT_scalar + ph_alphabyT_scalar
    KO_dev = 100.0_r64*abs(&
         (el_sigmaS_scalar - tot_alphabyT_scalar)/tot_alphabyT_scalar)

    ! We need to compute the RTA values again (that is relatively cheap)
    call calculate_transport_coeff('ph', 'T', crys%T, 1_i64, 0.0_r64, ph%ens, ph%vels, &
         crys%volume, ph%wvmesh, self%ph_response_T, sym, trans%ph_kappa, trans%dummy)
    call calculate_transport_coeff('ph', 'E', crys%T, 1_i64, 0.0_r64, ph%ens,  ph%vels, &
         crys%volume, ph%wvmesh, self%ph_response_E, sym, trans%ph_alphabyT, trans%dummy)
    call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, el%ens, &
         el%vels, crys%volume, el%wvmesh, self%el_response_T, sym, trans%el_kappa0, trans%el_sigmaS)
    call calculate_transport_coeff('el', 'E', crys%T, el%spindeg, el%chempot, el%ens, el%vels, &
         crys%volume, el%wvmesh, self%el_response_E, sym, trans%el_alphabyT, trans%el_sigma)
    trans%el_alphabyT = trans%el_alphabyT/crys%T
    trans%ph_alphabyT = trans%ph_alphabyT/crys%T

    ! Change to data output directory
    call chdir(trim(adjustl(Tdir)))
    call append2file_transport_tensor('drag_ph_kappa_', 0, trans%ph_kappa)
    call append2file_transport_tensor('drag_ph_alphabyT_', 0, trans%ph_alphabyT)
    call append2file_transport_tensor('drag_el_sigmaS_', 0, trans%el_sigmaS, el%bandlist)
    call append2file_transport_tensor('drag_el_sigma_', 0, trans%el_sigma, el%bandlist)
    call append2file_transport_tensor('drag_el_alphabyT_', 0, trans%el_alphabyT, el%bandlist)
    call append2file_transport_tensor('drag_el_kappa0_', 0, trans%el_kappa0, el%bandlist)
    ! Change back to cwd
    call chdir(trim(adjustl(num%cwd)))

    call print_message("Coupled electron-phonon transport:")
    call print_message("----------------------------------")

    if(this_image() == 1) then
       tableheader = "iter     k0_el[W/m/K]         sigmaS[A/m/K]         k_ph[W/m/K]"&
            //"         sigma[1/Ohm/m]         alpha_el/T[A/m/K]         alpha_ph/T[A/m/K]"&
            //"         KO dev.[%]"
       write(*,*) trim(tableheader)
    end if

    !These will be needed below
    allocate(I_drag(el%nwv, el%numbands, 3), I_diff(el%nwv, el%numbands, 3), &
         ph_drag_term_T(el%nwv, el%numbands, 3), ph_drag_term_E(el%nwv, el%numbands, 3), &
         ph_coherence_term_T(el%nwv, el%numbands, 3), ph_coherence_term_E(el%nwv, el%numbands, 3))

    ph_coherence_term_E = 0.0
    ph_coherence_term_T = 0.0

    !Start iterator
    do it_ph_occ = 1, num%maxiter       
       !Scheme: for each step of phonon response, fully iterate the electron response.

       !Iterate phonon response once
       call iterate_ph_occupations_eqn(crys%T, num, crys, ph, el, self%ph_rta_rates_ibz, &
            self%ph_field_term_T, self%ph_response_T, self%el_response_T, self%ph_coherence_T)
       call iterate_ph_occupations_eqn(crys%T, num, crys, ph, el, self%ph_rta_rates_ibz, &
            self%ph_field_term_E, self%ph_response_E, self%el_response_E, self%ph_coherence_E)

       !Calculate phonon transport coefficients
       call calculate_transport_coeff('ph', 'T', crys%T, 1_i64, 0.0_r64, ph%ens, ph%vels, &
            crys%volume, ph%wvmesh, self%ph_response_T, sym, trans%ph_kappa, trans%dummy)
       call calculate_transport_coeff('ph', 'E', crys%T, 1_i64, 0.0_r64, ph%ens, ph%vels, &
            crys%volume, ph%wvmesh, self%ph_response_E, sym, trans%ph_alphabyT, trans%dummy)
       trans%ph_alphabyT = trans%ph_alphabyT/crys%T

       !Calculate phonon drag term for the current phBTE iteration.
       call calculate_phonon_drag(num, el, ph, idc, widc, sym, self%el_rta_rates_ibz, &
            self%ph_response_E, ph_drag_term_E)
       call calculate_phonon_drag(num, el, ph, idc, widc, sym, self%el_rta_rates_ibz, &
            self%ph_response_T, ph_drag_term_T)

       !do it_ph_coh = 1, num%maxiter

       !DFG: turn this off for standard dragful BTEs limit
       !Calculate phonon coherence term that will enter the el BTE iteration cycle
       call calculate_ph_coh_term_of_el_BTE(&
            num, el, ph, idc, widc, sym, &
            self%el_rta_rates_ibz, self%ph_coherence_E, ph_coherence_term_E)
       call calculate_ph_coh_term_of_el_BTE(&
            num, el, ph, idc, widc, sym, &
            self%el_rta_rates_ibz, self%ph_coherence_T, ph_coherence_term_T)

       !Iterate phonon coherence equation
       call iterate_ph_coherence_eqn(num, crys, ph, el, &
            self%ph_rta_coherence_rates_ibz, self%ph_response_T, self%el_response_T, self%ph_coherence_T)
       call iterate_ph_coherence_eqn(num, crys, ph, el, &
            self%ph_rta_coherence_rates_ibz, self%ph_response_E, self%el_response_E, self%ph_coherence_E)

       !Iterate electron response all the way
       do it_el = 1, num%maxiter
          !E field:
          call iterate_el_occupations_eqn(num, el, crys, &
               self%el_rta_rates_ibz, self%el_field_term_E, self%el_response_E, &
               ph_drag_term_E, ph_coherence_term_E)

          !delT field:
          call iterate_el_occupations_eqn(num, el, crys, &
               self%el_rta_rates_ibz, self%el_field_term_T, self%el_response_T, &
               ph_drag_term_T, ph_coherence_term_T)
          !Enforce Kelvin-Onsager relation:
          !Fix "diffusion" part
          do icart = 1, 3
             I_diff(:,:,icart) = (el%ens(:,:) - el%chempot)/qe/crys%T*&
                  self%el_response_E(:,:,icart)
          end do

          !TODO: How to generalize the following to include the coherence term?

          !Correct "drag" part
          I_drag = self%el_response_T - I_diff
          call correct_I_drag(I_drag, trace(sum(trans%ph_alphabyT, dim = 1))/crys%dim, lambda)
          self%el_response_T = I_diff + lambda*I_drag

          !Calculate electron transport coefficients
          call calculate_transport_coeff('el', 'E', crys%T, el%spindeg, el%chempot, &
               el%ens, el%vels, crys%volume, el%wvmesh, self%el_response_E, sym, &
               trans%el_alphabyT, trans%el_sigma, Bfield = num%Bfield)
          trans%el_alphabyT = trans%el_alphabyT/crys%T

          !Calculate electron transport coefficients
          call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, &
               el%ens, el%vels, crys%volume, el%wvmesh, self%el_response_T, sym, &
               trans%el_kappa0, trans%el_sigmaS, Bfield = num%Bfield)

          !Calculate electron transport scalars
          el_kappa0_scalar = trace(sum(trans%el_kappa0, dim = 1))/crys%dim
          el_sigmaS_scalar = trace(sum(trans%el_sigmaS, dim = 1))/crys%dim
          el_sigma_scalar = trace(sum(trans%el_sigma, dim = 1))/crys%dim
          el_alphabyT_scalar = trace(sum(trans%el_alphabyT, dim = 1))/crys%dim

          !Check convergence
          if(converged(el_kappa0_scalar_old, el_kappa0_scalar, num%conv_thres) .and. &
               converged(el_sigmaS_scalar_old, el_sigmaS_scalar, num%conv_thres) .and. &
               converged(el_sigma_scalar_old, el_sigma_scalar, num%conv_thres) .and. &
               converged(el_alphabyT_scalar_old, el_alphabyT_scalar, num%conv_thres)) then
             exit
          else
             el_kappa0_scalar_old = el_kappa0_scalar
             el_sigmaS_scalar_old = el_sigmaS_scalar
             el_sigma_scalar_old = el_sigma_scalar
             el_alphabyT_scalar_old = el_alphabyT_scalar
          end if
       end do !el occupation iterator

       !Calculate phonon transport scalar
       ph_kappa_scalar = trace(sum(trans%ph_kappa, dim = 1))/crys%dim
       ph_alphabyT_scalar = trace(sum(trans%ph_alphabyT, dim = 1))/crys%dim

       if(it_ph_occ == 1) then
          !Print RTA band/branch resolved response functions
          ! Change to data output directory
          call chdir(trim(adjustl(Tdir)))
          call write2file_response('partdcpl_I0_', self%el_response_T, el%bandlist) !gradT, el
          call write2file_response('partdcpl_J0_', self%el_response_E, el%bandlist) !E, el
          ! Change back to cwd
          call chdir(trim(adjustl(num%cwd)))
       end if

       tot_alphabyT_scalar = el_alphabyT_scalar + ph_alphabyT_scalar
       KO_dev = 100.0_r64*abs(&
            (el_sigmaS_scalar - tot_alphabyT_scalar)/tot_alphabyT_scalar)

       if(this_image() == 1) then
          write(*,"(I3, A, 1E16.8, A, 1E16.8, A, 1E16.8, A, 1E16.8, &
               A, 1E16.8, A, 1E16.8, A, 1F6.3)") it_ph_occ, "     ", el_kappa0_scalar, &
               "      ", el_sigmaS_scalar, "     ", ph_kappa_scalar, &
               "    ", el_sigma_scalar, "        ", el_alphabyT_scalar, &
               "         ", ph_alphabyT_scalar, "           ", KO_dev
       end if

       !Print out band resolved transport coefficients
       ! Change to data output directory
       call chdir(trim(adjustl(Tdir)))
       call append2file_transport_tensor('drag_ph_kappa_', it_ph_occ, trans%ph_kappa)
       call append2file_transport_tensor('drag_ph_alphabyT_', it_ph_occ, trans%ph_alphabyT)
       call append2file_transport_tensor('drag_el_sigmaS_', it_ph_occ, trans%el_sigmaS, el%bandlist)
       call append2file_transport_tensor('drag_el_sigma_', it_ph_occ, trans%el_sigma, el%bandlist)
       call append2file_transport_tensor('drag_el_alphabyT_', it_ph_occ, trans%el_alphabyT, el%bandlist)
       call append2file_transport_tensor('drag_el_kappa0_', it_ph_occ, trans%el_kappa0, el%bandlist)
       ! Change back to cwd
       call chdir(trim(adjustl(num%cwd)))

!!$       !Iterate phonon coherence equation
!!$       call iterate_ph_coherence_eqn(num, crys, ph, el, &
!!$            self%ph_rta_coherence_rates_ibz, self%ph_response_T, self%el_response_T, self%ph_coherence_T)
!!$       call iterate_ph_coherence_eqn(num, crys, ph, el, &
!!$            self%ph_rta_coherence_rates_ibz, self%ph_response_E, self%el_response_E, self%ph_coherence_E)

       !Check convergence
       if(converged(ph_kappa_scalar_old, ph_kappa_scalar, num%conv_thres) .and. &
            converged(ph_alphabyT_scalar_old, ph_alphabyT_scalar, num%conv_thres)) then

          !Print converged band/branch resolved response functions
          ! Change to data output directory
          call chdir(trim(adjustl(Tdir)))
          call write2file_response('drag_F0_', self%ph_response_T) !gradT, ph
          call write2file_response('drag_I0_', self%el_response_T, el%bandlist) !gradT, el
          call write2file_response('drag_G0_', self%ph_response_E) !E, ph
          call write2file_response('drag_J0_', self%el_response_E, el%bandlist) !E, el
          ! Change back to cwd
          call chdir(trim(adjustl(num%cwd)))
          exit
       else
          ph_kappa_scalar_old = ph_kappa_scalar
          ph_alphabyT_scalar_old = ph_alphabyT_scalar
       end if
    end do

    !Don't need these anymore
    deallocate(I_drag, I_diff, ph_drag_term_T, ph_drag_term_E)
    deallocate(widc, idc, ksint)

    call t%end_timer('Coupled e-ph BTEs')

    sync all

  contains

    subroutine correct_I_drag(I_drag, constraint, lambda)
      !! Subroutine to find scaling correction to I_drag.

      real(r64), intent(in) :: I_drag(:, :, :), constraint
      real(r64), intent(out) :: lambda

      !Internal variables
      integer(i64) :: it, maxiter
      real(r64) :: a, b, sigmaS(size(I_drag(1, :, 1)), 3, 3),&
           thresh, sigmaS_scalar, dummy(size(I_drag(1, :, 1)), 3, 3)

      a = 0.0_r64 !lower bound
      b = 2.0_r64 !upper bound

      maxiter = 100 !maximum number of iterations to try
      thresh = 1.0e-6_r64 !convergence threshold

      do it = 1, maxiter
         lambda = 0.5_r64*(a + b)
         !Calculate electron transport coefficients
         call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, &
              el%ens, el%vels, crys%volume, el%wvmesh, lambda*I_drag, sym, &
              dummy, sigmaS)         
         sigmaS_scalar = trace(sum(sigmaS, dim = 1))/crys%dim

         if(abs(sigmaS_scalar - constraint) < thresh) then
            exit
         else if(abs(sigmaS_scalar) < abs(constraint)) then
            a = lambda
         else
            b = lambda
         end if
      end do
    end subroutine correct_I_drag

    subroutine correct_I_drag_expt(I_drag, constraint, lambda)
      !! Subroutine to find scaling correction to I_drag.
      !
      ! This generalized the previous one by considering each
      ! diagonal element separately. [NOT FULLY TESTED!]

      real(r64), intent(in) :: I_drag(:, :, :), constraint(3, 3)
      real(r64), intent(out) :: lambda(3)

      !Internal variables
      integer(i64) :: it, maxiter, j
      real(r64) :: a(3), b(3), sigmaS(size(I_drag(1, :, 1)), 3, 3), &
           thresh, sigmaS_mat(3, 3), dummy(size(I_drag(1, :, 1)), 3, 3)
      logical :: flag_conv

      a = [0, 0, 0]*0.0_r64 !lower bound
      b = [2, 2, 2]*1.0_r64 !upper bound

      maxiter = 100 !maximum number of iterations to try
      thresh = 1.0e-6_r64 !convergence threshold

      do it = 1, maxiter
         lambda = 0.5_r64*(a + b)

         !Calculate electron transport coefficients
         call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, &
              el%ens, el%vels, crys%volume, el%wvmesh, &
              I_drag*spread(spread(lambda, dim = 1, ncopies = size(I_drag, 1)), dim = 2, ncopies = size(I_drag, 2)), &
              sym, dummy, sigmaS)
         sigmaS_mat = sum(sigmaS, dim = 1)

         flag_conv = .true.
         do j = 1,3 
            if(abs(sigmaS_mat(j, j) - constraint(j, j)) < thresh) then
               cycle  
            else if(abs(sigmaS_mat(j, j)) < abs(constraint(j, j))) then
               a(j) = lambda(j)
            else
               b(j) = lambda(j)
            end if

            flag_conv = .false.
         end do

         if(flag_conv) exit
      end do
    end subroutine correct_I_drag_expt
  end subroutine dragfull_ephbtes

  subroutine iterate_el_occupations_eqn(num, el, crys, rta_rates_ibz, field_term, &
       response_el, ph_drag_term, ph_coherence_term)
    !! Subroutine to iterate the electron BTE one step.
    !! 
    !! T Temperature in K
    !! drag Is drag included?
    !! el Electron object
    !! sym Symmetry
    !! rta_rates_ibz Electron RTA scattering rates
    !! field_term Electron field coupling term
    !! response_el Electron response function
    !! ph_drag_term Phonon drag term
    !! ph_coherence_term Phonon drag term => H[gradT-field] or P[E-field]

    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    real(r64), intent(in) :: rta_rates_ibz(:, :), field_term(:, :, :)
    real(r64), intent(in), optional :: ph_drag_term(:, :, :), ph_coherence_term(:, :, :)
    real(r64), intent(inout) :: response_el(:, :, :)

    !Local variables
    integer(i64) :: nstates_irred, nprocs, chunk, istate, numbands, numbranches, &
         ik_ibz, m, ieq, ik_sym, ik_fbz, iproc, ikp, n, nk, num_active_images, &
         aux, aux2, aux3, aux4, start, end, nprocs_echimp, neg_ik_fbz
    integer(i64), allocatable :: istate_el(:), istate_ph(:), istate_el_echimp(:)

    real(r64) :: tau_ibz
    real(r64), allocatable :: Xphplus(:), Xphminus(:), Xchimp(:), &
         response_el_reduce(:, :, :)
    character(1024) :: filepath_Xphminus, filepath_Xphplus, filepath_Xechimp, tag

    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T

    !Number of electron bands
    numbands = size(rta_rates_ibz, 2)

    !Number of in-window FBZ wave vectors
    nk = size(field_term, 1)

    !Total number of IBZ states
    nstates_irred = size(rta_rates_ibz, 1)*numbands

    if(present(ph_drag_term) .or. present(ph_coherence_term)) then
       !Number of phonon branches
       numbranches = size(ph_drag_term, 2)
    end if

    !Allocate and initialize response reduction array
    allocate(response_el_reduce(nk, numbands, 3))
    response_el_reduce(:, :, :) = 0.0_r64

    !Divide electron states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over electron IBZ states
       do istate = start, end
          !Demux state index into band (m) and wave vector (ik_ibz) indices
          call demux_state(istate, numbands, m, ik_ibz)

          !Apply energy window to initial (IBZ blocks) electron
          if(abs(el%ens_irred(ik_ibz, m) - el%enref) > el%fsthick) cycle

          !RTA lifetime
          tau_ibz = 0.0_r64
          if(rta_rates_ibz(ik_ibz, m) /= 0.0_r64) then
             tau_ibz = 1.0_r64/rta_rates_ibz(ik_ibz, m)
          end if

          !The e-ph (population) bit:

          !Set X+ filename
          write(tag, '(I9)') istate
          filepath_Xphplus = trim(adjustl(num%Xdir))//'/Xplus.istate'//trim(adjustl(tag))

          !Read X+ from file
          call read_transition_probs_e(trim(adjustl(filepath_Xphplus)), nprocs, Xphplus, &
               istate_el, istate_ph)

          !Set X- filename
          write(tag, '(I9)') istate
          filepath_Xphminus = trim(adjustl(num%Xdir))//'/Xminus.istate'//trim(adjustl(tag))

          !Read X- from file
          call read_transition_probs_e(trim(adjustl(filepath_Xphminus)), nprocs, Xphminus)

          !The e-ch. imp. bit:

          !Read Xchimp from file
          if(num%elchimp) then
             !Set Xchimp filename
             write(tag, '(I9)') istate
             filepath_Xechimp = trim(adjustl(num%Xdir))//'/Xchimp.istate'//trim(adjustl(tag))
             call read_transition_probs_e(trim(adjustl(filepath_Xechimp)), nprocs_echimp, Xchimp, &
                  istate_el_echimp)
          end if

          !Sum over the number of equivalent k-points of the IBZ point
          do ieq = 1, el%nequiv(ik_ibz)
             ik_sym = el%ibz2fbz_map(ieq, ik_ibz, 1) !symmetry
             call binsearch(el%indexlist, el%ibz2fbz_map(ieq, ik_ibz, 2), ik_fbz)

             !The e-ph population bit:

             !Sum over scattering processes
             do iproc = 1, nprocs
                !Grab the final electron
                call demux_state(istate_el(iproc), numbands, n, ikp)

                !Self contribution:

                !Find image of final electron wave vector due to the current symmetry
                call binsearch(el%indexlist, el%equiv_map(ik_sym, ikp), aux)

                response_el_reduce(ik_fbz, m, :) = response_el_reduce(ik_fbz, m, :) + &
                     response_el(aux, n, :)*(Xphplus(iproc) + Xphminus(iproc))
             end do

             !Add charged impurity contribution to the self consistent term
             if(num%elchimp) then
                do iproc = 1, nprocs_echimp
                   !Grab the final electron
                   call demux_state(istate_el_echimp(iproc), numbands, n, ikp)

                   !Self contribution:
                   !Find image of final electron wave vector due to the current symmetry
                   call binsearch(el%indexlist, el%equiv_map(ik_sym, ikp), aux)

                   response_el_reduce(ik_fbz, m, :) = response_el_reduce(ik_fbz, m, :) + &
                        response_el(aux, n, :)*Xchimp(iproc)
                end do
             end if

             !Iterate BTE
             response_el_reduce(ik_fbz, m, :) = field_term(ik_fbz, m, :) + &
                  response_el_reduce(ik_fbz, m, :)*tau_ibz
          end do
       end do
    end if

    !Update the response function
    call co_sum(response_el_reduce)
    response_el = response_el_reduce

    if(present(ph_drag_term)) then
       !Drag contribution:
       response_el(:, :, :) = response_el(:, :, :) + ph_drag_term(:, :, :)
    end if

    if(present(ph_coherence_term)) then
       !Coherence contribution:
       response_el(:, :, :) = response_el(:, :, :) + ph_coherence_term(:, :, :)
    end if

    !Symmetrize response function
    do ik_fbz = 1, nk
       response_el(ik_fbz, :, :) = transpose(&
            matmul(el%symmetrizers(:, :, ik_fbz), transpose(response_el(ik_fbz, :, :))))
    end do
  end subroutine iterate_el_occupations_eqn

  subroutine iterate_ph_occupations_eqn(T, num, crys, ph, el, rta_rates_ibz, &
       field_term, response_ph, response_el, coherence_ph)
    !! Subroutine to iterate the phonon occupations equation one step.
    !! 
    !! T Temperature in K
    !! num Numerics object
    !! crys Crystal object
    !! ph Phonon object
    !! el Electron object
    !! rta_rates_ibz Phonon RTA scattering rates
    !! field_term Phonon field coupling term
    !! response_ph Phonon response function
    !! response_el Electron response function
    !! coherence_ph Phonon coherence term

    type(phonon), intent(in) :: ph
    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    real(r64), intent(in) :: T, rta_rates_ibz(:, :), field_term(:, :, :)
    real(r64), intent(in) :: response_el(:, :, :)
    complex(r64), intent(in) :: coherence_ph(:, :, :)
    real(r64), intent(inout) :: response_ph(:, :, :)

    !Local variables
    integer(i64) :: nstates_irred, chunk, istate1, numbranches, s1, &
         iq1_ibz, ieq, iq1_sym, iq1_fbz, iproc, iq2, s2, iq3, s3, nq, &
         num_active_images, numbands, ik, ikp, m, n, nprocs_phe, aux1, aux2, &
         nprocs_3ph_plus, nprocs_3ph_minus, start, end, nprocs_phcoh
    integer(i64), allocatable :: istate2_plus(:), istate3_plus(:), &
         istate2_minus(:), istate3_minus(:), istate_el1(:), istate_el2(:)
    real(r64) :: tau_ibz, Qcoh
    real(r64), allocatable :: Wp(:), Wm(:), Y(:), U(:), response_ph_reduce(:, :, :), &
         coherence_ph_real(:, :, :)
    character(len = 1024) :: filepath_Wm, filepath_Wp, filepath_Y, filepath_U, tag

    !Set output directory of transition probilities
    write(tag, "(E9.3)") T

    !Number of electron bands
    numbands = size(response_el(1,:,1))

    !Number of phonon branches
    numbranches = size(rta_rates_ibz(1,:))

    !Number of FBZ wave vectors
    nq = size(field_term(:,1,1))

    !Total number of IBZ states
    nstates_irred = size(rta_rates_ibz(:,1))*numbranches

    !Allocate and initialize response reduction array
    allocate(response_ph_reduce(nq, numbranches, 3))
    response_ph_reduce(:,:,:) = 0.0_r64

    !Allocate and set the real part of the coherence function
    !allocate(coherence_ph_real(size(coherence_ph, 1), numbranches, 3))
    !coherence_ph_real = real(coherence_ph)

    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    !Only work with the active images
    if(this_image() <= num_active_images) then       
       !Run over first phonon IBZ states
       do istate1 = start, end
          !Demux state index into branch (s) and wave vector (iq1_ibz) indices
          call demux_state(istate1, numbranches, s1, iq1_ibz)

          !Set file tag
          write(tag, '(I9)') istate1

          !RTA lifetime
          tau_ibz = 0.0_r64
          if(rta_rates_ibz(iq1_ibz, s1) /= 0.0_r64) then
             tau_ibz = 1.0_r64/rta_rates_ibz(iq1_ibz, s1)
          end if

          if(num%W_OTF) then
             call calculate_W3ph_OTF(ph, num, istate1, T, &
                  Wm, Wp, istate2_plus, istate3_plus, istate2_minus, istate3_minus)
             nprocs_3ph_plus = size(Wp); nprocs_3ph_minus = size(Wm)
          else
             !Set W+ filename
             filepath_Wp = trim(adjustl(num%Wdir))//'/Wp.istate'//trim(adjustl(tag))

             !Read W+ from file
             if(allocated(Wp)) deallocate(Wp)
             if(allocated(istate2_plus)) deallocate(istate2_plus)
             if(allocated(istate3_plus)) deallocate(istate3_plus)
             call read_transition_probs_e(trim(adjustl(filepath_Wp)), nprocs_3ph_plus, Wp, &
                  istate2_plus, istate3_plus)

             !Set W- filename
             filepath_Wm = trim(adjustl(num%Wdir))//'/Wm.istate'//trim(adjustl(tag))

             !Read W- from file
             if(allocated(Wm)) deallocate(Wm)
             if(allocated(istate2_minus)) deallocate(istate2_minus)
             if(allocated(istate3_minus)) deallocate(istate3_minus)
             call read_transition_probs_e(trim(adjustl(filepath_Wm)), nprocs_3ph_minus, Wm, &
                  istate2_minus, istate3_minus)
          end if

          !if(present(response_el)) then
          if(num%Y_OTF) then
             call calculate_Y_OTF(el, ph, num, crys, istate1, T, Y, istate_el1, istate_el2)
             nprocs_phe = size(Y)
          else
             !Set Y filename
             filepath_Y = trim(adjustl(num%Ydir))//'/Y.istate'//trim(adjustl(tag))

             !Read Y from file
             if(allocated(Y)) deallocate(Y)
             if(allocated(istate_el1)) deallocate(istate_el1)
             if(allocated(istate_el2)) deallocate(istate_el2)
             call read_transition_probs_e(trim(adjustl(filepath_Y)), nprocs_phe, Y, &
                  istate_el1, istate_el2)
          end if
          !end if

          !Set U filename
          filepath_U = trim(adjustl(num%Ydir))//'/U.istate'//trim(adjustl(tag))

          !Read U from file
          if(allocated(U)) deallocate(U)
          call read_transition_probs_e(trim(adjustl(filepath_U)), nprocs_phcoh, U)

          Qcoh = el%spindeg*sum(U)
          
          !Sum over the number of equivalent q-points of the IBZ point
          do ieq = 1, ph%nequiv(iq1_ibz)
             iq1_sym = ph%ibz2fbz_map(ieq, iq1_ibz, 1) !symmetry
             iq1_fbz = ph%ibz2fbz_map(ieq, iq1_ibz, 2) !image due to symmetry

             !Sum over scattering processes
             !Self contribution from plus processes:
             do iproc = 1, nprocs_3ph_plus
                !Grab 2nd and 3rd phonons
                call demux_state(istate2_plus(iproc), numbranches, s2, iq2)
                call demux_state(istate3_plus(iproc), numbranches, s3, iq3)

                response_ph_reduce(iq1_fbz, s1, :) = response_ph_reduce(iq1_fbz, s1, :) + &
                     Wp(iproc)*(response_ph(ph%equiv_map(iq1_sym, iq3), s3, :) - &
                     response_ph(ph%equiv_map(iq1_sym, iq2), s2, :))
             end do

             !Self contribution from minus processes:
             do iproc = 1, nprocs_3ph_minus
                !Grab 2nd and 3rd phonons
                call demux_state(istate2_minus(iproc), numbranches, s2, iq2)
                call demux_state(istate3_minus(iproc), numbranches, s3, iq3)

                response_ph_reduce(iq1_fbz, s1, :) = response_ph_reduce(iq1_fbz, s1, :) + &
                     0.5_r64*Wm(iproc)*(response_ph(ph%equiv_map(iq1_sym, iq3), s3, :) + &
                     response_ph(ph%equiv_map(iq1_sym, iq2), s2, :))
             end do

             !Drag contribution:
             !if(present(response_el)) then
             do iproc = 1, nprocs_phe
                !Grab initial and final electron states
                call demux_state(istate_el1(iproc), numbands, m, ik)
                call demux_state(istate_el2(iproc), numbands, n, ikp)

                !Find image of electron wave vector due to the current symmetry
                call binsearch(el%indexlist, el%equiv_map(iq1_sym, ik), aux1)
                call binsearch(el%indexlist, el%equiv_map(iq1_sym, ikp), aux2)

                response_ph_reduce(iq1_fbz, s1, :) = response_ph_reduce(iq1_fbz, s1, :) + &
                     el%spindeg*Y(iproc)*(response_el(aux2, n, :) - response_el(aux1, m, :))
             end do
             !end if             

             !Coherence contribution:
             response_ph_reduce(iq1_fbz, s1, :) = response_ph_reduce(iq1_fbz, s1, :) + &
                  Qcoh*real(coherence_ph(iq1_fbz, s1, :))

             !Iterate BTE
             response_ph_reduce(iq1_fbz, s1, :) = field_term(iq1_fbz, s1, :) + &
                  response_ph_reduce(iq1_fbz, s1, :)*tau_ibz          
          end do
       end do
    end if

    !Update the response function
    call co_sum(response_ph_reduce)
    response_ph = response_ph_reduce

    !Symmetrize response function
    do iq1_fbz = 1, nq
       response_ph(iq1_fbz,:,:)=transpose(&
            matmul(ph%symmetrizers(:,:,iq1_fbz), transpose(response_ph(iq1_fbz,:,:))))
    end do
  end subroutine iterate_ph_occupations_eqn

  subroutine iterate_ph_coherence_eqn(num, crys, ph, el, &
       rta_coherence_rates_ibz, response_ph, response_el, coherence_ph)
    !! Subroutine to calculate the phonon coherence equation.
    !! 
    !! num Numerics object
    !! crys Crystal object
    !! ph Phonon object
    !! el Electron object
    !! rta_coherence_rates_ibz Phonon RTA coherence rates
    !! response_ph Phonon response function
    !! response_el Electron response function
    !! coherence_ph Phonon coherence term

    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(phonon), intent(in) :: ph
    type(electron), intent(in) :: el
    real(r64), intent(in) :: rta_coherence_rates_ibz(:, :)
    real(r64), intent(in) :: response_el(:, :, :)
    real(r64), intent(in) :: response_ph(:, :, :)
    complex(r64), intent(inout) :: coherence_ph(:, :, :)

    !Local variables
    integer(i64) :: el_nstates_irred, chunk, istate, numbranches, s, &
         ik_ibz, m, ieq, ik_sym, ik_fbz, ikp_fbz_rot, iproc, ikp, n, nk, &
         iq_fbz, iq_ibz, nq, numbands, neg_iq_fbz, q_indvec(3), &
         num_active_images, start, end, nprocs
    integer(i64), allocatable :: istate_el(:), istate_ph(:)
    real(r64) :: ph_en, coherence_rate, prefactor_denom, occ_fac
    complex(r64) :: prefactor
    real(r64), allocatable :: Xphplus(:), Xphminus(:)
    complex(r64), allocatable :: coherence_ph_reduce(:, :, :)
    character(len = 1024) :: filepath_Xphplus, filepath_Xphminus, tag

    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T

    !Number of electron bands
    numbands = el%numbands

    !Number of in-window FBZ wave vectors
    nk = el%nwv

    !Number of phonon FBZ wave vectors
    nq = ph%nwv

    !Number phonon bands
    numbranches = ph%numbands

    !Total number of IBZ electron states
    el_nstates_irred = el%nwv_irred*numbands

    !Allocate and initialize response reduction array
    allocate(coherence_ph_reduce(nq, numbranches, 3))
    coherence_ph_reduce(:,:,:) = 0.0_r64

    !Divide electron states among images
    call distribute_points(el_nstates_irred, chunk, start, end, num_active_images)

    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over electron IBZ states
       do istate = start, end
          !Demux state index into band (m) and wave vector (ik_ibz) indices
          call demux_state(istate, numbands, m, ik_ibz)

          !Apply energy window to initial (IBZ blocks) electron
          if(abs(el%ens_irred(ik_ibz, m) - el%enref) > el%fsthick) cycle

          !Read the state resolved transition rates and the states involved
          !in the e-ph scattering:

          !Set X+ filename
          write(tag, '(I9)') istate
          filepath_Xphplus = trim(adjustl(num%Xdir))//'/Xplus.istate'//trim(adjustl(tag))

          !Read X+ from file
          call read_transition_probs_e(trim(adjustl(filepath_Xphplus)), nprocs, &
               Xphplus, istate_el, istate_ph)

          !Set X- filename
          write(tag, '(I9)') istate
          filepath_Xphminus = trim(adjustl(num%Xdir))//'/Xminus.istate'//trim(adjustl(tag))

          !Read X- from file
          call read_transition_probs_e(trim(adjustl(filepath_Xphminus)), nprocs, &
               Xphminus)

          !Sum over the number of equivalent k-points of the IBZ point
          do ieq = 1, el%nequiv(ik_ibz)
             ik_sym = el%ibz2fbz_map(ieq, ik_ibz, 1) !symmetry
             call binsearch(el%indexlist, el%ibz2fbz_map(ieq, ik_ibz, 2), ik_fbz)

             !Sum over scattering processes
             do iproc = 1, nprocs
                !Grab the final electron
                call demux_state(istate_el(iproc), numbands, n, ikp)

                !Find image of final electron wave vector due to the current symmetry
                call binsearch(el%indexlist, el%equiv_map(ik_sym, ikp), ikp_fbz_rot)

                !Recall that phonons that are not on the coarser q-mesh were tagged with a negative index.
                !Below, I only care about those q-vectors that live on the coarser q-mesh.
                if(istate_ph(iproc) >= 0) then
                   call demux_state(istate_ph(iproc), numbranches, s, iq_fbz)

                   !Compute the prefactor
                   iq_ibz = ph%fbz2ibz_map(iq_fbz)
                   coherence_rate = rta_coherence_rates_ibz(iq_ibz, s)
                   !Note that I don't save the ph%ens_irred to save space. 
                   ph_en = ph%ens(iq_fbz, s)
                   
                   occ_fac = Bose(ph_en, crys%T)
                   occ_fac = occ_fac*(1.0_r64 + occ_fac)

                   prefactor_denom = (coherence_rate)**2 + &
                        4.0_r64*(ph_en/hbar_eVps*occ_fac)**2

                   !Here accumulate contribution to \mathbf{R}_{\lambda}
                   !if(prefactor /= complex_zero) then
                   if(prefactor_denom /= 0.0_r64 .and. ph_en > 0.0_r64) then
                      prefactor = &
                           (coherence_rate - 2.0_r64*ph_en/hbar_eVps*oneI*occ_fac)&
                           /prefactor_denom

                      !from the 1st term of the RHS
                      coherence_ph_reduce(iq_fbz, s, :) = coherence_ph_reduce(iq_fbz, s, :) + &
                           prefactor*(Xphplus(iproc) - Xphminus(iproc))*&
                           (response_el(ik_fbz, m, :) - response_el(ikp_fbz_rot, n, :))

                      !Below I'll need the negative q Umklapped
                      call demux_vector(iq_fbz, q_indvec, ph%wvmesh, 0_i64)

                      !Find index of -q after Umklapping
                      neg_iq_fbz = mux_vector(modulo(-q_indvec, ph%wvmesh), ph%wvmesh, 0_i64)
                      
                      !Now accumulate the 2nd term.
                      !(This does not assume odd parity of F and G.)
                      coherence_ph_reduce(iq_fbz, s, :) = coherence_ph_reduce(iq_fbz, s, :) + &
                           prefactor*(Xphplus(iproc)*response_ph(iq_fbz, s, :) + &
                           Xphminus(iproc)*response_ph(neg_iq_fbz, s, :))
                   end if
                end if

!!$                !Get q = k' - k (represented w.r.t. electron k-grid)
!!$                !q_vec_wrt_elmesh = vec(el%indexlist(ikp_fbz_rot), el%wvmesh, crys%reclattvecs)
!!$
!!$                !Check if the q_vec is in the FBZ q-list.
!!$                !We have to do this since the phonon and electon
!!$                !wave vector meshes are in general different
!!$                compatible_iq_frac = 0
!!$                do iq_dummy = 1, ph%nwv
!!$                   if(all(q_vec_wrt_elmesh%frac == ph%wavevecs(iq_dummy, :))) then
!!$                      compatible_iq_frac = iq_dummy
!!$                   else
!!$                      cycle
!!$                   end if
!!$                end do

!!$                if(compatible_iq_frac > 0) then
!!$                   !Compute the prefactor
!!$                   !TODO
!!$                   
!!$                   !Here accumulate contribution to \mathbf{R}_{\lambda}
!!$                   !from the 1st term of the RHS
!!$                   R_reduce(compatible_iq_frac, s, :) = R_reduce(ik_fbz, m, :) + &
!!$                        prefactor*(Xphplus(iproc) - Xphminus(iproc))*&
!!$                        (response_el(ik_fbz, m, :) - response_el(ikp_fbz_rot, n, :))
!!$
!!$                   !And, similarly, accumulate the 2nd term
!!$                   R_reduce(compatible_iq_frac, s, :) = R_reduce(ik_fbz, m, :) + &
!!$                        prefactor*(Xphplus(iproc) - Xphminus(iproc))*&
!!$                        (response_ph(ik_fbz, m, :) - response_el(ikp_fbz_rot, n, :))
!!$                else
!!$                   cycle
!!$                end if
             end do
          end do
       end do
    end if

    !Update the response function
    call co_sum(coherence_ph_reduce)
    coherence_ph = coherence_ph_reduce

    !Symmetrize response function
    do iq_fbz = 1, nq
       coherence_ph(iq_fbz, :, :)=transpose(&
            matmul(ph%symmetrizers(:, :, iq_fbz), &
            transpose(coherence_ph(iq_fbz, :, :))))
    end do
  end subroutine iterate_ph_coherence_eqn

  subroutine calculate_phonon_drag(num, el, ph, idc, widc, sym, rta_rates_ibz, &
       response_ph, ph_drag_term)
    !! Subroutine to calculate the phonon drag term that enters the electron
    !! occupations equation.
    !! 
    !! num Numerics object
    !! el Electron object
    !! ph Phonon object
    !! idc corners in the coarse mesh for the refined mesh.
    !! widc weights of the corners for interpolation of values in corse mesh to the refined one
    !! sym Symmetry
    !! rta_rates_ibz Electron RTA scattering rates
    !! response_ph Phonon response function
    !! ph_drag_term Phonon drag term

    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    type(symmetry), intent(in) :: sym
    integer(i64), intent(in) :: idc(:, :)
    real(r64), intent(in) :: rta_rates_ibz(:, :), response_ph(:, :, :), widc(:, :)
    real(r64), intent(out) :: ph_drag_term(:, :, :)

    !Local variables
    integer(i64) :: nstates_irred, nprocs, chunk, istate, numbands, numbranches, &
         ik_ibz, m, ieq, ik_sym, ik_fbz, iproc, iq, s, nk, num_active_images, &
         fineq_indvec(3), fineq_indvec_im(3), start, end, iq2inter, neg_iq, q_indvec(3)
    integer(i64), allocatable :: istate_el(:), istate_ph(:)
    real(r64) :: tau_ibz, ForG(3), ForG_minusq(3)
    real(r64), allocatable :: Xplus(:), Xminus(:), ph_drag_term_reduce(:, :, :)
    character(1024) :: filepath_Xminus, filepath_Xplus, tag

    !Number of electron bands
    numbands = el%numbands

    !Number of in-window FBZ wave vectors
    nk = el%nwv

    !Total number of IBZ states
    nstates_irred = el%nwv_irred*numbands

    !Number of phonon branches
    numbranches = ph%numbands

    !Allocate and initialize response reduction array
    allocate(ph_drag_term_reduce(nk, numbands, 3))
    ph_drag_term_reduce(:, :, :) = 0.0_r64

    !Divide electron states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over electron IBZ states
       do istate = start, end
          !Demux state index into band (m) and wave vector (ik_ibz) indices
          call demux_state(istate, numbands, m, ik_ibz)

          !Apply energy window to initial (IBZ blocks) electron
          if(abs(el%ens_irred(ik_ibz, m) - el%enref) > el%fsthick) cycle

          !RTA lifetime
          tau_ibz = 0.0_r64
          if(rta_rates_ibz(ik_ibz, m) /= 0.0_r64) then
             tau_ibz = 1.0_r64/rta_rates_ibz(ik_ibz, m)
          end if

          !Set X+ filename
          write(tag, '(I9)') istate
          filepath_Xplus = trim(adjustl(num%Xdir))//'/Xplus.istate'//trim(adjustl(tag))

          !Read X+ from file
          call read_transition_probs_e(trim(adjustl(filepath_Xplus)), nprocs, Xplus, &
               istate_el, istate_ph)

          !Set X- filename
          write(tag, '(I9)') istate
          filepath_Xminus = trim(adjustl(num%Xdir))//'/Xminus.istate'//trim(adjustl(tag))

          !Read X- from file
          call read_transition_probs_e(trim(adjustl(filepath_Xminus)), nprocs, Xminus)

          !Sum over the number of equivalent k-points of the IBZ point
          do ieq = 1, el%nequiv(ik_ibz)
             ik_sym = el%ibz2fbz_map(ieq, ik_ibz, 1) !symmetry
             call binsearch(el%indexlist, el%ibz2fbz_map(ieq, ik_ibz, 2), ik_fbz)

             !Sum over scattering processes
             do iproc = 1, nprocs
                if(istate_ph(iproc) < 0) then !This phonon is on the (fine) electron mesh
                   call demux_state(-istate_ph(iproc), numbranches, s, iq)
                   iq = -iq !Keep the negative tag
                else !This phonon is on the phonon mesh
                   call demux_state(istate_ph(iproc), numbranches, s, iq)
                end if

                !Drag contribution:             
                if(iq < 0) then !Need to interpolate on this point
                   !Calculate the fine mesh wave vector, 0-based index vector
                   call demux_vector(-iq, fineq_indvec, el%wvmesh, 0_i64)

                   !Find image of phonon wave vector due to the current symmetry
                   fineq_indvec_im = modulo( &
                        nint(matmul(sym%qrotations(:, :, ik_sym), fineq_indvec)), el%wvmesh)

                   !Interpolate response function on this wave vector using precomputed tabulated weights
                   !and points. I note that response_ph(:, s, :) is not contiguous in memory.
                   iq2inter = mux_vector(fineq_indvec_im, el%wvmesh, 0_i64)
                   call interpolate_using_precomputed(idc(iq2inter, :), widc(iq2inter, :),&
                        response_ph(:, s, :), ForG(:))

                   !The -q part:
                   !Find image of phonon wave vector due to the current symmetry
                   fineq_indvec_im = modulo(-fineq_indvec_im, el%wvmesh)
                   
                   !Interpolate response function on this wave vector using precomputed tabulated weights
                   !and points. I note that response_ph(:, s, :) is not contiguous in memory.
                   iq2inter = mux_vector(fineq_indvec_im, el%wvmesh, 0_i64)
                   call interpolate_using_precomputed(idc(iq2inter, :), widc(iq2inter, :),&
                        response_ph(:, s, :), ForG_minusq(:))
                else
                   !F(q) or G(q)
                   ForG(:) = response_ph(ph%equiv_map(ik_sym, iq), s, :)

                   !The -q part:
                   !Find image of phonon wave vector due to the current symmetry
                   call demux_vector(ph%equiv_map(ik_sym, iq), q_indvec, ph%wvmesh, 0_i64)

                   !Find index of -q after Umklapping
                   neg_iq = mux_vector(modulo(-q_indvec, ph%wvmesh), ph%wvmesh, 0_i64)

                   ForG_minusq(:) = response_ph(neg_iq, s, :)
                end if
                
!!$                !Here we use the fact that F(-q) = -F(q) and G(-q) = -G(q)
!!$                ph_drag_term_reduce(ik_fbz, m, :) = ph_drag_term_reduce(ik_fbz, m, :) - &
!!$                     ForG(:)*(Xplus(iproc) + Xminus(iproc))

                !Let's not assume the odd parity since the coherence term is even.       
                ph_drag_term_reduce(ik_fbz, m, :) = ph_drag_term_reduce(ik_fbz, m, :) - &
                     ForG(:)*Xplus(iproc) + ForG_minusq(:)*Xminus(iproc)
             end do

             !Multiply life time factor 
             ph_drag_term_reduce(ik_fbz, m, :) = ph_drag_term_reduce(ik_fbz, m, :)*tau_ibz
          end do
       end do
    end if

    !Reduce from all images
    call co_sum(ph_drag_term_reduce)
    ph_drag_term = ph_drag_term_reduce
  end subroutine calculate_phonon_drag

  subroutine calculate_ph_coh_term_of_el_BTE(&
       num, el, ph, coarse_mesh_corners, weights_coarse_mesh_corners, sym, &
       rta_rates_ibz, coherence_ph, ph_coherence_term)
    !! Computes the phonon coherence term that enters the electron
    !! occupations equation.
    !!
    !! ph_coherence Phonon coherence
    !! num Numerics object
    !! el Electron object
    !! ph Phonon object
    !! coarse_mesh_corners corners in the coarse mesh for the refined mesh.
    !! weights_coarse_mesh_corners weights of the corners for interpolation of values in corse mesh to the refined one
    !! sym Symmetry
    !! rta_rates_ibz Electron RTA scattering rates
    !! coherence_ph Phonon coherence function (H[delT] or P[E])
    !! ph_coherence_term Phonon coherence term

    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    type(symmetry), intent(in) :: sym
    integer(i64), intent(in) :: coarse_mesh_corners(:, :)
    real(r64), intent(in) :: rta_rates_ibz(:, :), weights_coarse_mesh_corners(:, :)
    complex(r64), intent(in) :: coherence_ph(:, :, :)
    real(r64), intent(out) :: ph_coherence_term(:, :, :)

    !Local variables
    integer(i64) :: nstates_irred, nprocs_phcoh, &
         chunk, istate, numbands, numbranches, &
         ik_ibz, m, ieq, ik_sym, ik_fbz, iproc, iq, s, nk, num_active_images, &
         fineq_indvec(3), start, end, iq2inter
    integer(i64), allocatable :: istate_el_phcoh(:), istate_ph_phcoh(:)
    real(r64) :: tau_ibz, HorP(3)
    !real(r64), allocatable :: Omegaplus(:), Omegaminus(:), &
    real(r64), allocatable :: Xplus(:), Xminus(:), &
         ph_coherence_term_reduce(:, :, :), coherence_ph_real(:, :, :)
    !character(1024) :: filepath_Omegaminus, filepath_Omegaplus, tag
    character(1024) :: filepath_Xminus, filepath_Xplus, tag

    !Number of electron bands
    numbands = el%numbands

    !Number of in-window FBZ wave vectors
    nk = el%nwv

    !Total number of IBZ states
    nstates_irred = el%nwv_irred*numbands

    !Number of phonon branches
    numbranches = ph%numbands

    !Allocate and initialize response reduction array
    allocate(ph_coherence_term_reduce(nk, numbands, 3))
    ph_coherence_term_reduce = 0.0_r64

    !Allocate and set the real part of the coherence function
    allocate(coherence_ph_real(size(coherence_ph, 1), numbranches, 3))
    coherence_ph_real = real(coherence_ph)

    !Divide electron states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over electron IBZ states
       do istate = start, end
          !Demux state index into band (m) and wave vector (ik_ibz) indices
          call demux_state(istate, numbands, m, ik_ibz)

          !Apply energy window to initial (IBZ blocks) electron
          if(abs(el%ens_irred(ik_ibz, m) - el%enref) > el%fsthick) cycle

          !RTA lifetime
          tau_ibz = 0.0_r64
          if(rta_rates_ibz(ik_ibz, m) /= 0.0_r64) then
             tau_ibz = 1.0_r64/rta_rates_ibz(ik_ibz, m)
          end if

!!$          !Set Omega+ filename
!!$          write(tag, '(I9)') istate
!!$          filepath_Omegaplus = trim(adjustl(num%Xdir))//'/Omegaplus.istate'//trim(adjustl(tag))
          !Set X+ filename
          write(tag, '(I9)') istate
          filepath_Xplus = trim(adjustl(num%Xdir))//'/Xplus.istate'//trim(adjustl(tag))

!!$          !Read Omega+ from file
!!$          call read_transition_probs_e(trim(adjustl(filepath_Omegaplus)), nprocs_phcoh, Omegaplus, &
!!$               istate_el_phcoh, istate_ph_phcoh)
          !Read X+ from file
          call read_transition_probs_e(trim(adjustl(filepath_Xplus)), nprocs_phcoh, Xplus, &
               istate_el_phcoh, istate_ph_phcoh)

!!$          !Set Omega- filename
!!$          write(tag, '(I9)') istate
!!$          filepath_Omegaminus = trim(adjustl(num%Xdir))//'/Omegaminus.istate'//trim(adjustl(tag))
          !Set X- filename
          write(tag, '(I9)') istate
          filepath_Xminus = trim(adjustl(num%Xdir))//'/Xminus.istate'//trim(adjustl(tag))

!!$          !Read Omega- from file
!!$          call read_transition_probs_e(trim(adjustl(filepath_Omegaminus)), nprocs_phcoh, Omegaminus)
          !Read X- from file
          call read_transition_probs_e(trim(adjustl(filepath_Xminus)), nprocs_phcoh, Xminus)

          !Sum over the number of equivalent k-points of the IBZ point
          do ieq = 1, el%nequiv(ik_ibz)
             ik_sym = el%ibz2fbz_map(ieq, ik_ibz, 1) !symmetry
             call binsearch(el%indexlist, el%ibz2fbz_map(ieq, ik_ibz, 2), ik_fbz)

             !Sum over scattering processes
             do iproc = 1, nprocs_phcoh
                if(istate_ph_phcoh(iproc) < 0) then !This phonon is on the (fine) electron mesh
                   call demux_state(-istate_ph_phcoh(iproc), numbranches, s, iq)
                   iq = -iq !Keep the negative tag
                else !This phonon is on the phonon mesh
                   call demux_state(istate_ph_phcoh(iproc), numbranches, s, iq)
                end if

                !Coherence contribution:             
                if(iq < 0) then !Need to interpolate on this point
                   !Calculate the fine mesh wave vector, 0-based index vector
                   call demux_vector(-iq, fineq_indvec, el%wvmesh, 0_i64)

                   !Find image of phonon wave vector due to the current symmetry
                   fineq_indvec = modulo( &
                        nint(matmul(sym%qrotations(:, :, ik_sym), fineq_indvec)), el%wvmesh)

                   !Interpolate response function on this wave vector using precomputed tabulated weights
                   !and points. I note that coherence_ph_real(:, s, :) is not contiguous in memory.
                   iq2inter = mux_vector(fineq_indvec,el%wvmesh, 0_i64)
                   call interpolate_using_precomputed(&
                        coarse_mesh_corners(iq2inter,:), &
                        weights_coarse_mesh_corners(iq2inter,:),&
                        coherence_ph_real(:, s, :), HorP(:))
                else
                   !H(q) or P(q)
                   HorP(:) = coherence_ph_real(ph%equiv_map(ik_sym, iq), s, :)
                end if

                !Note that below we use the fact that H and P are even in wave vector
                !This follows from the definition of coherence: Eq. 24 of Stefanucci & Perfetto SciPost 2023.
!!$                ph_coherence_term_reduce(ik_fbz, m, :) = &
!!$                     ph_coherence_term_reduce(ik_fbz, m, :) - &
!!$                     HorP(:)*(Omegaplus(iproc) + Omegaminus(iproc))
                ph_coherence_term_reduce(ik_fbz, m, :) = &
                     ph_coherence_term_reduce(ik_fbz, m, :) - &
                     HorP(:)*(Xplus(iproc) + Xminus(iproc))
             end do

             !Multiply life time factor 
             ph_coherence_term_reduce(ik_fbz, m, :) = &
                  ph_coherence_term_reduce(ik_fbz, m, :)*tau_ibz
          end do
       end do
    end if

    !Reduce from all images
    call co_sum(ph_coherence_term_reduce)
    ph_coherence_term = ph_coherence_term_reduce
  end subroutine calculate_ph_coh_term_of_el_BTE

  pure logical function converged(oldval, newval, thres)
    !! Function to check if newval is the same as oldval

    real(r64), intent(in) :: oldval, newval, thres

    converged = .False.

    if(newval == oldval) then
       converged = .True.
    else if(oldval /= 0.0_r64) then
       if(abs(newval - oldval)/abs(oldval) < thres) converged = .True.
    end if
  end function converged

end module SEPE_module
