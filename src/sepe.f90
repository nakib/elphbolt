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
  !! Semiconductor Electron-Phonon Equation (SEPE) a la Stefanucci & Perfetto.

  use precision, only: r64, i64
  use params, only: qe, kB, hbar_eVps
  use misc, only: Bose
!!$  use misc, only: print_message, exit_with_message, write2file_rank2_real, &
!!$       distribute_points, demux_state, binsearch, interpolate, demux_vector, mux_vector, &
!!$       trace, subtitle, append2file_transport_tensor, write2file_response, &
!!$       linspace, readfile_response, write2file_spectral_tensor, subtitle, timer, &
!!$       twonorm, write2file_rank1_real, precompute_interpolation_corners_and_weights, &
!!$       interpolate_using_precomputed, Jacobian, cross_product, qdist
  use numerics_module, only: numerics
  use crystal_module, only: crystal
!!$  use nano_module, only: nanostructure
  use symmetry_module, only: symmetry
  use phonon_module, only: phonon
  use electron_module, only: electron
!!$  use interactions, only: calculate_ph_rta_rates, read_transition_probs_e, &
!!$       calculate_el_rta_rates, calculate_bound_scatt_rates, calculate_thinfilm_scatt_rates, &
!!$       calculate_4ph_rta_rates, calculate_W3ph_OTF, calculate_Y_OTF

  implicit none
  
  private
  public sepe

  type sepe
     !! Data and procedures related to the BTE.

     real(r64), allocatable :: ph_rta_rates_phe_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to ph-e interactions.
     real(r64), allocatable :: ph_rta_rates_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ.
     real(r64), allocatable :: ph_field_term_T(:,:,:)
     !! Phonon field coupling term for gradT field on the FBZ.
     real(r64), allocatable :: ph_response_T(:,:,:)
     !! Phonon response function for gradT field on the FBZ.
     real(r64), allocatable :: ph_field_term_E(:,:,:)
     !! Phonon field coupling term for E field on the FBZ.
     real(r64), allocatable :: ph_response_E(:,:,:)
     !! Phonon response function for E field on the FBZ.
     real(r64), allocatable :: ph_coherence(:, :)
     !! Phonon coherence
     
     real(r64), allocatable :: el_rta_rates_echimp_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to charged impurity scattering.
     real(r64), allocatable :: el_rta_rates_bound_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to boundary scattering.
     real(r64), allocatable :: el_rta_rates_eph_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to e-ph interactions.
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

contains

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

    !Local variables
    character(1024) :: tag, Tdir

    call subtitle("Calculating SEPEs...")

    !call print_message("Only the trace-averaged transport coefficients are printed below:")

    !Create output folder tagged by temperature and create it
    write(tag, "(E9.3)") crys%T
    Tdir = trim(adjustl(num%cwd))//'/T'//trim(adjustl(tag))
    if(this_image() == 1) then
       call system('mkdir -p '//trim(adjustl(Tdir)))
    end if
    sync all
    
    !Electron RTA
    call dragless_el_eqn(Tdir, self, num, crys, sym, el, ph)
  end subroutine sepe_driver

  subroutine dragless_el_eqn(Tdir, self, num, crys, sym, el, ph)
    !! Electron transport equation with phonons at equilibrium
    
    class(sepe), intent(inout) :: self !Mutation alert!
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    character(*), intent(in) :: Tdir

    !Locals
    real(r64) :: el_kappa0_scalar, el_kappa0_scalar_old, el_alphabyT_scalar, el_alphabyT_scalar_old, &
         el_sigma_scalar, el_sigma_scalar_old, el_sigmaS_scalar, el_sigmaS_scalar_old, &
         sepe_term(ph%nwv, ph%numbands)
!!$    type(timer) :: t
    integer :: it_el, icart, s, iq
!!$    type(transport_coeffs) :: trans

!!$    call trans%initialize_el(el%numbands)

    call t%start_timer('Iterative electron sector of SEPE')

    call print_message("Dragless electron transport:")
    call print_message("-----------------------------")

    !Restart with RTA solution
    self%el_response_T = self%el_field_term_T
    self%el_response_E = self%el_field_term_E

    !Calculate the "sepe term": 1 + theta/n0
    do s = 1, ph%numbands
       do iq = 1, ph%nwv
          sepe_term(iq, s) = 1.0_r64 + &
               self%ph_coherence(iq, s)/Bose(ph%ens(iq, s), crys%T)
       end do
    end do

!!$    if(this_image() == 1) then
!!$       write(*,*) "iter    k0_el[W/m/K]        sigmaS[A/m/K]", &
!!$            "         sigma[1/Ohm/m]      alpha_el/T[A/m/K]"
!!$    end if

    do it_el = 1, num%maxiter
       !E field:
       call iterate_el_eqn(num, el, crys, &
            self%el_rta_rates_ibz, self%el_field_term_E, self%el_response_E, sepe_term)

!!$       !Calculate electron transport coefficients
!!$       call calculate_transport_coeff('el', 'E', crys%T, el%spindeg, el%chempot, &
!!$            el%ens, el%vels, crys%volume, el%wvmesh, self%el_response_E, sym, &
!!$            trans%el_alphabyT, trans%el_sigma, Bfield = num%Bfield)
!!$       trans%el_alphabyT = trans%el_alphabyT/crys%T

       !delT field:
       call iterate_el_eqn(num, el, crys, &
            self%el_rta_rates_ibz, self%el_field_term_T, self%el_response_T, sepe_term)

       !Enforce Kelvin-Onsager relation
       do icart = 1, 3
          self%el_response_T(:,:,icart) = (el%ens(:,:) - el%chempot)/qe/crys%T*&
               self%el_response_E(:,:,icart)
       end do

!!$       call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, &
!!$            el%ens, el%vels, crys%volume, el%wvmesh, self%el_response_T, sym, &
!!$            trans%el_kappa0, trans%el_sigmaS, Bfield = num%Bfield)
!!$
!!$       !Calculate and print electron transport scalars
!!$       el_kappa0_scalar = trace(sum(trans%el_kappa0, dim = 1))/crys%dim
!!$       el_sigmaS_scalar = trace(sum(trans%el_sigmaS, dim = 1))/crys%dim
!!$       el_sigma_scalar = trace(sum(trans%el_sigma, dim = 1))/crys%dim
!!$       el_alphabyT_scalar = trace(sum(trans%el_alphabyT, dim = 1))/crys%dim
!!$       if(this_image() == 1) then
!!$          write(*,"(I3, A, 1E16.8, A, 1E16.8, A, 1E16.8, A, 1E16.8)") it_el, &
!!$               "    ", el_kappa0_scalar, "     ", el_sigmaS_scalar, &
!!$               "     ", el_sigma_scalar, "     ", el_alphabyT_scalar
!!$       end if
!!$
!!$       !Print out band resolved transport coefficients
!!$       ! Change to data output directory
!!$       call chdir(trim(adjustl(Tdir)))
!!$       call append2file_transport_tensor('nodrag_el_sigmaS_', it_el, trans%el_sigmaS, el%bandlist)
!!$       call append2file_transport_tensor('nodrag_el_sigma_', it_el, trans%el_sigma, el%bandlist)
!!$       call append2file_transport_tensor('nodrag_el_alphabyT_', it_el, trans%el_alphabyT, el%bandlist)
!!$       call append2file_transport_tensor('nodrag_el_kappa0_', it_el, trans%el_kappa0, el%bandlist)
!!$       ! Change back to cwd
!!$       call chdir(trim(adjustl(num%cwd)))
!!$
!!$       !Check convergence
!!$       if(converged(el_kappa0_scalar_old, el_kappa0_scalar, num%conv_thres) .and. &
!!$            converged(el_sigmaS_scalar_old, el_sigmaS_scalar, num%conv_thres) .and. &
!!$            converged(el_sigma_scalar_old, el_sigma_scalar, num%conv_thres) .and. &
!!$            converged(el_alphabyT_scalar_old, el_alphabyT_scalar, num%conv_thres)) then
!!$
!!$          !Print converged band resolved response functions
!!$          ! Change to data output directory
!!$          call chdir(trim(adjustl(Tdir)))
!!$          call write2file_response('nodrag_I0_', self%el_response_T, el%bandlist) !gradT, el
!!$          call write2file_response('nodrag_J0_', self%el_response_E, el%bandlist) !E, el
!!$          ! Change back to cwd
!!$          call chdir(trim(adjustl(num%cwd)))
!!$
!!$          exit
!!$       else
!!$          el_kappa0_scalar_old = el_kappa0_scalar
!!$          el_sigmaS_scalar_old = el_sigmaS_scalar
!!$          el_sigma_scalar_old = el_sigma_scalar
!!$          el_alphabyT_scalar_old = el_alphabyT_scalar
!!$       end if
    end do
!!$
!!$    call t%end_timer('Iterative dragless electron transport')

    sync all
  end subroutine dragless_el_eqn

  subroutine iterate_el_eqn(num, el, crys, rta_rates_ibz, field_term, &
       response_el, sepe_term)
    !! Subroutine to calculate the Fan-Migdal term
    !! 
    !! T Temperature in K
    !! drag Is drag included?
    !! el Electron object
    !! num Numerics object
    !! crys Crystal object
    !! rta_rates_ibz Electron RTA scattering rates
    !! field_term Electron field coupling term
    !! response_el Electron response function
    !! sepe_term Scaling factor of electron transition rates due to phonon coherence
    
    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    real(r64), intent(in) :: rta_rates_ibz(:,:), field_term(:,:,:)
    real(r64), intent(in) :: sepe_term(:,:,:)
    real(r64), intent(inout) :: response_el(:,:,:)

    !Local variables
    integer(i64) :: nstates_irred, nprocs, chunk, istate, numbands, numbranches, &
         ik_ibz, m, ieq, ik_sym, ik_fbz, iproc, ikp, n, nk, num_active_images, aux, &
         start, end, neg_ik_fbz
    integer :: i, j
    integer(i64), allocatable :: istate_el(:), istate_ph(:), istate_el_echimp(:)
    real(r64) :: tau_ibz
    real(r64), allocatable :: Xplus(:), Xminus(:),  Xchimp(:), response_el_reduce(:,:,:)
    character(1024) :: filepath_Xminus, filepath_Xplus, tag
    
    !Set output directory of transition probilities
    write(tag, "(E9.3)") crys%T
    
    !Number of electron bands
    numbands = size(rta_rates_ibz(1,:))

    !Number of in-window FBZ wave vectors
    nk = size(field_term(:,1,1))

    !Total number of IBZ states
    nstates_irred = size(rta_rates_ibz(:,1))*numbands
    
    !Allocate and initialize response reduction array
    allocate(response_el_reduce(nk, numbands, 3))
    response_el_reduce(:,:,:) = 0.0_r64
    
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

          !X+ -> X+(1 + theta_q/n0_q)
          Xplus = Xplus*sepe_term

          !Set X- filename
          write(tag, '(I9)') istate
          filepath_Xminus = trim(adjustl(num%Xdir))//'/Xminus.istate'//trim(adjustl(tag))

          !Read X- from file
          call read_transition_probs_e(trim(adjustl(filepath_Xminus)), nprocs, Xminus)

          !X- -> X-(1 + theta_q/n0_q)
          Xminus = Xminus*sepe_term

          !Sum over the number of equivalent k-points of the IBZ point
          do ieq = 1, el%nequiv(ik_ibz)
             ik_sym = el%ibz2fbz_map(ieq, ik_ibz, 1) !symmetry
             call binsearch(el%indexlist, el%ibz2fbz_map(ieq, ik_ibz, 2), ik_fbz)

             !Sum over scattering processes
             do iproc = 1, nprocs
                !Grab the final electron and, if needed, the interacting phonon
                call demux_state(istate_el(iproc), numbands, n, ikp)

                !Self contribution:

                !Find image of final electron wave vector due to the current symmetry
                call binsearch(el%indexlist, el%equiv_map(ik_sym, ikp), aux)

                response_el_reduce(ik_fbz, m, :) = response_el_reduce(ik_fbz, m, :) + &
                     response_el(aux, n, :)*(Xplus(iproc) + Xminus(iproc))
             end do
                  
             !Iterate BTE
             response_el_reduce(ik_fbz, m, :) = field_term(ik_fbz, m, :) + &
                  response_el_reduce(ik_fbz, m, :)*tau_ibz
          end do
       end do
    end if
    
    !Update the response function
    call co_sum(response_el_reduce)
    response_el = response_el_reduce
  end subroutine iterate_el_eqn

!!$  !TODO Here set the phonon coherence
!!$  subroutine calculate_ph_coherenece
!!$  end subroutine calculate_ph_coherenece
end module SEPE_module
