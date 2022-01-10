! Copyright (C) 2020- Nakib Haider Protik <nakib.haider.protik@gmail.com>
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

module bte_module
  !! Module containing type and procedures related to the solution of the
  !! Boltzmann transport equation (BTE).

  use params, only: dp, k8, qe, kB
  use misc, only: print_message, exit_with_message, write2file_rank2_real, &
       distribute_points, demux_state, binsearch, interpolate, demux_vector, &
       trace, subtitle, append2file_transport_tensor, write2file_response, &
       linspace, readfile_response, write2file_spectral_tensor, subtitle, timer
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use phonon_module, only: phonon
  use electron_module, only: electron
  use interactions, only: calculate_ph_rta_rates, read_transition_probs_e, &
       calculate_el_rta_rates
  use bz_sums, only: calculate_transport_coeff, calculate_spectral_transport_coeff

  implicit none

  private
  public bte

  type bte
     !! Data and procedures related to the BTE.

     real(dp), allocatable :: ph_rta_rates_iso_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to isotope scattering.
     real(dp), allocatable :: ph_rta_rates_subs_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to substitution scattering.
     real(dp), allocatable :: ph_rta_rates_bound_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ due to boundary scattering.
     real(dp), allocatable :: ph_rta_rates_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ.
     real(dp), allocatable :: ph_field_term_T(:,:,:)
     !! Phonon field coupling term for gradT field on the FBZ.
     real(dp), allocatable :: ph_response_T(:,:,:)
     !! Phonon response function for gradT field on the FBZ.
     real(dp), allocatable :: ph_field_term_E(:,:,:)
     !! Phonon field coupling term for E field on the FBZ.
     real(dp), allocatable :: ph_response_E(:,:,:)
     !! Phonon response function for E field on the FBZ.
     
     real(dp), allocatable :: el_rta_rates_echimp_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to charged impurity scattering.
     real(dp), allocatable :: el_rta_rates_bound_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ due to boundary scattering.
     real(dp), allocatable :: el_rta_rates_ibz(:,:)
     !! Electron RTA scattering rates on the IBZ.
     real(dp), allocatable :: el_field_term_T(:,:,:)
     !! Electron field coupling term for gradT field on the FBZ.
     real(dp), allocatable :: el_response_T(:,:,:)
     !! Electron response function for gradT field on the FBZ.
     real(dp), allocatable :: el_field_term_E(:,:,:)
     !! Electron field coupling term for E field on the FBZ.
     real(dp), allocatable :: el_response_E(:,:,:)
     !! Electron response function for E field on the FBZ.
   contains

     procedure :: solve_bte, post_process
     
  end type bte

contains

  subroutine solve_bte(bt, num, crys, sym, ph, el)
    !! Subroutine to solve the BTE
    !!
    !! bt BTE object
    !! num Numerics object
    !! crys Crystal object
    !! sym Symmertry object
    !! ph Phonon object
    !! el Electron object
    
    class(bte), intent(inout) :: bt
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(phonon), intent(in) :: ph
    type(electron), intent(in), optional :: el

    !Local variables
    character(len = 1024) :: tag, Tdir, tableheader
    integer(k8) :: iq, ik, it_ph, it_el, icart
    real(dp), allocatable :: rates_3ph(:,:), rates_phe(:,:), rates_eph(:,:), &
         I_diff(:,:,:), I_drag(:,:,:), ph_kappa(:,:,:), ph_alphabyT(:,:,:), &
         el_sigma(:,:,:), el_sigmaS(:,:,:), el_alphabyT(:,:,:), el_kappa0(:,:,:), &
         dummy(:,:,:), ph_drag_term_T(:,:,:), ph_drag_term_E(:,:,:)
    real(dp) :: ph_kappa_scalar, ph_kappa_scalar_old, el_sigma_scalar, el_sigma_scalar_old, &
         el_sigmaS_scalar, el_sigmaS_scalar_old, el_kappa0_scalar, el_kappa0_scalar_old, &
         ph_alphabyT_scalar, ph_alphabyT_scalar_old, el_alphabyT_scalar, el_alphabyT_scalar_old, &
         KO_dev, tot_alphabyT_scalar, lambda
    type(timer) :: t

    call subtitle("Calculating transport...")

    call print_message("Only the trace-averaged transport coefficients are printed below:")

    !Create output folder tagged by temperature and create it
    write(tag, "(E9.3)") crys%T
    Tdir = trim(adjustl(num%cwd))//'/T'//trim(adjustl(tag))
    if(this_image() == 1) then
       call system('mkdir -p '//trim(adjustl(Tdir)))
    end if
    sync all

    if(.not. num%onlyebte) then
       call t%start_timer('RTA e BTE')
       
       !Allocate phonon transport coefficients
       allocate(ph_kappa(ph%numbands, 3, 3), ph_alphabyT(ph%numbands, 3, 3), &
            dummy(ph%numbands, 3, 3))

       !Allocate total RTA scattering rates
       allocate(bt%ph_rta_rates_ibz(ph%nwv_irred, ph%numbands))

       !Calculate RTA scattering rates
       if(num%phe) then
          call calculate_ph_rta_rates(rates_3ph, rates_phe, num, crys, ph, el)
       else
          call calculate_ph_rta_rates(rates_3ph, rates_phe, num, crys, ph)
       end if
       
       !Matthiessen's rule
       bt%ph_rta_rates_ibz = rates_3ph + rates_phe + &
            bt%ph_rta_rates_iso_ibz + bt%ph_rta_rates_subs_ibz + &
            bt%ph_rta_rates_bound_ibz

       !gradT field:
       ! Calculate field term (gradT=>F0)
       call calculate_field_term('ph', 'T', ph%nequiv, ph%ibz2fbz_map, &
            crys%T, 0.0_dp, ph%ens, ph%vels, bt%ph_rta_rates_ibz, bt%ph_field_term_T)

       ! Symmetrize field term
       do iq = 1, ph%nwv
          bt%ph_field_term_T(iq,:,:)=transpose(&
               matmul(ph%symmetrizers(:,:,iq),transpose(bt%ph_field_term_T(iq,:,:))))
       end do

       ! RTA solution of BTE
       allocate(bt%ph_response_T(ph%nwv, ph%numbands, 3))
       bt%ph_response_T = bt%ph_field_term_T

       ! Calculate transport coefficient
       call calculate_transport_coeff('ph', 'T', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
            crys%volume, ph%wvmesh, bt%ph_response_T, sym, ph_kappa, dummy)
       !---------------------------------------------------------------------------------!

       !E field:
       ! Calculate field term (E=>G0)
       call calculate_field_term('ph', 'E', ph%nequiv, ph%ibz2fbz_map, &
            crys%T, 0.0_dp, ph%ens, ph%vels, bt%ph_rta_rates_ibz, bt%ph_field_term_E)

       ! RTA solution of BTE
       allocate(bt%ph_response_E(ph%nwv, ph%numbands, 3))
       bt%ph_response_E = bt%ph_field_term_E

       ! Calculate transport coefficient
       call calculate_transport_coeff('ph', 'E', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
            crys%volume, ph%wvmesh, bt%ph_response_E, sym, ph_alphabyT, dummy)
       ph_alphabyT = ph_alphabyT/crys%T
       !---------------------------------------------------------------------------------!

       !Change to data output directory
       call chdir(trim(adjustl(Tdir)))

       !Write T-dependent RTA scattering rates to file
       call write2file_rank2_real('ph.W_rta_3ph', rates_3ph)
       call write2file_rank2_real('ph.W_rta_phe', rates_phe)
       call write2file_rank2_real('ph.W_rta', bt%ph_rta_rates_ibz)

       !Change back to cwd
       call chdir(trim(adjustl(num%cwd)))

       !Calculate and print transport scalars
       !gradT:
       ph_kappa_scalar = trace(sum(ph_kappa, dim = 1))/crys%dim
       !E:
       ph_alphabyT_scalar = trace(sum(ph_alphabyT, dim = 1))/crys%dim

       if(.not. num%drag .and. this_image() == 1) then
          write(*,*) "iter    k_ph[W/m/K]"
          write(*,"(I3, A, 1E16.8)") 0, "    ", ph_kappa_scalar
       end if
       
       ph_kappa_scalar_old = ph_kappa_scalar
       ph_alphabyT_scalar_old = ph_alphabyT_scalar
       
       ! Append RTA coefficients in no-drag files
       ! Change to data output directory
       call chdir(trim(adjustl(Tdir)))
       call append2file_transport_tensor('nodrag_ph_kappa_', 0_k8, ph_kappa)

       ! Print RTA band/branch resolved response functions
       call write2file_response('RTA_F0_', bt%ph_response_T) !gradT, ph

       ! Change back to cwd
       call chdir(trim(adjustl(num%cwd)))

       call t%end_timer('RTA e BTE')
    end if

    if(.not. num%onlyphbte) then
       call t%start_timer('RTA ph BTE')
       
       !Allocate electron transport coefficients
       allocate(el_sigma(el%numbands, 3, 3), el_sigmaS(el%numbands, 3, 3), &
            el_alphabyT(el%numbands, 3, 3), el_kappa0(el%numbands, 3, 3))
       
       !Calculate RTA scattering rates
       call calculate_el_rta_rates(rates_eph, bt%el_rta_rates_echimp_ibz, num, crys, el)

       !Allocate total RTA scattering rates
       allocate(bt%el_rta_rates_ibz(el%nwv_irred, el%numbands))

       !Matthiessen's rule
       bt%el_rta_rates_ibz = rates_eph + bt%el_rta_rates_echimp_ibz + &
            bt%el_rta_rates_bound_ibz

       !gradT field:
       ! Calculate field term (gradT=>I0)
       call calculate_field_term('el', 'T', el%nequiv, el%ibz2fbz_map, &
            crys%T, el%chempot, el%ens, el%vels, bt%el_rta_rates_ibz, &
            bt%el_field_term_T, el%indexlist)

       ! Symmetrize field term
       do ik = 1, el%nwv
          bt%el_field_term_T(ik,:,:)=transpose(&
               matmul(el%symmetrizers(:,:,ik),transpose(bt%el_field_term_T(ik,:,:))))
       end do

       ! RTA solution of BTE
       allocate(bt%el_response_T(el%nwv, el%numbands, 3))
       bt%el_response_T = bt%el_field_term_T

       ! Calculate transport coefficient
       call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, el%wvmesh, bt%el_response_T, sym, el_kappa0, el_sigmaS)

       !E field:
       ! Calculate field term (E=>J0)
       call calculate_field_term('el', 'E', el%nequiv, el%ibz2fbz_map, &
            crys%T, el%chempot, el%ens, el%vels, bt%el_rta_rates_ibz, &
            bt%el_field_term_E, el%indexlist)

       ! Symmetrize field term
       do ik = 1, el%nwv
          bt%el_field_term_E(ik,:,:)=transpose(&
               matmul(el%symmetrizers(:,:,ik),transpose(bt%el_field_term_E(ik,:,:))))
       end do

       ! RTA solution of BTE
       allocate(bt%el_response_E(el%nwv, el%numbands, 3))
       bt%el_response_E = bt%el_field_term_E

       ! Calculate transport coefficient
       call calculate_transport_coeff('el', 'E', crys%T, el%spindeg, el%chempot, el%ens, el%vels, &
            crys%volume, el%wvmesh, bt%el_response_E, sym, el_alphabyT, el_sigma)
       el_alphabyT = el_alphabyT/crys%T
       !--!

       !Change to data output directory
       call chdir(trim(adjustl(Tdir)))

       !Write RTA scattering rates to file
       call write2file_rank2_real('el.W_rta_eph', rates_eph)

       !Write e-chimp RTA scattering rates to file
       call write2file_rank2_real('el.W_rta_echimp', bt%el_rta_rates_echimp_ibz)
       
       !Change back to cwd
       call chdir(trim(adjustl(num%cwd)))

       !Calculate and print transport scalars
       !gradT:
       el_kappa0_scalar = trace(sum(el_kappa0, dim = 1))/crys%dim
       el_sigmaS_scalar = trace(sum(el_sigmaS, dim = 1))/crys%dim
       
       !E:
       el_sigma_scalar = trace(sum(el_sigma, dim = 1))/crys%dim
       el_alphabyT_scalar = trace(sum(el_alphabyT, dim = 1))/crys%dim

       if(.not. num%drag .and. this_image() == 1) then
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
       call append2file_transport_tensor('nodrag_el_sigmaS_', 0_k8, el_sigmaS, el%bandlist)
       call append2file_transport_tensor('nodrag_el_sigma_', 0_k8, el_sigma, el%bandlist)
       call append2file_transport_tensor('nodrag_el_alphabyT_', 0_k8, el_alphabyT, el%bandlist)
       call append2file_transport_tensor('nodrag_el_kappa0_', 0_k8, el_kappa0, el%bandlist)

       ! Print RTA band/branch resolved response functions
       call write2file_response('RTA_I0_', bt%el_response_T, el%bandlist) !gradT, el
       call write2file_response('RTA_J0_', bt%el_response_E, el%bandlist) !E, el

       ! Change back to cwd
       call chdir(trim(adjustl(num%cwd)))

       call t%end_timer('RTA ph BTE')
    end if
    
    if(num%drag) then !Coupled BTEs
       call t%start_timer('Coupled e-ph BTEs')
       
       tot_alphabyT_scalar = el_alphabyT_scalar + ph_alphabyT_scalar
       KO_dev = 100.0_dp*abs(&
            (el_sigmaS_scalar - tot_alphabyT_scalar)/tot_alphabyT_scalar)

       call print_message("RTA solution:")
       call print_message("-------------")
       if(this_image() == 1) then
          tableheader = "iter     k0_el[W/m/K]         sigmaS[A/m/K]         k_ph[W/m/K]"&
               //"         sigma[1/Ohm/m]         alpha_el/T[A/m/K]         alpha_ph/T[A/m/K]"&
               //"         KO dev.[%]"
          write(*,*) trim(tableheader)
       end if
       !RTA
       if(this_image() == 1) then
          write(*,"(I3, A, 1E16.8, A, 1E16.8, A, 1E16.8, A, 1E16.8, &
               A, 1E16.8, A, 1E16.8, A, 1F6.3)") 0, "     ", el_kappa0_scalar, &
               "      ", el_sigmaS_scalar, "     ", ph_kappa_scalar, &
               "    ", el_sigma_scalar, "        ", el_alphabyT_scalar, &
               "         ", ph_alphabyT_scalar, "           ", KO_dev
       end if
       
       !Append RTA coefficients in drag files
       ! Change to data output directory
       call chdir(trim(adjustl(Tdir)))
       call append2file_transport_tensor('drag_ph_kappa_', 0_k8, ph_kappa)
       call append2file_transport_tensor('drag_ph_alphabyT_', 0_k8, ph_alphabyT)
       call append2file_transport_tensor('drag_el_sigmaS_', 0_k8, el_sigmaS, el%bandlist)
       call append2file_transport_tensor('drag_el_sigma_', 0_k8, el_sigma, el%bandlist)
       call append2file_transport_tensor('drag_el_alphabyT_', 0_k8, el_alphabyT, el%bandlist)
       call append2file_transport_tensor('drag_el_kappa0_', 0_k8, el_kappa0, el%bandlist)
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
            ph_drag_term_T(el%nwv, el%numbands, 3), ph_drag_term_E(el%nwv, el%numbands, 3))
       
       !Start iterator
       do it_ph = 1, num%maxiter       
          !Scheme: for each step of phonon response, fully iterate the electron response.

          !Iterate phonon response once
          call iterate_bte_ph(crys%T, .True., num, ph, el, bt%ph_rta_rates_ibz, &
               bt%ph_field_term_T, bt%ph_response_T, bt%el_response_T)
          call iterate_bte_ph(crys%T, .True., num, ph, el, bt%ph_rta_rates_ibz, &
               bt%ph_field_term_E, bt%ph_response_E, bt%el_response_E)

          !Calculate phonon transport coefficients
          call calculate_transport_coeff('ph', 'T', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
               crys%volume, ph%wvmesh, bt%ph_response_T, sym, ph_kappa, dummy)
          call calculate_transport_coeff('ph', 'E', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
               crys%volume, ph%wvmesh, bt%ph_response_E, sym, ph_alphabyT, dummy)
          ph_alphabyT = ph_alphabyT/crys%T

          !Calculate phonon drag term for the current phBTE iteration.
          call calculate_phonon_drag(crys%T, num, el, ph, sym, bt%el_rta_rates_ibz, &
               bt%ph_response_E, ph_drag_term_E)
          call calculate_phonon_drag(crys%T, num, el, ph, sym, bt%el_rta_rates_ibz, &
               bt%ph_response_T, ph_drag_term_T)
          
          !Iterate electron response all the way
          do it_el = 1, num%maxiter
             !E field:
             call iterate_bte_el(crys%T, .True., num, el, ph, sym, &
                  bt%el_rta_rates_ibz, bt%el_field_term_E, bt%el_response_E, ph_drag_term_E)

             !Calculate electron transport coefficients
             call calculate_transport_coeff('el', 'E', crys%T, el%spindeg, el%chempot, &
                  el%ens, el%vels, crys%volume, el%wvmesh, bt%el_response_E, sym, &
                  el_alphabyT, el_sigma)
             el_alphabyT = el_alphabyT/crys%T

             !delT field:
             call iterate_bte_el(crys%T, .True., num, el, ph, sym, &
                  bt%el_rta_rates_ibz, bt%el_field_term_T, bt%el_response_T, ph_drag_term_T)
             !Enforce Kelvin-Onsager relation:
             !Fix "diffusion" part
             do icart = 1, 3
                I_diff(:,:,icart) = (el%ens(:,:) - el%chempot)/qe/crys%T*&
                     bt%el_response_E(:,:,icart)
             end do
             !Correct "drag" part
             I_drag = bt%el_response_T - I_diff
             call correct_I_drag(I_drag, trace(sum(ph_alphabyT, dim = 1))/crys%dim, lambda)
             bt%el_response_T = I_diff + lambda*I_drag

             !Calculate electron transport coefficients
             call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, &
                  el%ens, el%vels, crys%volume, el%wvmesh, bt%el_response_T, sym, &
                  el_kappa0, el_sigmaS)

             !Calculate electron transport scalars
             el_kappa0_scalar = trace(sum(el_kappa0, dim = 1))/crys%dim
             el_sigmaS_scalar = trace(sum(el_sigmaS, dim = 1))/crys%dim
             el_sigma_scalar = trace(sum(el_sigma, dim = 1))/crys%dim
             el_alphabyT_scalar = trace(sum(el_alphabyT, dim = 1))/crys%dim

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
          end do

          !Calculate phonon transport scalar
          ph_kappa_scalar = trace(sum(ph_kappa, dim = 1))/crys%dim
          ph_alphabyT_scalar = trace(sum(ph_alphabyT, dim = 1))/crys%dim

          if(it_ph == 1) then
             !Print RTA band/branch resolved response functions
             ! Change to data output directory
             call chdir(trim(adjustl(Tdir)))
             call write2file_response('partdcpl_I0_', bt%el_response_T, el%bandlist) !gradT, el
             call write2file_response('partdcpl_J0_', bt%el_response_E, el%bandlist) !E, el
             ! Change back to cwd
             call chdir(trim(adjustl(num%cwd)))
          end if

          tot_alphabyT_scalar = el_alphabyT_scalar + ph_alphabyT_scalar
          KO_dev = 100.0_dp*abs(&
               (el_sigmaS_scalar - tot_alphabyT_scalar)/tot_alphabyT_scalar)

          if(this_image() == 1) then
             write(*,"(I3, A, 1E16.8, A, 1E16.8, A, 1E16.8, A, 1E16.8, &
                  A, 1E16.8, A, 1E16.8, A, 1F6.3)") it_ph, "     ", el_kappa0_scalar, &
                  "      ", el_sigmaS_scalar, "     ", ph_kappa_scalar, &
                  "    ", el_sigma_scalar, "        ", el_alphabyT_scalar, &
                  "         ", ph_alphabyT_scalar, "           ", KO_dev
          end if

          !Print out band resolved transport coefficients
          ! Change to data output directory
          call chdir(trim(adjustl(Tdir)))
          call append2file_transport_tensor('drag_ph_kappa_', it_ph, ph_kappa)
          call append2file_transport_tensor('drag_ph_alphabyT_', it_ph, ph_alphabyT)
          call append2file_transport_tensor('drag_el_sigmaS_', it_ph, el_sigmaS, el%bandlist)
          call append2file_transport_tensor('drag_el_sigma_', it_ph, el_sigma, el%bandlist)
          call append2file_transport_tensor('drag_el_alphabyT_', it_ph, el_alphabyT, el%bandlist)
          call append2file_transport_tensor('drag_el_kappa0_', it_ph, el_kappa0, el%bandlist)
          ! Change back to cwd
          call chdir(trim(adjustl(num%cwd)))
          
          !Check convergence
          if(converged(ph_kappa_scalar_old, ph_kappa_scalar, num%conv_thres) .and. &
               converged(ph_alphabyT_scalar_old, ph_alphabyT_scalar, num%conv_thres)) then

             !Print converged band/branch resolved response functions
             ! Change to data output directory
             call chdir(trim(adjustl(Tdir)))
             call write2file_response('drag_F0_', bt%ph_response_T) !gradT, ph
             call write2file_response('drag_I0_', bt%el_response_T, el%bandlist) !gradT, el
             call write2file_response('drag_G0_', bt%ph_response_E) !E, ph
             call write2file_response('drag_J0_', bt%el_response_E, el%bandlist) !E, el
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

       call t%end_timer('Coupled e-ph BTEs')
    end if !drag
    
    if(num%onlyphbte .or. num%drag) then !Phonon BTE
       call t%start_timer('Iterative decoupled ph BTE')
       
       call print_message("Decoupled phonon transport:")
       call print_message("---------------------------")

       !Restart with RTA solution
       bt%ph_response_T = bt%ph_field_term_T

       if(this_image() == 1) then
          write(*,*) "iter    k_ph[W/m/K]"
       end if

       do it_ph = 1, num%maxiter
          call iterate_bte_ph(crys%T, .False., num, ph, el, bt%ph_rta_rates_ibz, &
               bt%ph_field_term_T, bt%ph_response_T)

          !Calculate phonon transport coefficients
          call calculate_transport_coeff('ph', 'T', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
               crys%volume, ph%wvmesh, bt%ph_response_T, sym, ph_kappa, dummy)

          !Calculate and print phonon transport scalar
          ph_kappa_scalar = trace(sum(ph_kappa, dim = 1))/crys%dim
          if(this_image() == 1) then
             write(*,"(I3, A, 1E16.8)") it_ph, "    ", ph_kappa_scalar
          end if

          !Print out branch resolved transport coefficients
          ! Change to data output directory
          call chdir(trim(adjustl(Tdir)))
          call append2file_transport_tensor('nodrag_ph_kappa_', it_ph, ph_kappa)
          ! Change back to cwd
          call chdir(trim(adjustl(num%cwd)))
          
          if(converged(ph_kappa_scalar_old, ph_kappa_scalar, num%conv_thres)) then
             !Print converged branch resolved response functions
             ! Change to data output directory
             call chdir(trim(adjustl(Tdir)))
             call write2file_response('nodrag_F0_', bt%ph_response_T) !gradT, ph
             ! Change back to cwd
             call chdir(trim(adjustl(num%cwd)))
             
             exit
          else
             ph_kappa_scalar_old = ph_kappa_scalar
          end if
       end do

       call t%end_timer('Iterative decoupled ph BTE')
    end if

    if(num%onlyebte .or. num%drag) then !Electron BTE
       call t%start_timer('Iterative decoupled e BTE')
       
       call print_message("Decoupled electron transport:")
       call print_message("-----------------------------")

       !Restart with RTA solution
       bt%el_response_T = bt%el_field_term_T
       bt%el_response_E = bt%el_field_term_E
       
       if(this_image() == 1) then
          write(*,*) "iter    k0_el[W/m/K]        sigmaS[A/m/K]", &
               "         sigma[1/Ohm/m]      alpha_el/T[A/m/K]"
       end if

       do it_el = 1, num%maxiter
          !E field:
          call iterate_bte_el(crys%T, .False., num, el, ph, sym,&
               bt%el_rta_rates_ibz, bt%el_field_term_E, bt%el_response_E)

          !Calculate electron transport coefficients
          call calculate_transport_coeff('el', 'E', crys%T, el%spindeg, el%chempot, &
               el%ens, el%vels, crys%volume, el%wvmesh, bt%el_response_E, sym, &
               el_alphabyT, el_sigma)
          el_alphabyT = el_alphabyT/crys%T

          !delT field:
          call iterate_bte_el(crys%T, .False., num, el, ph, sym,&
               bt%el_rta_rates_ibz, bt%el_field_term_T, bt%el_response_T, bt%ph_response_T)
          !Enforce Kelvin-Onsager relation
          do icart = 1, 3
             bt%el_response_T(:,:,icart) = (el%ens(:,:) - el%chempot)/qe/crys%T*&
                  bt%el_response_E(:,:,icart)
          end do

          call calculate_transport_coeff('el', 'T', crys%T, el%spindeg, el%chempot, &
               el%ens, el%vels, crys%volume, el%wvmesh, bt%el_response_T, sym, &
               el_kappa0, el_sigmaS)

          !Calculate and print electron transport scalars
          el_kappa0_scalar = trace(sum(el_kappa0, dim = 1))/crys%dim
          el_sigmaS_scalar = trace(sum(el_sigmaS, dim = 1))/crys%dim
          el_sigma_scalar = trace(sum(el_sigma, dim = 1))/crys%dim
          el_alphabyT_scalar = trace(sum(el_alphabyT, dim = 1))/crys%dim
          if(this_image() == 1) then
             write(*,"(I3, A, 1E16.8, A, 1E16.8, A, 1E16.8, A, 1E16.8)") it_el, &
                  "    ", el_kappa0_scalar, "     ", el_sigmaS_scalar, &
                  "     ", el_sigma_scalar, "     ", el_alphabyT_scalar
          end if

          !Print out band resolved transport coefficients
          ! Change to data output directory
          call chdir(trim(adjustl(Tdir)))
          call append2file_transport_tensor('nodrag_el_sigmaS_', it_el, el_sigmaS, el%bandlist)
          call append2file_transport_tensor('nodrag_el_sigma_', it_el, el_sigma, el%bandlist)
          call append2file_transport_tensor('nodrag_el_alphabyT_', it_el, el_alphabyT, el%bandlist)
          call append2file_transport_tensor('nodrag_el_kappa0_', it_el, el_kappa0, el%bandlist)
          ! Change back to cwd
          call chdir(trim(adjustl(num%cwd)))
          
          !Check convergence
          if(converged(el_kappa0_scalar_old, el_kappa0_scalar, num%conv_thres) .and. &
               converged(el_sigmaS_scalar_old, el_sigmaS_scalar, num%conv_thres) .and. &
               converged(el_sigma_scalar_old, el_sigma_scalar, num%conv_thres) .and. &
               converged(el_alphabyT_scalar_old, el_alphabyT_scalar, num%conv_thres)) then

             !Print converged band resolved response functions
             ! Change to data output directory
             call chdir(trim(adjustl(Tdir)))
             call write2file_response('nodrag_I0_', bt%el_response_T, el%bandlist) !gradT, el
             call write2file_response('nodrag_J0_', bt%el_response_E, el%bandlist) !E, el
             ! Change back to cwd
             call chdir(trim(adjustl(num%cwd)))
             
             exit
          else
             el_kappa0_scalar_old = el_kappa0_scalar
             el_sigmaS_scalar_old = el_sigmaS_scalar
             el_sigma_scalar_old = el_sigma_scalar
             el_alphabyT_scalar_old = el_alphabyT_scalar
          end if
       end do

       call t%end_timer('Iterative decoupled e BTE')
    end if

  contains

    subroutine correct_I_drag(I_drag, constraint, lambda)
      !! Subroutine to find scaling correction to I_drag.

      real(dp), intent(in) :: I_drag(:,:,:), constraint
      real(dp), intent(out) :: lambda

      !Internal variables
      integer(k8) :: it, maxiter
      real(dp) :: a, b, sigmaS(size(I_drag(1,:,1)), 3, 3),&
           thresh, sigmaS_scalar, dummy(size(I_drag(1,:,1)), 3, 3)
      
      a = 0.0_dp !lower bound
      b = 2.0_dp !upper bound
      maxiter = 100
      thresh = 1.0e-6_dp
      do it = 1, maxiter
         lambda = 0.5_dp*(a + b)
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
  end subroutine solve_bte

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
    integer(k8), intent(in) :: nequiv(:), ibz2fbz_map(:,:,:)
    real(dp), intent(in) :: T, chempot, ens(:,:), vels(:,:,:), rta_rates_ibz(:,:)
    real(dp), allocatable, intent(out) :: field_term(:,:,:)
    integer(k8), intent(in), optional :: el_indexlist(:)

    !Local variables
    integer(k8) :: ik_ibz, ik_fbz, ieq, ib, nk_ibz, nk, nbands, pow, &
         chunk, num_active_images
    integer(k8), allocatable :: start[:], end[:]
    real(dp), allocatable :: field_term_reduce(:,:,:)[:]
    real(dp) :: A
    logical :: trivial_case

    !Set constant and power of energy depending on species and field type
    if(species == 'ph') then
       A = 1.0_dp/T
       pow = 1
       if(chempot /= 0.0_dp) then
          call exit_with_message("Phonon chemical potential non-zero in calculate_field_term. Exiting.")
       end if
    else if(species == 'el') then
       if(field == 'T') then
          A = 1.0_dp/T
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
    field_term(:,:,:) = 0.0_dp

    !No field-coupling case
    trivial_case = species == 'ph' .and. field == 'E'

    if(.not. trivial_case) then
       !Allocate start and end coarrays
       allocate(start[*], end[*])

       !Divide IBZ states among images
       call distribute_points(nk_ibz, chunk, start, end, num_active_images)

       !Allocate and initialize field term coarrays
       allocate(field_term_reduce(nk, nbands, 3)[*])
       field_term_reduce(:,:,:) = 0.0_dp

       !Work the active images only:
       do ik_ibz = start, end
          do ieq = 1, nequiv(ik_ibz)
             if(species == 'ph') then
                ik_fbz = ibz2fbz_map(ieq, ik_ibz, 2)
             else
                !Find index of electron in indexlist
                call binsearch(el_indexlist, ibz2fbz_map(ieq, ik_ibz, 2), ik_fbz)
             end if
             do ib = 1, nbands
                if(rta_rates_ibz(ik_ibz, ib) /= 0.0_dp) then
                   field_term_reduce(ik_fbz, ib, :) = A*vels(ik_fbz, ib, :)*&
                        (ens(ik_fbz, ib) - chempot)**pow/rta_rates_ibz(ik_ibz, ib)
                end if
             end do
          end do
       end do
       
       !Reduce field term coarrays
       !Units:
       ! nm.eV/K for phonons, gradT-field
       ! nm.eV/K for electrons, gradT-field
       ! nm.C for electrons, E-field
       call co_sum(field_term_reduce)
       field_term = field_term_reduce
    end if
  end subroutine calculate_field_term

  subroutine iterate_bte_ph(T, drag, num, ph, el, rta_rates_ibz, &
       field_term, response_ph, response_el)
    !! Subroutine to iterate the phonon BTE one step.
    !! 
    !! T Temperature in K
    !! drag Is drag included?
    !! ph Phonon object
    !! rta_rates_ibz Phonon RTA scattering rates
    !! field_term Phonon field coupling term
    !! response_ph Phonon response function
    !! response_el Electron response function

    type(phonon), intent(in) :: ph
    type(electron), intent(in) :: el
    type(numerics), intent(in) :: num
    logical, intent(in) :: drag
    real(dp), intent(in) :: T, rta_rates_ibz(:,:), field_term(:,:,:)
    real(dp), intent(in), optional :: response_el(:,:,:)
    real(dp), intent(inout) :: response_ph(:,:,:)

    !Local variables
    integer(k8) :: nstates_irred, chunk, istate1, numbranches, s1, &
         iq1_ibz, ieq, iq1_sym, iq1_fbz, iproc, iq2, s2, iq3, s3, nq, &
         num_active_images, numbands, ik, ikp, m, n, nprocs_phe, aux1, aux2, &
         nprocs_3ph_plus, nprocs_3ph_minus
    integer(k8), allocatable :: istate2_plus(:), istate3_plus(:), &
         istate2_minus(:), istate3_minus(:), istate_el1(:), istate_el2(:), &
         start[:], end[:]
    real(dp) :: tau_ibz
    real(dp), allocatable :: Wp(:), Wm(:), Y(:), response_ph_reduce(:,:,:)[:]
    character(len = 1024) :: filepath_Wm, filepath_Wp, filepath_Y, tag

    !Set output directory of transition probilities
    write(tag, "(E9.3)") T

    if(drag .and. .not. present(response_el)) then
       call exit_with_message("For drag in phonon BTE, must provide electron response. Exiting.")
    end if
    
    if(drag) then
       !Number of electron bands
       numbands = size(response_el(1,:,1))
    end if
    
    !Number of phonon branches
    numbranches = size(rta_rates_ibz(1,:))

    !Number of FBZ wave vectors
    nq = size(field_term(:,1,1))
    
    !Total number of IBZ states
    nstates_irred = size(rta_rates_ibz(:,1))*numbranches
    
    !Allocate coarrays
    allocate(start[*], end[*])
    allocate(response_ph_reduce(nq, numbranches, 3)[*])

    !Initialize coarray
    response_ph_reduce(:,:,:) = 0.0_dp
    
    !Divide phonon states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
    
    !Run over first phonon IBZ states
    do istate1 = start, end
       !Demux state index into branch (s) and wave vector (iq1_ibz) indices
       call demux_state(istate1, numbranches, s1, iq1_ibz)

       !RTA lifetime
       tau_ibz = 0.0_dp
       if(rta_rates_ibz(iq1_ibz, s1) /= 0.0_dp) then
          tau_ibz = 1.0_dp/rta_rates_ibz(iq1_ibz, s1)
       end if

       !Set W+ filename
       write(tag, '(I9)') istate1
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

       if(drag) then
          !Set Y filename
          filepath_Y = trim(adjustl(num%Ydir))//'/Y.istate'//trim(adjustl(tag))

          !Read Y from file
          if(allocated(Y)) deallocate(Y)
          if(allocated(istate_el1)) deallocate(istate_el1)
          if(allocated(istate_el2)) deallocate(istate_el2)
          call read_transition_probs_e(trim(adjustl(filepath_Y)), nprocs_phe, Y, &
               istate_el1, istate_el2)
       end if

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
                  0.5_dp*Wm(iproc)*(response_ph(ph%equiv_map(iq1_sym, iq3), s3, :) + &
                  response_ph(ph%equiv_map(iq1_sym, iq2), s2, :))
          end do

          !Drag contribution:
          
          if(drag) then
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
          end if

          !Iterate BTE
          response_ph_reduce(iq1_fbz, s1, :) = field_term(iq1_fbz, s1, :) + &
               response_ph_reduce(iq1_fbz, s1, :)*tau_ibz          
       end do
    end do

    !Update the response function
    call co_sum(response_ph_reduce)
    response_ph = response_ph_reduce

    !Symmetrize response function
    do iq1_fbz = 1, nq
       response_ph(iq1_fbz,:,:)=transpose(&
            matmul(ph%symmetrizers(:,:,iq1_fbz),transpose(response_ph(iq1_fbz,:,:))))
    end do
  end subroutine iterate_bte_ph

  subroutine iterate_bte_el(T, drag, num, el, ph, sym, rta_rates_ibz, field_term, &
       response_el, ph_drag_term)
    !! Subroutine to iterate the electron BTE one step.
    !! 
    !! T Temperature in K
    !! drag Is drag included?
    !! el Electron object
    !! ph Phonon object
    !! sym Symmetry
    !! rta_rates_ibz Electron RTA scattering rates
    !! field_term Electron field coupling term
    !! response_el Electron response function
    !! ph_drag_term Phonon drag term
    
    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    type(symmetry), intent(in) :: sym
    logical, intent(in) :: drag
    real(dp), intent(in) :: T, rta_rates_ibz(:,:), field_term(:,:,:)
    real(dp), intent(in), optional :: ph_drag_term(:,:,:)
    real(dp), intent(inout) :: response_el(:,:,:)

    !Local variables
    integer(k8) :: nstates_irred, nprocs, chunk, istate, numbands, numbranches, &
         ik_ibz, m, ieq, ik_sym, ik_fbz, iproc, ikp, n, nk, num_active_images, aux
    integer(k8), allocatable :: istate_el(:), istate_ph(:), start[:], end[:]
    real(dp) :: tau_ibz
    real(dp), allocatable :: Xplus(:), Xminus(:), response_el_reduce(:,:,:)[:]
    character(1024) :: filepath_Xminus, filepath_Xplus, tag

    !Set output directory of transition probilities
    write(tag, "(E9.3)") T

    if(drag .and. .not. present(ph_drag_term)) then
       call exit_with_message("For drag in electron BTE, must provide phonon drag term. Exiting.")
    end if
    
    !Number of electron bands
    numbands = size(rta_rates_ibz(1,:))

    !Number of in-window FBZ wave vectors
    nk = size(field_term(:,1,1))

    !Total number of IBZ states
    nstates_irred = size(rta_rates_ibz(:,1))*numbands

    if(drag) then
       !Number of phonon branches
       !numbranches = size(response_ph(1,:,1))
       numbranches = size(ph_drag_term(1,:,1))
    end if

    !Allocate coarrays
    allocate(start[*], end[*])
    allocate(response_el_reduce(nk, numbands, 3)[*])

    !Initialize coarray
    response_el_reduce(:,:,:) = 0.0_dp
    
    !Divide electron states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    !Run over electron IBZ states
    do istate = start, end
       !Demux state index into band (m) and wave vector (ik_ibz) indices
       call demux_state(istate, numbands, m, ik_ibz)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(el%ens_irred(ik_ibz, m) - el%enref) > el%fsthick) cycle

       !RTA lifetime
       tau_ibz = 0.0_dp
       if(rta_rates_ibz(ik_ibz, m) /= 0.0_dp) then
          tau_ibz = 1.0_dp/rta_rates_ibz(ik_ibz, m)
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

    !Update the response function
    call co_sum(response_el_reduce)
    response_el = response_el_reduce

    if(drag) then
       !Drag contribution:
       response_el(:,:,:) = response_el(:,:,:) + ph_drag_term(:,:,:)
    end if
       
    !Symmetrize response function
    do ik_fbz = 1, nk
       response_el(ik_fbz,:,:)=transpose(&
            matmul(el%symmetrizers(:,:,ik_fbz),transpose(response_el(ik_fbz,:,:))))
    end do
  end subroutine iterate_bte_el

  subroutine calculate_phonon_drag(T, num, el, ph, sym, rta_rates_ibz, response_ph, ph_drag_term)
    !! Subroutine to calculate the phonon drag term.
    !! 
    !! T Temperature in K
    !! num Numerics object
    !! el Electron object
    !! ph Phonon object
    !! sym Symmetry
    !! rta_rates_ibz Electron RTA scattering rates
    !! response_ph Phonon response function
    !! ph_drag_term Phonon drag term

    type(electron), intent(in) :: el
    type(phonon), intent(in) :: ph
    type(numerics), intent(in) :: num
    type(symmetry), intent(in) :: sym
    real(dp), intent(in) :: T, rta_rates_ibz(:,:), response_ph(:,:,:)
    real(dp), intent(out) :: ph_drag_term(:,:,:)

    !Local variables
    integer(k8) :: nstates_irred, nprocs, chunk, istate, numbands, numbranches, &
         ik_ibz, m, ieq, ik_sym, ik_fbz, iproc, iq, s, nk, num_active_images, &
         ipol, fineq_indvec(3)
    integer(k8), allocatable :: istate_el(:), istate_ph(:), start[:], end[:]
    real(dp) :: tau_ibz, ForG(3)
    real(dp), allocatable :: Xplus(:), Xminus(:), ph_drag_term_reduce(:,:,:)[:]
    character(1024) :: filepath_Xminus, filepath_Xplus, tag
    
    !Number of electron bands
    numbands = el%numbands

    !Number of in-window FBZ wave vectors
    nk = el%nwv

    !Total number of IBZ states
    nstates_irred = el%nwv_irred*numbands
    
    !Number of phonon branches
    numbranches = ph%numbands

    !Allocate coarrays
    allocate(start[*], end[*])
    allocate(ph_drag_term_reduce(nk, numbands, 3)[*])
    
    !Initialize coarray
    ph_drag_term_reduce(:,:,:) = 0.0_dp

    !Divide electron states among images
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    !Run over electron IBZ states
    do istate = start, end
       !Demux state index into band (m) and wave vector (ik_ibz) indices
       call demux_state(istate, numbands, m, ik_ibz)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(el%ens_irred(ik_ibz, m) - el%enref) > el%fsthick) cycle

       !RTA lifetime
       tau_ibz = 0.0_dp
       if(rta_rates_ibz(ik_ibz, m) /= 0.0_dp) then
          tau_ibz = 1.0_dp/rta_rates_ibz(ik_ibz, m)
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
                call demux_vector(-iq, fineq_indvec, el%wvmesh, 0_k8)

                !Find image of phonon wave vector due to the current symmetry
                fineq_indvec = modulo( &
                     nint(matmul(sym%qrotations(:, :, ik_sym), fineq_indvec)), el%wvmesh)

                !Interpolate response function on this wave vector
                do ipol = 1, 3
                   call interpolate(ph%wvmesh, el%mesh_ref_array, response_ph(:, s, ipol), &
                        fineq_indvec, ForG(ipol))
                end do
             else
                !F(q) or G(q)
                ForG(:) = response_ph(ph%equiv_map(ik_sym, iq), s, :)
             end if
             !Here we use the fact that F(-q) = -F(q) and G(-q) = -G(q)
             ph_drag_term_reduce(ik_fbz, m, :) = ph_drag_term_reduce(ik_fbz, m, :) - &
                  ForG(:)*(Xplus(iproc) + Xminus(iproc))
          end do
          
          !Multiply life time factor 
          ph_drag_term_reduce(ik_fbz, m, :) = ph_drag_term_reduce(ik_fbz, m, :)*tau_ibz
       end do
    end do

    !Reduce from all images
    call co_sum(ph_drag_term_reduce)
    ph_drag_term = ph_drag_term_reduce
  end subroutine calculate_phonon_drag

  pure logical function converged(oldval, newval, thres)
    !! Function to check if newval is the same as oldval

    real(dp), intent(in) :: oldval, newval, thres

    converged = .False.

    if(newval == oldval) then
       converged = .True.
    else if(oldval /= 0.0_dp) then
       if(abs(newval - oldval)/abs(oldval) < thres) converged = .True.
    end if
  end function converged

  subroutine post_process(bt, num, crys, sym, ph, el)
    !! Subroutine to post-process results of the BTEs.

    class(bte), intent(inout) :: bt
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(symmetry), intent(in) :: sym
    type(phonon), intent(in) :: ph
    type(electron), intent(in), optional :: el

    !Local variables
    real(dp), allocatable :: ph_en_grid(:), el_en_grid(:), ph_kappa(:,:,:,:), dummy(:,:,:,:), &
         el_kappa0(:,:,:,:), el_sigmaS(:,:,:,:), el_sigma(:,:,:,:), el_alphabyT(:,:,:,:), &
         ph_alphabyT(:,:,:,:)
    character(len = 1024) :: tag, Tdir
    integer(k8) :: ie

    !Calculate electron and/or phonon sampling energy grid
    call linspace(ph_en_grid, num%ph_en_min, num%ph_en_max, num%ph_en_num)
    call linspace(el_en_grid, num%el_en_min, num%el_en_max, num%el_en_num)

    !Write energy grids to file
    if(this_image() == 1) then
       open(1, file = "ph.en_grid", status = "replace")
       do ie = 1, num%ph_en_num
          write(1, "(E20.10)") ph_en_grid(ie)
       end do
       close(1)
       open(1, file = "el.en_grid", status = "replace")
       do ie = 1, num%el_en_num
          write(1, "(E20.10)") el_en_grid(ie)
       end do
       close(1)
    end if

    !Change to T-dependent directory
    write(tag, "(E9.3)") crys%T
    Tdir = trim(adjustl(num%cwd))//'/T'//trim(adjustl(tag))
    call chdir(trim(adjustl(Tdir)))

    !Decoupled electron BTE
    if(.not. num%onlyphbte) then
       call print_message("Decoupled electron BTE:")
       call print_message("---------------------")

       !gradT:
       call print_message("gradT field:")

       ! RTA:
       call print_message(" Calculating RTA electron kappa0 and sigmaS...")

       !  Allocate response function
       allocate(bt%el_response_T(el%nwv, el%numbands, 3))

       !  Read response function
       call readfile_response('RTA_I0_', bt%el_response_T, el%bandlist)

       !  Allocate spectral transport coefficients
       allocate(el_kappa0(el%numbands, 3, 3, num%el_en_num), &
            el_sigmaS(el%numbands, 3, 3, num%el_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(el, 'T', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, bt%el_response_T, el_en_grid, num%tetrahedra, sym, &
            el_kappa0, el_sigmaS)

       !  Write spectral electron kappa
       call write2file_spectral_tensor('RTA_el_kappa0_spectral_', el_kappa0, el%bandlist)

       !  Write spectral electron sigmaS
       call write2file_spectral_tensor('RTA_el_sigmaS_spectral_', el_sigmaS, el%bandlist)
       !------------------------------------------------------------------!

       ! Iterated:
       call print_message(" Calculating iterated electron kappa0 and sigmaS...")

       !  Read response function
       call readfile_response('nodrag_I0_', bt%el_response_T, el%bandlist)

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(el, 'T', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, bt%el_response_T, el_en_grid, num%tetrahedra, sym, &
            el_kappa0, el_sigmaS)

       !  Write spectral electron kappa
       call write2file_spectral_tensor('nodrag_iterated_el_kappa0_spectral_', el_kappa0, el%bandlist)

       !  Write spectral electron sigmaS
       call write2file_spectral_tensor('nodrag_iterated_el_sigmaS_spectral_', el_sigmaS, el%bandlist)

       !  Release memory
       deallocate(bt%el_response_T, el_kappa0, el_sigmaS)
       !------------------------------------------------------------------!

       !E:
       call print_message("E field:")

       ! RTA:
       call print_message(" Calculating RTA electron sigma and alpha/T...")

       !  Allocate response function
       allocate(bt%el_response_E(el%nwv, el%numbands, 3))

       !  Read response function
       call readfile_response('RTA_J0_', bt%el_response_E, el%bandlist)

       !  Allocate spectral transport coefficients
       allocate(el_sigma(el%numbands, 3, 3, num%el_en_num), &
            el_alphabyT(el%numbands, 3, 3, num%el_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(el, 'E', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, bt%el_response_E, el_en_grid, num%tetrahedra, sym, &
            el_alphabyT, el_sigma)

       !  Write spectral electron alpha/T
       call write2file_spectral_tensor('RTA_el_alphabyT_spectral_', el_alphabyT, el%bandlist)

       !  Write spectral electron sigma
       call write2file_spectral_tensor('RTA_el_sigma_spectral_', el_sigma, el%bandlist)
       !------------------------------------------------------------------!

       ! Iterated:
       call print_message(" Calculating iterated electron sigma and alpha/T...")

       !  Read response function
       call readfile_response('nodrag_J0_', bt%el_response_E, el%bandlist)

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(el, 'E', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, bt%el_response_E, el_en_grid, num%tetrahedra, sym, &
            el_alphabyT, el_sigma)

       !  Write spectral electron alpha/T
       call write2file_spectral_tensor('nodrag_iterated_el_alphabyT_spectral_', el_alphabyT, el%bandlist)

       !  Write spectral electron sigma
       call write2file_spectral_tensor('nodrag_iterated_el_sigma_spectral_', el_sigma, el%bandlist)

       !  Release memory
       deallocate(bt%el_response_E, el_alphabyT, el_sigma)
    end if

    !Decoupled phonon BTE
    if(.not. num%onlyebte) then
       call print_message("Decoupled phonon BTE:")
       call print_message("---------------------")

       !gradT:
       call print_message("gradT field:")

       ! RTA:
       call print_message(" Calculating RTA phonon kappa...")

       !  Allocate response function
       allocate(bt%ph_response_T(ph%nwv, ph%numbands, 3))

       !  Read response function
       call readfile_response('RTA_F0_', bt%ph_response_T)

       !  Allocate spectral transport coefficients
       allocate(ph_kappa(ph%numbands, 3, 3, num%ph_en_num), &
            dummy(ph%numbands, 3, 3, num%ph_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(ph, 'T', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
            crys%volume, bt%ph_response_T, ph_en_grid, num%tetrahedra, sym, ph_kappa, dummy)

       !  Write spectral phonon kappa
       call write2file_spectral_tensor('RTA_ph_kappa_spectral_', ph_kappa)
       !------------------------------------------------------------------!

       ! Iterated:
       call print_message(" Calculating iterated phonon kappa...")

       !  Read response function
       call readfile_response('nodrag_F0_', bt%ph_response_T)

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(ph, 'T', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
            crys%volume, bt%ph_response_T, ph_en_grid, num%tetrahedra, sym, ph_kappa, dummy)

       !  Write spectral phonon kappa
       call write2file_spectral_tensor('nodrag_iterated_ph_kappa_spectral_', ph_kappa)

       !  Release memory
       deallocate(bt%ph_response_T, ph_kappa, dummy)
    end if

    !Partially decoupled electron BTE
    if(num%drag) then
       call print_message("Partially decoupled electron BTE:")
       call print_message("---------------------")

       !gradT:
       call print_message("gradT field:")

       ! Iterated:
       call print_message(" Calculating iterated electron kappa0 and sigmaS...")

       !  Allocate response function
       allocate(bt%el_response_T(el%nwv, el%numbands, 3))

       !  Read response function
       call readfile_response('partdcpl_I0_', bt%el_response_T, el%bandlist)

       !  Allocate spectral transport coefficients
       allocate(el_kappa0(el%numbands, 3, 3, num%el_en_num), &
            el_sigmaS(el%numbands, 3, 3, num%el_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(el, 'T', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, bt%el_response_T, el_en_grid, num%tetrahedra, sym, &
            el_kappa0, el_sigmaS)

       !  Write spectral electron kappa
       call write2file_spectral_tensor('partdcpl_iterated_el_kappa0_spectral_', el_kappa0, el%bandlist)

       !  Write spectral electron sigmaS
       call write2file_spectral_tensor('partdcpl_iterated_el_sigmaS_spectral_', el_sigmaS, el%bandlist)

       !  Release memory
       deallocate(bt%el_response_T, el_kappa0, el_sigmaS)
       !------------------------------------------------------------------!

       !E:
       call print_message("E field:")

       ! Iterated:
       call print_message(" Calculating iterated electron sigma and alpha/T...")

       !  Allocate response function
       allocate(bt%el_response_E(el%nwv, el%numbands, 3))

       !  Read response function
       call readfile_response('partdcpl_J0_', bt%el_response_E, el%bandlist)

       !  Allocate spectral transport coefficients
       allocate(el_sigma(el%numbands, 3, 3, num%el_en_num), &
            el_alphabyT(el%numbands, 3, 3, num%el_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(el, 'E', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, bt%el_response_E, el_en_grid, num%tetrahedra, sym, &
            el_alphabyT, el_sigma)

       !  Write spectral electron alpha/T
       call write2file_spectral_tensor('partdcpl_iterated_el_alphabyT_spectral_', el_alphabyT, el%bandlist)

       !  Write spectral electron sigma
       call write2file_spectral_tensor('partdcpl_iterated_el_sigma_spectral_', el_sigma, el%bandlist)

       !  Release memory
       deallocate(bt%el_response_E, el_alphabyT, el_sigma)
    end if

    !Coupled electron BTE
    if(num%drag) then
       call print_message("Coupled electron BTE:")
       call print_message("---------------------")

       !gradT:
       call print_message("gradT field:")

       ! Iterated:
       call print_message(" Calculating iterated electron kappa0 and sigmaS...")

       !  Allocate response function
       allocate(bt%el_response_T(el%nwv, el%numbands, 3))

       !  Read response function
       call readfile_response('drag_I0_', bt%el_response_T, el%bandlist)

       !  Allocate spectral transport coefficients
       allocate(el_kappa0(el%numbands, 3, 3, num%el_en_num), &
            el_sigmaS(el%numbands, 3, 3, num%el_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(el, 'T', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, bt%el_response_T, el_en_grid, num%tetrahedra, sym, &
            el_kappa0, el_sigmaS)

       !  Write spectral electron kappa
       call write2file_spectral_tensor('drag_iterated_el_kappa0_spectral_', el_kappa0, el%bandlist)

       !  Write spectral electron sigmaS
       call write2file_spectral_tensor('drag_iterated_el_sigmaS_spectral_', el_sigmaS, el%bandlist)

       !  Release memory
       deallocate(bt%el_response_T, el_kappa0, el_sigmaS)
       !------------------------------------------------------------------!

       !E:
       call print_message("E field:")

       ! Iterated:
       call print_message(" Calculating iterated electron sigma and alpha/T...")

       !  Allocate response function
       allocate(bt%el_response_E(el%nwv, el%numbands, 3))

       !  Read response function
       call readfile_response('drag_J0_', bt%el_response_E, el%bandlist)

       !  Allocate spectral transport coefficients
       allocate(el_sigma(el%numbands, 3, 3, num%el_en_num), &
            el_alphabyT(el%numbands, 3, 3, num%el_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(el, 'E', crys%T, el%spindeg, el%chempot, el%ens, &
            el%vels, crys%volume, bt%el_response_E, el_en_grid, num%tetrahedra, sym, &
            el_alphabyT, el_sigma)

       !  Write spectral electron alpha/T
       call write2file_spectral_tensor('drag_iterated_el_alphabyT_spectral_', el_alphabyT, el%bandlist)

       !  Write spectral electron sigma
       call write2file_spectral_tensor('drag_iterated_el_sigma_spectral_', el_sigma, el%bandlist)

       !  Release memory
       deallocate(bt%el_response_E, el_alphabyT, el_sigma)
    end if

    !Coupled phonon BTE
    if(num%drag) then
       call print_message("Coupled phonon BTE:")
       call print_message("---------------------")

       !gradT:
       call print_message("gradT field:")

       ! Iterated:
       call print_message(" Calculating iterated phonon kappa...")

       !  Allocate response function
       allocate(bt%ph_response_T(ph%nwv, ph%numbands, 3))

       !  Read response function
       call readfile_response('drag_F0_', bt%ph_response_T)

       !  Allocate spectral transport coefficients
       allocate(ph_kappa(ph%numbands, 3, 3, num%ph_en_num), &
            dummy(ph%numbands, 3, 3, num%ph_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(ph, 'T', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
            crys%volume, bt%ph_response_T, ph_en_grid, num%tetrahedra, sym, ph_kappa, dummy)

       !  Write spectral phonon kappa
       call write2file_spectral_tensor('drag_iterated_ph_kappa_spectral_', ph_kappa)

       !  Release memory
       deallocate(bt%ph_response_T, ph_kappa)

       !E:
       call print_message("E field:")

       ! Iterated:
       call print_message(" Calculating iterated phonon alpha/T...")

       !  Allocate response function
       allocate(bt%ph_response_E(ph%nwv, ph%numbands, 3))

       !  Read response function
       call readfile_response('drag_G0_', bt%ph_response_E)

       !  Allocate spectral transport coefficients
       allocate(ph_alphabyT(ph%numbands, 3, 3, num%ph_en_num))

       !  Calculate spectral function
       call calculate_spectral_transport_coeff(ph, 'E', crys%T, 1_k8, 0.0_dp, ph%ens, ph%vels, &
            crys%volume, bt%ph_response_E, ph_en_grid, num%tetrahedra, sym, ph_alphabyT, dummy)

       !  Write spectral phonon kappa
       call write2file_spectral_tensor('drag_iterated_ph_alphabyT_spectral_', ph_alphabyT)

       !  Release memory
       deallocate(bt%ph_response_E, ph_alphabyT, dummy)
       !------------------------------------------------------------------!
    end if

    !Return to working directory
    call chdir(trim(adjustl(num%cwd)))

    sync all
  end subroutine post_process
end module bte_module
