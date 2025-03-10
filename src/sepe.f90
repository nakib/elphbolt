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
  use interactions, only: read_transition_probs_e

  implicit none
  
  private
  public sepe

  type sepe
     !! Data and procedures related to the BTE.

!!$     real(r64), allocatable :: ph_rta_rates_phe_ibz(:,:)
!!$     !! Phonon RTA scattering rates on the IBZ due to ph-e interactions.
!!$     real(r64), allocatable :: ph_rta_rates_ibz(:,:)
!!$     !! Phonon RTA scattering rates on the IBZ.
!!$     real(r64), allocatable :: ph_field_term_T(:,:,:)
!!$     !! Phonon field coupling term for gradT field on the FBZ.
!!$     real(r64), allocatable :: ph_response_T(:,:,:)
!!$     !! Phonon response function for gradT field on the FBZ.
!!$     real(r64), allocatable :: ph_field_term_E(:,:,:)
!!$     !! Phonon field coupling term for E field on the FBZ.
!!$     real(r64), allocatable :: ph_response_E(:,:,:)
!!$     !! Phonon response function for E field on the FBZ.
     complex(r64), allocatable :: ph_coherence_T(:, :)
     !! Phonon coherence term for gradT field on the FBZ
     complex(r64), allocatable :: ph_coherence_E(:, :)
     !! Phonon coherence term for E field on the FBZ
     
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
    !real(r64) :: sepe_term(ph%nwv, ph%numbands)
    integer :: it_el, s

    call t%start_timer('Iterative electron sector of SEPE')

    call print_message("Dragless electron transport:")
    call print_message("-----------------------------")

    !Restart with RTA solution
    self%el_response_E = self%el_field_term_E

    do it_coh = 1, num%maxiter
       call calculate_ph_coherenece(self%el_response_E, self%ph_coherence)

!!$       !Calculate the "sepe term": 1 + theta/n0
!!$       do s = 1, ph%numbands
!!$          sepe_term(:, s) = 1.0_r64 + &
!!$               self%ph_coherence(:, s)/Bose(ph%ens(:, s), crys%T)
!!$       end do

       do it_el = 1, num%maxiter
          !E field:
          call iterate_el_eqn(num, el, crys, &
               self%el_rta_rates_ibz, self%el_field_term_E, sepe_term, self%el_response_E)

          !TODO Set up convergence criterion to exit iteration loop
       end do !electron iterator

       !TODO Set up convergence criterion to exit iteration loop
    end do !phonon coherence iterator

    sync all
  end subroutine dragless_el_eqn

  subroutine iterate_el_eqn(num, el, crys, rta_rates_ibz, field_term, &
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
         aux, aux2, aux3, aux4, start, end, nprocs_echimp, neg_ik_fbz, nprocs_phcoh
    integer(i64), allocatable :: istate_el(:), istate_ph(:), istate_el_echimp(:), &
         istate_el_phcoh(:), istate_ph_phcoh(:)

    real(r64) :: tau_ibz
    real(r64), allocatable :: Xphplus(:), Xphminus(:), Xchimp(:), &
         response_el_reduce(:, :, :), Omegaplus(:), Omegaminus(:)
    character(1024) :: filepath_Xphminus, filepath_Xphplus, filepath_Xechimp, tag, &
         filepath_Omegaminus, filepath_Omegaplus

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
          
          !The e-ph (coherence) bit:
          
          !Set Omega+ filename
          write(tag, '(I9)') istate
          filepath_Omegaplus = trim(adjustl(num%Xdir))//'/Omegaplus.istate'//trim(adjustl(tag))

          !Read Omega+ from file
          call read_transition_probs_e(trim(adjustl(filepath_Omegaplus)), nprocs_phcoh, Omegaplus, &
               istate_el_phcoh, istate_ph_phcoh)

          !Set Omega- filename
          write(tag, '(I9)') istate
          filepath_Omegaminus = trim(adjustl(num%Xdir))//'/Omegaminus.istate'//trim(adjustl(tag))

          !Read Omega- from file
          call read_transition_probs_e(trim(adjustl(filepath_Omegaminus)), nprocs_phcoh, Omegaminus)

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

             !The e-ph coherence bit:
             
             !Sum over scattering processes
             do iproc = 1, nprocs_phcoh
                !Grab the final electron
                call demux_state(istate_el_phcoh(iproc), numbands, n, ikp)

                !Find image of final electron wave vector due to the current symmetry
                call binsearch(el%indexlist, el%equiv_map(ik_sym, ikp), aux)

                !(Note that below we use the fact that H and P are even in wave vector)
                response_el_reduce(ik_fbz, m, :) = response_el_reduce(ik_fbz, m, :) + &
                     response_el(aux, n, :)*(-Omegaplus(iproc) - Omegaminus(iproc))
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
  end subroutine iterate_el_eqn

  subroutine calculate_phonon_drag(num, el, ph, idc, widc, sym, rta_rates_ibz, response_ph, &
       ph_drag_term)
    !! Subroutine to calculate the phonon drag term.
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
         fineq_indvec(3), start, end, iq2inter
    integer(i64), allocatable :: istate_el(:), istate_ph(:)
    real(r64) :: tau_ibz, ForG(3)
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
    ph_drag_term_reduce(:,:,:) = 0.0_r64

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
                   fineq_indvec = modulo( &
                        nint(matmul(sym%qrotations(:, :, ik_sym), fineq_indvec)), el%wvmesh)

                   !Interpolate response function on this wave vector using precomputed tabulated weights
                   !and points. I note that response_ph(:, s, :) is not contiguous in memory.
                   iq2inter = mux_vector(fineq_indvec,el%wvmesh, 0_i64)
                   call interpolate_using_precomputed(idc(iq2inter,:), widc(iq2inter,:),&
                        response_ph(:, s, :), ForG(:))
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
    end if

    !Reduce from all images
    call co_sum(ph_drag_term_reduce)
    ph_drag_term = ph_drag_term_reduce
  end subroutine calculate_phonon_drag
  
!!$  subroutine calculate_ph_coherenece(ph, ph_coherence_term)
!!$    !! Computes the phonon coherence.
!!$    !!
!!$    !! ph_coherence Phonon coherence
!!$    
!!$    real(r64), intent(inout) :: ph_coherence_term(:, :)
!!$    
!!$    !Number of electron bands
!!$    numbands = el%numbands
!!$
!!$    !Number of in-window FBZ wave vectors
!!$    nk = el%nwv
!!$
!!$    !Total number of IBZ states
!!$    nstates_irred = el%nwv_irred*numbands
!!$    
!!$    !Number of phonon branches
!!$    numbranches = ph%numbands
!!$
!!$    !Allocate and initialize response reduction array
!!$    allocate(ph_coherence_term_reduce(nk, numbands, 3))
!!$    ph_drag_term_reduce(:,:,:) = 0.0_r64
!!$
!!$    !Divide electron states among images
!!$    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
!!$
!!$    !Only work with the active images
!!$    if(this_image() <= num_active_images) then
!!$       !Run over electron IBZ states
!!$       do istate = start, end
!!$          !Demux state index into band (m) and wave vector (ik_ibz) indices
!!$          call demux_state(istate, numbands, m, ik_ibz)
!!$
!!$          !Apply energy window to initial (IBZ blocks) electron
!!$          if(abs(el%ens_irred(ik_ibz, m) - el%enref) > el%fsthick) cycle
!!$
!!$          !RTA lifetime
!!$          tau_ibz = 0.0_r64
!!$          if(rta_rates_ibz(ik_ibz, m) /= 0.0_r64) then
!!$             tau_ibz = 1.0_r64/rta_rates_ibz(ik_ibz, m)
!!$          end if
!!$
!!$          !Set X+ filename
!!$          write(tag, '(I9)') istate
!!$          filepath_Xplus = trim(adjustl(num%Xdir))//'/Xplus.istate'//trim(adjustl(tag))
!!$
!!$          !Read X+ from file
!!$          call read_transition_probs_e(trim(adjustl(filepath_Xplus)), nprocs, Xplus, &
!!$               istate_el, istate_ph)
!!$
!!$          !Set X- filename
!!$          write(tag, '(I9)') istate
!!$          filepath_Xminus = trim(adjustl(num%Xdir))//'/Xminus.istate'//trim(adjustl(tag))
!!$
!!$          !Read X- from file
!!$          call read_transition_probs_e(trim(adjustl(filepath_Xminus)), nprocs, Xminus)
!!$
!!$          !Sum over the number of equivalent k-points of the IBZ point
!!$          do ieq = 1, el%nequiv(ik_ibz)
!!$             ik_sym = el%ibz2fbz_map(ieq, ik_ibz, 1) !symmetry
!!$             call binsearch(el%indexlist, el%ibz2fbz_map(ieq, ik_ibz, 2), ik_fbz)
!!$
!!$             !Sum over scattering processes
!!$             do iproc = 1, nprocs
!!$                if(istate_ph(iproc) < 0) then !This phonon is on the (fine) electron mesh
!!$                   call demux_state(-istate_ph(iproc), numbranches, s, iq)
!!$                   iq = -iq !Keep the negative tag
!!$                else !This phonon is on the phonon mesh
!!$                   call demux_state(istate_ph(iproc), numbranches, s, iq)
!!$                end if
!!$
!!$                !Drag contribution:             
!!$                if(iq < 0) then !Need to interpolate on this point
!!$                   !Calculate the fine mesh wave vector, 0-based index vector
!!$                   call demux_vector(-iq, fineq_indvec, el%wvmesh, 0_i64)
!!$
!!$                   !Find image of phonon wave vector due to the current symmetry
!!$                   fineq_indvec = modulo( &
!!$                        nint(matmul(sym%qrotations(:, :, ik_sym), fineq_indvec)), el%wvmesh)
!!$
!!$                   !Interpolate response function on this wave vector using precomputed tabulated weights
!!$                   !and points. I note that response_ph(:, s, :) is not contiguous in memory.
!!$                   iq2inter = mux_vector(fineq_indvec,el%wvmesh, 0_i64)
!!$                   call interpolate_using_precomputed(idc(iq2inter,:), widc(iq2inter,:),&
!!$                        response_ph(:, s, :), ForG(:))
!!$                else
!!$                   !F(q) or G(q)
!!$                   ForG(:) = response_ph(ph%equiv_map(ik_sym, iq), s, :)
!!$                end if
!!$                !Here we use the fact that F(-q) = -F(q) and G(-q) = -G(q)
!!$                ph_drag_term_reduce(ik_fbz, m, :) = ph_drag_term_reduce(ik_fbz, m, :) - &
!!$                     ForG(:)*(Xplus(iproc) + Xminus(iproc))
!!$             end do
!!$
!!$             !Multiply life time factor 
!!$             ph_drag_term_reduce(ik_fbz, m, :) = ph_drag_term_reduce(ik_fbz, m, :)*tau_ibz
!!$          end do
!!$       end do
!!$    end if
!!$
!!$    !Reduce from all images
!!$    call co_sum(ph_drag_term_reduce)
!!$    ph_drag_term = ph_drag_term_reduce
!!$  end subroutine calculate_ph_coherenece
end module SEPE_module
