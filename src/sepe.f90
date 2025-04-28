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
  use misc, only: Bose, mux_vector, binsearch, timer, subtitle, &
       print_message, distribute_points, demux_state
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
    real(r64), allocatable :: ph_drag_term_E(:, :, :), ph_coherence_term_E(:, :, :)
    type(timer) :: t
    integer :: it_el, s, it_ph_coh

    call t%start_timer('Iterative electron sector of SEPE')

    call print_message("Dragless electron transport:")
    call print_message("-----------------------------")

    !Restart with RTA solution
    self%el_response_E = self%el_field_term_E

    do it_ph_coh = 1, num%maxiter
       !call calculate_ph_coherenece(self%el_response_E, self%ph_coherence_E)

!!$       !Calculate the "sepe term": 1 + theta/n0
!!$       do s = 1, ph%numbands
!!$          sepe_term(:, s) = 1.0_r64 + &
!!$               self%ph_coherence(:, s)/Bose(ph%ens(:, s), crys%T)
!!$       end do

       !TODO start phonon iterator
       !call iterate_bte_ph(crys%T, num, crys, ph, el, self%ph_rta_rates_ibz, &
       !     self%ph_field_term_E, self%ph_response_E, self%el_response_E)
       
       do it_el = 1, num%maxiter
          !E field:
          call iterate_el_occupations_eqn(num, el, crys, self%el_rta_rates_ibz, self%el_field_term_E, &
               self%el_response_E, ph_drag_term_E, ph_coherence_term_E)
          
          !TODO Set up convergence criterion to exit iteration loop
       end do !electron iterator

       !TODO Set up convergence criterion to exit iteration loop
    end do !phonon coherence iterator

    sync all
  end subroutine dragless_el_eqn

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
    real(r64) :: tau_ibz
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
    allocate(coherence_ph_real(size(coherence_ph, 1), numbranches, 3))
    coherence_ph_real = real(coherence_ph)
    
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
          
          !Read Y from file
          if(allocated(U)) deallocate(U)
          call read_transition_probs_e(trim(adjustl(filepath_Y)), nprocs_phcoh, U)

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
             do iproc = 1, nprocs_phcoh
                response_ph_reduce(iq1_fbz, s1, :) = response_ph_reduce(iq1_fbz, s1, :) - &
                     el%spindeg*U(iproc)*coherence_ph_real(iq1_fbz, s1, :)
             end do
             
             !Iterate BTE
             response_ph_reduce(iq1_fbz, s1, :) = field_term(iq1_fbz, s1, :) + &
                  response_ph_reduce(iq1_fbz, s1, :)*tau_ibz          
          end do
       end do
    end if

    !Update the response function
    sync all
    call co_sum(response_ph_reduce)
    sync all
    response_ph = response_ph_reduce

    !Symmetrize response function
    do iq1_fbz = 1, nq
       response_ph(iq1_fbz,:,:)=transpose(&
            matmul(ph%symmetrizers(:,:,iq1_fbz),transpose(response_ph(iq1_fbz,:,:))))
    end do
  end subroutine iterate_ph_occupations_eqn
  
!!$  subroutine iterate_ph_coherence_eqn(T, num, crys, ph, el, rta_rates_ibz, &
!!$       field_term, response_ph, response_el, coherence_ph)
!!$    !! Subroutine to calculate the phonon coherence equation.
!!$    !! 
!!$    !! T Temperature in K
!!$    !! num Numerics object
!!$    !! crys Crystal object
!!$    !! ph Phonon object
!!$    !! el Electron object
!!$    !! rta_rates_ibz Phonon RTA scattering rates
!!$    !! field_term Phonon field coupling term
!!$    !! response_ph Phonon response function
!!$    !! response_el Electron response function
!!$    !! coherence_ph Phonon coherence term
!!$
!!$    type(phonon), intent(in) :: ph
!!$    type(electron), intent(in) :: el
!!$    type(numerics), intent(in) :: num
!!$    type(crystal), intent(in) :: crys
!!$    real(r64), intent(in) :: T, rta_rates_ibz(:, :), field_term(:, :, :)
!!$    real(r64), intent(in) :: response_el(:, :, :)
!!$    real(r64), intent(in) :: response_ph(:, :, :)
!!$    complex(r64), intent(inout) :: coherence_ph(:, :, :)
!!$
!!$    !Local variables
!!$    integer(i64) :: nstates_irred, chunk, istate1, numbranches, s1, &
!!$         iq1_ibz, ieq, iq1_sym, iq1_fbz, iproc, iq2, s2, iq3, s3, nq, &
!!$         num_active_images, numbands, ik, ikp, m, n, nprocs_phe, aux1, aux2, &
!!$         nprocs_3ph_plus, nprocs_3ph_minus, start, end, nprocs_phcoh
!!$    integer(i64), allocatable :: istate2_plus(:), istate3_plus(:), &
!!$         istate2_minus(:), istate3_minus(:), istate_el1(:), istate_el2(:)
!!$    real(r64) :: tau_ibz
!!$    real(r64), allocatable :: Y(:), U(:), response_ph_reduce(:, :, :), &
!!$         coherence_ph_real(:, :, :)
!!$    character(len = 1024) :: filepath_Wm, filepath_Wp, filepath_Y, filepath_U, tag
!!$
!!$    !Set output directory of transition probilities
!!$    write(tag, "(E9.3)") T
!!$    
!!$    !Number of electron bands
!!$    numbands = size(response_el(1,:,1))
!!$    
!!$    !Number of phonon branches
!!$    numbranches = size(rta_rates_ibz(1,:))
!!$
!!$    !Number of FBZ wave vectors
!!$    nq = size(field_term(:,1,1))
!!$    
!!$    !Total number of IBZ states
!!$    nstates_irred = size(rta_rates_ibz(:,1))*numbranches
!!$    
!!$    !Allocate and initialize response reduction array
!!$    allocate(response_ph_reduce(nq, numbranches, 3))
!!$    response_ph_reduce(:,:,:) = 0.0_r64
!!$    
!!$    !Divide phonon states among images
!!$    call distribute_points(nstates_irred, chunk, start, end, num_active_images)
!!$
!!$    !Only work with the active images
!!$    if(this_image() <= num_active_images) then
!!$       !TODO
!!$    end if
!!$
!!$    !Update the response function
!!$    call co_sum(response_ph_reduce)
!!$    response_ph = response_ph_reduce
!!$
!!$    !Symmetrize response function
!!$    do iq1_fbz = 1, nq
!!$       response_ph(iq1_fbz,:,:)=transpose(&
!!$            matmul(ph%symmetrizers(:,:,iq1_fbz),transpose(response_ph(iq1_fbz,:,:))))
!!$    end do    
!!$  end subroutine iterate_ph_coherence_eqn

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
    real(r64), allocatable :: Omegaplus(:), Omegaminus(:), &
         ph_coherence_term_reduce(:, :, :), coherence_ph_real(:, :, :)
    character(1024) :: filepath_Omegaminus, filepath_Omegaplus, tag
    
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
                   !and points. I note that response_ph(:, s, :) is not contiguous in memory.
                   iq2inter = mux_vector(fineq_indvec,el%wvmesh, 0_i64)
                   call interpolate_using_precomputed(&
                        coarse_mesh_corners(iq2inter,:), &
                        weights_coarse_mesh_corners(iq2inter,:),&
                        coherence_ph_real(:, s, :), HorP(:))
                else
                   !H(q) or P(q)
                   HorP(:) = coherence_ph_real(ph%equiv_map(ik_sym, iq), s, :)
                end if

                !(Note that below we use the fact that H and P are even in wave vector)
                ph_coherence_term_reduce(ik_fbz, m, :) = &
                     ph_coherence_term_reduce(ik_fbz, m, :) - &
                     HorP(:)*(Omegaplus(iproc) + Omegaminus(iproc))
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
  
end module SEPE_module
