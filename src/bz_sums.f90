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

module bz_sums
  !! Module containing the procedures to do Brillouin zone sums.

  use params, only: dp, k8, kB, qe, pi, hbar_eVps, perm0
  use misc, only: exit_with_message, print_message, write2file_rank2_real, &
       distribute_points, Bose, Fermi
  use phonon_module, only: phonon
  use electron_module, only: electron
  use crystal_module, only: crystal
  use delta, only: delta_fn_tetra, delta_fn_triang
  use symmetry_module, only: symmetry, symmetrize_3x3_tensor

  implicit none

  public calculate_dos, calculate_transport_coeff, calculate_qTF
  private calculate_el_dos, calculate_ph_dos_iso

  interface calculate_dos
     module procedure :: calculate_el_dos, calculate_ph_dos_iso
  end interface calculate_dos
  
contains
  
  subroutine calculate_qTF(crys, el)
    !! Calculate Thomas-Fermi screening wavevector in the simple electron-gas model.
    ! qTF**2 = spindeg*e^2*beta/nptq/vol_pcell/perm0/kappainf*Sum_{BZ}f0_{k}(1-f0_{k})

    type(crystal), intent(inout) :: crys
    type(electron), intent(in) :: el

    !Local variables
    real(dp) :: beta, fFD
    integer(k8) :: ib, ik

    beta = 1.0_dp/kB/crys%T/qe !1/J
    crys%qTF=0.d0

    if(crys%polar) then
       call print_message("Calculating Thomas-Fermi screening...")
       
       do ib = 1, el%numbands
          do ik = 1, el%nk
             fFD = Fermi(el%ens(ik, ib), el%chempot, crys%T)
             crys%qTF = crys%qTF + fFD*(1.0_dp - fFD)
          end do
       end do
       crys%qTF = sqrt(1.0e9_dp*crys%qTF*el%spindeg*beta*qe**2/product(el%kmesh)&
            /crys%volume/perm0/crys%epsiloninf) !nm^-1

       if(this_image() == 1) then
          write(*, "(A, 1E16.8, A)") ' Thomas-Fermi screening wave length = ', crys%qTF, ' 1/nm'
       end if
    end if
  end subroutine calculate_qTF
  
  subroutine calculate_el_dos(el, usetetra)
    !! Calculate the density of states (DOS) in units of 1/energy. 
    !! The DOS will be evaluates on the IBZ mesh energies.
    !!
    !! el Electron data type
    !! usetetra Use the tetrahedron method for delta functions?

    type(electron), intent(inout) :: el
    logical, intent(in) :: usetetra
    
    !Local variables
    integer(k8) :: ik, ib, ikp, ibp, im, chunk, counter, num_active_images
    integer(k8), allocatable :: start[:], end[:]
    real(dp) :: e, delta
    real(dp), allocatable :: dos_chunk(:,:)[:]

    call print_message("Calculating electron density of states...")
    
    !Allocate start and end coarrays
    allocate(start[*], end[*])
    
    !Divide wave vectors among images
    call distribute_points(el%nk_irred, chunk, start, end, num_active_images)

    !Allocate small work variable chunk for each image
    allocate(dos_chunk(chunk, el%numbands)[*])
    
    !Allocate dos
    allocate(el%dos(el%nk_irred, el%numbands))

    !Initialize dos arrays
    el%dos(:,:) = 0.0_dp
    dos_chunk(:,:) = 0.0_dp

    counter = 0
    do ik = start, end !Run over IBZ wave vectors
       !Increase counter
       counter = counter + 1
       do ib = 1, el%numbands !Run over wave vectors   
          !Grab sample energy from the IBZ
          e = el%ens_irred(ik, ib) 

          do ikp = 1, el%nk !Sum over FBZ wave vectors
             do ibp = 1, el%numbands !Sum over wave vectors
                if(usetetra) then
                   !Evaluate delta[E(iq,ib) - E(iq',ib')]
                   delta = delta_fn_tetra(e, ikp, ibp, el%kmesh, el%tetramap, &
                        el%tetracount, el%tetra_evals)
                else
                   delta = delta_fn_triang(e, ikp, ibp, el%kmesh, el%triangmap, &
                        el%triangcount, el%triang_evals)
                end if
                !Sum over delta function
                dos_chunk(counter, ib) = dos_chunk(counter, ib) + delta
             end do
          end do
       end do
    end do
    !Multiply with spin degeneracy factor
    dos_chunk(:,:) = el%spindeg*dos_chunk(:,:)
    sync all
    
    !Collect dos_chunks into dos
    do im = 1, num_active_images
       el%dos(start[im]:end[im], :) = dos_chunk(:,:)[im]
    end do
    sync all

    !Write dos to file
    call write2file_rank2_real(el%prefix // '.dos', el%dos)

    sync all
  end subroutine calculate_el_dos

  subroutine calculate_ph_dos_iso(ph, usetetra, gfactors, atomtypes, W_phiso, phiso)
    !! Calculate the phonon density of states (DOS) in units of 1/energy and,
    !! optionally, the phonon-isotope scattering rates.
    !!
    !! The DOS and isotope scattering rates will be evaluates on the IBZ mesh energies.
    !!
    !! ph Phonon data type
    !! usetetra Use the tetrahedron method for delta functions?

    type(phonon), intent(inout) :: ph
    logical, intent(in) :: usetetra, phiso
    real(dp), intent(in) :: gfactors(:)
    integer(k8), intent(in) :: atomtypes(:)
    real(dp), intent(out), allocatable :: W_phiso(:,:)
    
    !Local variables
    integer(k8) :: iq, ib, iqp, ibp, im, chunk, counter, num_active_images, &
         pol, a, numatoms
    integer(k8), allocatable :: start[:], end[:]
    real(dp) :: e, delta, aux
    real(dp), allocatable :: dos_chunk(:,:)[:], W_phiso_chunk(:,:)[:]
    
    call print_message("Calculating phonon density of states and (if needed) isotope scattering...")

    !Number of basis atoms
    numatoms = size(atomtypes)
    
    !Allocate start and end coarrays
    allocate(start[*], end[*])
    
    !Divide wave vectors among images
    call distribute_points(ph%nq_irred, chunk, start, end, num_active_images)
    
    !Allocate small work variable chunk for each image
    allocate(dos_chunk(chunk, ph%numbranches)[*])
    if(phiso) allocate(W_phiso_chunk(chunk, ph%numbranches)[*])
    
    !Allocate dos and W_phiso
    allocate(ph%dos(ph%nq_irred, ph%numbranches))
    allocate(W_phiso(ph%nq_irred, ph%numbranches))

    !Initialize arrays and coarrays
    ph%dos(:,:) = 0.0_dp
    dos_chunk(:,:) = 0.0_dp
    W_phiso(:,:) = 0.0_dp
    if(phiso) W_phiso_chunk(:,:) = 0.0_dp

    counter = 0
    do iq = start, end !Run over IBZ wave vectors
       !Increase counter
       counter = counter + 1
       do ib = 1, ph%numbranches !Run over wave vectors   
          !Grab sample energy from the IBZ
          e = ph%ens(ph%indexlist_irred(iq), ib) 

          do iqp = 1, ph%nq !Sum over FBZ wave vectors
             do ibp = 1, ph%numbranches !Sum over wave vectors
                !Evaluate delta[E(iq,ib) - E(iq',ib')]
                if(usetetra) then
                   delta = delta_fn_tetra(e, iqp, ibp, ph%qmesh, ph%tetramap, &
                        ph%tetracount, ph%tetra_evals)
                else
                   delta = delta_fn_triang(e, iqp, ibp, ph%qmesh, ph%triangmap, &
                        ph%triangcount, ph%triang_evals)
                end if
                !Sum over delta function
                dos_chunk(counter, ib) = dos_chunk(counter, ib) + delta

                if(phiso) then
                   !Calculate phonon-isotope scattering in the Tamura model
                   do a = 1, numatoms
                      pol = (a - 1)*3
                      aux = (abs(dot_product(&
                           ph%evecs(ph%indexlist_irred(iq), ib, pol + 1 : pol + 3), &
                           ph%evecs(iqp, ibp, pol + 1 : pol + 3))))**2
                      W_phiso_chunk(counter, ib) = W_phiso_chunk(counter, ib) + &
                           delta*aux*gfactors(atomtypes(a))*e**2
                   end do
                end if
             end do
          end do
       end do
    end do
    if(phiso) W_phiso_chunk = W_phiso_chunk*0.5_dp*pi/hbar_eVps !THz
    sync all
    
    !Collect chunks into full array
    do im = 1, num_active_images
       ph%dos(start[im]:end[im], :) = dos_chunk(:,:)[im]
       if(phiso) W_phiso(start[im]:end[im], :) = W_phiso_chunk(:,:)[im]
    end do
    sync all

    !Write to file
    call write2file_rank2_real(ph%prefix // '.dos', ph%dos)
    call write2file_rank2_real(ph%prefix // '.W_rta_phiso', W_phiso)

    sync all
  end subroutine calculate_ph_dos_iso

  subroutine calculate_transport_coeff(species, field, T, deg, chempot, ens, vels, &
       volume, mesh, response, sym, trans_coeff_hc, trans_coeff_cc)
    !! Subroutine to calculate transport coefficients.
    !!
    !! species Type of particle
    !! field Type of field
    !! T Temperature in K
    !! deg Degeneracy
    !! chempot Chemical potential in eV
    !! ens FBZ energies in eV
    !! vels FBZ velocities in Km/s
    !! volume Primitive cell volume in nm^3
    !! mesh Wave vector grid
    !! response FBZ response function
    !! trans_coeff_hc Heat current coefficient
    !! trans_coeff_cc Charge current coefficient

    character(len = 2), intent(in) :: species
    character(len = 1), intent(in) :: field
    integer(k8), intent(in) :: mesh(3), deg
    real(dp), intent(in) :: T, chempot, ens(:,:), vels(:,:,:), volume, response(:,:,:)
    type(symmetry), intent(in) :: sym
    real(dp), intent(out) :: trans_coeff_hc(:,:,:), trans_coeff_cc(:,:,:)
    ! Above, h(c)c = heat(charge) current
    
    !Local variables
    integer(k8) :: ik, ib, icart, nk, nbands, pow_hc, pow_cc
    real(dp) :: dist_factor, e, v, fac, A_hc, A_cc
    
    nk = size(ens(:,1))
    nbands = size(ens(1,:))

    !Common multiplicative factor
    fac = 1.0e21/kB/T/volume/product(mesh) 
    
    !Do checks related to particle and field type
    if(species == 'ph') then
       if(chempot /= 0.0_dp) then
          call exit_with_message(&
               "Phonon chemical potential non-zero in calculate_transport_coefficient. Exiting.")
       end if
       if(field == 'T') then
          A_hc = qe*fac
          pow_hc = 1
          A_cc = 0.0_dp
          pow_cc = 0
       else if(field == 'E') then
          A_hc = -fac
          pow_hc = 1
          A_cc = 0.0_dp
          pow_cc = 0
       else
          call exit_with_message("Unknown field type in calculate_transport_coefficient. Exiting.")
       end if
    else if(species == 'el') then       
       if(field == 'T') then
          A_cc = -deg*qe*fac
          pow_cc = 0
          A_hc = deg*qe*fac
          pow_hc = 1
       else if(field == 'E') then
          A_cc = deg*fac
          pow_cc = 0
          A_hc = -A_cc
          pow_hc = 1
       else
          call exit_with_message("Unknown field type in calculate_transport_coefficient. Exiting.")
       end if
    else
       call exit_with_message("Unknown particle species in calculate_transport_coefficient. Exiting.")
    end if
    
    trans_coeff_hc = 0.0_dp
    trans_coeff_cc = 0.0_dp
    do ik = 1, nk
       do ib = 1, nbands
          e = ens(ik, ib)
          if(species == 'ph') then
             if(e == 0.0_dp) cycle !Ignore zero energies phonons
             dist_factor = Bose(e, T)
             dist_factor = dist_factor*(1.0_dp + dist_factor)
          else
             dist_factor = Fermi(e, chempot, T)
             dist_factor = dist_factor*(1.0_dp - dist_factor)
          end if
          do icart = 1, 3
             v = vels(ik, ib, icart)
             trans_coeff_hc(ib, icart, :) = trans_coeff_hc(ib, icart, :) + &
                  (e - chempot)**pow_hc*dist_factor*v*response(ik, ib, :)
             trans_coeff_cc(ib, icart, :) = trans_coeff_cc(ib, icart, :) + &
                  (e - chempot)**pow_cc*dist_factor*v*response(ik, ib, :)
          end do
       end do
    end do
    !Units:
    ! W/m/K for thermal conductivity
    ! 1/Omega/m for charge conductivity
    ! V/K for thermopower
    ! A/m/K for alpha/T
    trans_coeff_hc = A_hc*trans_coeff_hc
    trans_coeff_cc = A_cc*trans_coeff_cc

    !Symmetrize transport tensor
    do ib = 1, nbands
       call symmetrize_3x3_tensor(trans_coeff_hc(ib, :, :), sym%crotations)
       call symmetrize_3x3_tensor(trans_coeff_cc(ib, :, :), sym%crotations)
    end do
  end subroutine calculate_transport_coeff
end module bz_sums
