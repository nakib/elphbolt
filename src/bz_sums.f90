module bz_sums
  !! Module containing the procedures to do Brillouin zone sums.

  use params, only: dp, k8, kB, qe
  use misc, only: exit_with_message, print_message, write2file_rank2_real, &
       distribute_points, Bose, Fermi
  use phonon_module, only: phonon
  use electron_module, only: electron
  use delta, only: delta_fn_tetra

  implicit none

  public calculate_dos, calculate_transport_coeff, calculate_chempot
  private calculate_el_dos, calculate_ph_dos_iso

  interface calculate_dos
     module procedure :: calculate_el_dos, calculate_ph_dos_iso
  end interface calculate_dos
  
contains

  subroutine calculate_chempot(el, T, vol)
    !! Subroutine to calculate the chemical potential for a
    !! given carrier concentration.

    type(electron), intent(inout) :: el
    real(dp), intent(in) :: T, vol

    !Local variables
    real(dp) :: a, b, aux, const, absconc, signconc, mu, thresh
    integer(k8) :: ib, ik, it, ngrid, maxiter

    call print_message("Calculating chemical potential...")

    !Total number of points in full mesh
    ngrid = product(el%kmesh)

    !Normalization and units factor
    const = el%spindeg/dble(ngrid)/vol/(1.0e-21_dp)

    !Maximum number of iterations
    maxiter = 1000

    !Convergence threshold
    thresh = 1.0e-12_dp

    !Absolute value and sign of concentration
    absconc = abs(el%conc)
    signconc = sign(1.0_dp, el%conc)

    a = el%enref - 20.0_dp !guess lower bound
    b = el%enref + 20.0_dp !guess upper bound
    do it = 1, maxiter
       mu = 0.5_dp*(a + b)
       aux = 0.0_dp
       do ib = 1, el%numbands
          do ik = 1, el%nk
             aux = aux + Fermi(el%ens(ik, ib), mu, T)
          end do
       end do
       aux = aux*const !cm^-3
       if(abs(aux - absconc)/absconc < thresh) then
          exit
       else if(aux < absconc) then
          a = mu
       else
          b = mu
       end if
    end do
    el%chempot = mu

    if(abs(aux - absconc)/absconc > thresh) then
       call exit_with_message(&
            "Could not converge to correct chemical potential. Exiting.")
    end if

    if(this_image() == 1) then
       print*, 'Carrier concentration =', signconc*aux, ' cm^-3'
       print*, 'Chemical potential = ', el%chempot, ' eV'
    end if
  end subroutine calculate_chempot
  
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

                   !Sum over delta function
                   dos_chunk(counter, ib) = dos_chunk(counter, ib) + delta

                   !
                   !TODO need to implement Gaussian broadening
                   !
                end if
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

  subroutine calculate_ph_dos_iso(ph, usetetra)
    !! Calculate the phonon density of states (DOS) in units of 1/energy and,
    !! optionally, the phonon-isotope scattering rates.
    !!
    !! The DOS and isotopr scattering rates will be evaluates on the IBZ mesh energies.
    !!
    !! ph Phonon data type
    !! usetetra Use the tetrahedron method for delta functions?

    type(phonon), intent(inout) :: ph
    logical, intent(in) :: usetetra
    
    !Local variables
    integer(k8) :: iq, ib, iqp, ibp, im, chunk, counter, num_active_images
    integer(k8), allocatable :: start[:], end[:]
    real(dp) :: e, delta
    real(dp), allocatable :: dos_chunk(:,:)[:]

    call print_message("Calculating phonon density of states...")
    
    !Allocate start and end coarrays
    allocate(start[*], end[*])
    
    !Divide wave vectors among images
    call distribute_points(ph%nq_irred, chunk, start, end, num_active_images)
    
    !Allocate small work variable chunk for each image
    allocate(dos_chunk(chunk, ph%numbranches)[*])
    
    !Allocate dos
    allocate(ph%dos(ph%nq_irred, ph%numbranches))

    !Initialize dos arrays
    ph%dos(:,:) = 0.0_dp
    dos_chunk(:,:) = 0.0_dp

    counter = 0
    do iq = start, end !Run over IBZ wave vectors
       !Increase counter
       counter = counter + 1
       do ib = 1, ph%numbranches !Run over wave vectors   
          !Grab sample energy from the IBZ
          e = ph%ens(ph%indexlist_irred(iq), ib) 
          
          do iqp = 1, ph%nq !Sum over FBZ wave vectors
             do ibp = 1, ph%numbranches !Sum over wave vectors
                if(usetetra) then
                   !Evaluate delta[E(iq,ib) - E(iq',ib')]
                   delta = delta_fn_tetra(e, iqp, ibp, ph%qmesh, ph%tetramap, &
                        ph%tetracount, ph%tetra_evals)

                   !Sum over delta function
                   dos_chunk(counter, ib) = dos_chunk(counter, ib) + delta
                   
                   !
                   !TODO need to implement Gaussian broadening
                   !
                end if
             end do
          end do
       end do
    end do

    sync all
    
    !Collect dos_chunks into dos
    do im = 1, num_active_images
       ph%dos(start[im]:end[im], :) = dos_chunk(:,:)[im]
    end do
    sync all

    !Write dos to file
    call write2file_rank2_real(ph%prefix // '.dos', ph%dos)

    sync all
  end subroutine calculate_ph_dos_iso

  subroutine calculate_transport_coeff(species, field, T, deg, chempot, ens, vels, &
       volume, mesh, response)
    !! Subroutine to calculate transport coefficients.
    !!
    !! species Type of particle
    !! field Type of field
    !! T Temperature in K
    !! deg Degeneracy
    !! chempot Chemical potential in eV
    !! ens FBZ energies
    !! vels FBZ velocities
    !! volume Primitive cell volume
    !! mesh Wave vector grid
    !! response FBZ response function

    character(len = 2), intent(in) :: species
    character(len = 1), intent(in) :: field
    integer(k8), intent(in) :: mesh(3), deg
    real(dp), intent(in) :: T, chempot, ens(:,:), vels(:,:,:), volume, response(:,:,:)

    !Local variables
    integer(k8) :: ik, ib, icart, nk, nbands, pow
    real(dp) :: dist_factor, e, v, A, trans_coeff(3,3)
    
    nk = size(ens(:,1))
    nbands = size(ens(1,:))
    
    !Do checks related to particle and field type
    if(species == 'ph') then
       if(chempot /= 0.0_dp) then
          call exit_with_message(&
               "Phonon chemical potential non-zero in calculate_transport_coefficient. Exiting.")
       end if
       if(field == 'T') then !thermal conductivity
          A = qe*1.0e21/kB/T/volume/product(mesh)
          pow = 1
       else if(field == 'E') then
          !A = 
          !pow = 
       else
          call exit_with_message("Unknown field type in calculate_transport_coefficient. Exiting.")
       end if
    else if(species == 'el') then
       if(field == 'T') then
          !A = 
          !pow = 
       else if(field == 'E') then
          A = deg*1.0e21/kB/T/volume/product(mesh)
          pow = 0
       else
          call exit_with_message("Unknown field type in calculate_transport_coefficient. Exiting.")
       end if
    else
       call exit_with_message("Unknown particle species in calculate_transport_coefficient. Exiting.")
    end if
    
    trans_coeff(:,:) = 0.0_dp
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
             trans_coeff(icart, :) = trans_coeff(icart, :) + &
                  (e - chempot)**pow*dist_factor*v*response(ik, ib, :)
          end do
       end do
    end do
    !Units:
    ! W/m/K for thermal conductivity
    ! 1/Omega/m for charge conductivity
    ! V/K for thermopower
    trans_coeff = A*trans_coeff
    sync all

    if(this_image() == 1) print*, trans_coeff
  end subroutine calculate_transport_coeff
end module bz_sums
