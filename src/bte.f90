module bte_module
  !! Module containing type and procedures related to the solution of the
  !! Boltzmann transport equation (BTE).

  use params, only: dp, k4, qe
  use misc, only: print_message, exit_with_message, write2file_rank2_real, &
       distribute_points
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use phonon_module, only: phonon
  use interactions, only: calculate_ph_rta_rates
  use bz_sums, only: calculate_transport_coeff

  implicit none

  private
  public bte

  type bte
     !! Data and procedures related to the BTE.

     real(dp), allocatable :: ph_rta_rates_ibz(:,:)
     !! Phonon RTA scattering rates on the IBZ.
     real(dp), allocatable :: ph_field_term(:,:,:)
     !! Phonon field coupling term on the FBZ.

   contains

     procedure :: solve_rta_ph
     
  end type bte

contains

  subroutine solve_rta_ph(bt, num, crys, ph)
    !! Subroutine to calculate the RTA solution of the phonon BTE.

    class(bte), intent(out) :: bt
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys
    type(phonon), intent(in) :: ph

    !Local variables
    character(len = 1024) :: tag, Tdir

    call print_message("Solving ph BTE in the RTA...")
    
    !Create output folder tagged by temperature and change into it
    write(tag, "(E9.3)") crys%T
    Tdir = trim(adjustl(num%cwd))//'/T'//trim(adjustl(tag))
    if(this_image() == 1) then
       call system('mkdir '//trim(adjustl(Tdir)))
    end if
    sync all
    
    !Calculate RTA scattering rates
    call calculate_ph_rta_rates(ph, num, crys, bt%ph_rta_rates_ibz)
    
    !Calculate RTA term F0
    call calculate_field_term('ph', 'T', ph%nequiv, ph%ibz2fbz_map, &
         crys%T, 0.0_dp, ph%ens, ph%vels, bt%ph_rta_rates_ibz, bt%ph_field_term)

    !Calculate transport coefficient
    call calculate_transport_coeff('ph', 'T', crys%T, 0.0_dp, ph%ens, ph%vels, &
         crys%volume, ph%qmesh, bt%ph_field_term)

    !Change to data output directory
    call chdir(trim(adjustl(Tdir)))

    !Write RTA scattering rates to file
    call write2file_rank2_real('ph.W_rta', bt%ph_rta_rates_ibz)
  end subroutine solve_rta_ph

  subroutine calculate_field_term(species, field, nequiv, ibz2fbz_map, &
       T, chempot, ens, vels, rta_rates_ibz, field_term)
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

    character(len = 2), intent(in) :: species
    character(len = 1), intent(in) :: field
    integer(k4), intent(in) :: nequiv(:), ibz2fbz_map(:,:,:)
    real(dp), intent(in) :: T, chempot, ens(:,:), vels(:,:,:), rta_rates_ibz(:,:)
    real(dp), allocatable, intent(out) :: field_term(:,:,:)

    !Local variables
    integer(k4) :: ik_ibz, ik_fbz, ieq, ib, nk_ibz, nk, nbands, im, chunk
    integer(k4), allocatable :: start[:], end[:]
    real(dp), allocatable :: field_term_reduce(:,:,:)[:]
    real(dp) :: A
    logical :: trivial_case

    !Set constant depending on species and field type
    if(species == 'ph') then
       A = 1.0_dp/T

       if(chempot /= 0.0_dp) then
          call exit_with_message("Phonon chemical potential non-zero in calculate_field_term. Exiting.")
       end if
    else if(species == 'el') then
       if(field == 'T') then
          A = 1.0_dp/T
       else if(field == 'E') then
          A = qe
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
       call distribute_points(nk_ibz, chunk, start, end)

       !Allocate and initialize field term coarrays
       allocate(field_term_reduce(nk, nbands, 3)[*])
       field_term_reduce(:,:,:) = 0.0_dp
       
       do ik_ibz = start, end
          do ieq = 1, nequiv(ik_ibz)
             ik_fbz = ibz2fbz_map(ieq, ik_ibz, 2)
             do ib = 1, nbands
                if(rta_rates_ibz(ik_ibz, ib) /= 0.0_dp) then
                   field_term_reduce(ik_fbz, ib, :) = A*vels(ik_fbz, ib, :)*&
                        (ens(ik_fbz, ib) - chempot)/rta_rates_ibz(ik_ibz, ib)
                end if
             end do
          end do
       end do
       
       sync all

       !Reduce field term coarrays
       do im = 1, num_images()
          field_term(:,:,:) = field_term(:,:,:) + field_term_reduce(:,:,:)[im]
       end do
       sync all
    end if
  end subroutine calculate_field_term
end module bte_module
