module bte_module
  !! Module containing type and procedure related to solution of the
  !! Boltzmann transport equation (BTE).

  use params, only: dp
  use misc, only: print_message, write2file_rank2_real
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use phonon_module, only: phonon
  use interactions, only: calculate_ph_rta_rates

  implicit none

  private
  public bte

  type bte
     !! Data and procedures related to the BTE.

     real(dp), allocatable :: ph_rta_rates(:,:)
     !! Phonon RTA scattering rates.

   contains

     procedure :: solve_rta_ph
     
  end type bte

contains

  subroutine solve_rta_ph(bt, num, crys, ph)
    !! Obtain the RTA solution of the phonon BTE.

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
    call calculate_ph_rta_rates(ph, num, crys, bt%ph_rta_rates)

    !Write RTA scattering rates to file
    call write2file_rank2_real('ph.W_rta', bt%ph_rta_rates)
    
    !TODO Calculate RTA term F0 

    !TODO Calculate transport coefficients
  end subroutine solve_rta_ph
end module bte_module
