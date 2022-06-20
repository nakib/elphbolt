! Copyright (C) 2022- Nakib Haider Protik <nakib.haider.protik@gmail.com>
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

module MigEl_sc_module
  !! Module containing types and procedures related to the
  !! Migdal-Eliashberg (MigEl) solver environment.
  
  use params, only: dp, k8, pi, kB
  use misc, only: subtitle, print_message, exit_with_message, write2file_rank1_real, &
       twonorm, distribute_points
  use numerics_module, only: numerics
  use electron_module, only: electron
  use wannier_module, only: epw_wannier
  use eliashberg, only: calculate_isotropic_lambda
  
  implicit none

  private
  public MigEl_sc

  type MigEl_sc
     !! Data and procedures related to the Migdal-Eliashberg equations solver environment.
     
     integer(k8) :: numqp
     !! Number of point on quasiparticle energy grid
     real(dp), allocatable :: qp_ens(:)
     !! Uniform quasiparticle energy mesh
     integer(k8) :: qp_cutoff
     !! Quasiparticle energy cutoff (factor that multiplies the highest phonon energy)
     integer(k8) :: nummatsubara !*
     !! Number of points on Matsubara mesh
     integer(k8) :: nummatsubara_upper !*
     !! Number of points on upper plane Matsubara mesh
     integer(k8) :: matsubara_cutoff
     !! Matsubara energy cutoff (factor of highest phonon energy)
     real(dp), allocatable :: bose_matsubara_ens(:) !*
     !! Uniform Bosonic Matsubara mesh
     real(dp), allocatable :: fermi_matsubara_ens(:) !*
     !! Number of point on phonon energy grid
     real(dp), allocatable :: omegas(:)
     !! Uniform Fermionic Matsubara mesh
     integer(k8) :: numomega
     !! Uniform phonon energy mesh
     real(dp) :: omegalog
     !! Logarithmic average of phonon energy
     real(dp) :: iso_lambda0
     !! Standard, isotropic e-ph coupling
     real(dp) :: domega
     !! Uniform bosonic mesh energy difference
     real(dp) :: Tstart, Tend, dT
     !! Temperature sweep: start, end, difference
     real(dp) :: mustar
     !! Dimensionless Coulomb pseudopotential parameter
     real(dp) :: MAD_Tc !*
     !! Superconducting transition temperature in the McMillan-Allen-Dynes (MAD) theory
     real(dp) :: BCS_delta !*
     !! Superconducting gap from the BCS theory using the MAD Tc
     logical :: isotropic
     !! Use isotropic approximation?

   contains
     procedure, public :: initialize, calculate_MAD_theory, calculate_iso_MigEl_theory
     procedure, private :: generate_real_ens_meshes, generate_matsubara_meshes
  end type MigEl_sc

contains

  subroutine initialize(self, max_ph_en)
    !! Read input file and setup the T-independent part of the MigEl environment.
    
    class(MigEl_sc), intent(out) :: self
    real(dp), intent(in) :: max_ph_en

    !Local variables
    integer(k8) :: qp_cutoff, matsubara_cutoff
    real(dp) :: domega, Tstart, Tend, &
         dT, mustar
    logical(dp) :: isotropic
    
    namelist /superconductivity/ domega, matsubara_cutoff, qp_cutoff, &
         Tstart, Tend, dT, mustar, isotropic

    call subtitle("Setting up Migdal-Eliashberg solver environment...")

    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Read superconductivity-related information
    qp_cutoff = 0_k8
    matsubara_cutoff = 0_k8
    domega = 0.0_dp
    Tstart = 0.0_dp
    Tend = 0.0_dp
    dT = 0.0_dp
    mustar = 0.0_dp
    isotropic = .false.
    read(1, nml = superconductivity)
    if(domega*qp_cutoff*matsubara_cutoff*Tstart* &
         Tend*dT == 0) then
       call exit_with_message('Bad input(s) provided for superconductivity solver.')
    end if
    self%qp_cutoff = qp_cutoff
    self%matsubara_cutoff = matsubara_cutoff
    self%domega = domega
    self%Tstart = Tstart
    self%Tend = Tend
    self%dT = dT
    self%mustar = mustar
    self%isotropic = isotropic
    
    !Set up meshes
    call generate_real_ens_meshes(self, max_ph_en)
    
    !Print out information.
    if(this_image() == 1) then
       write(*, "(A, 1E16.8, A)") "Matsubara energy cut-off = ", self%matsubara_cutoff*max_ph_en, " eV"
       write(*, "(A, 1E16.8, A)") "Quasiparticle energy cut-off = ", self%qp_cutoff*max_ph_en, " eV"
       write(*, "(A, 1E16.8, A)") "Bosonic energy mesh spacing = ", self%domega, " eV"
       write(*, "(A, 1E16.8, A)") "Migdal-Eliashberg solver first temperature = ", self%Tstart, " K"
       write(*, "(A, 1E16.8, A)") "Migdal-Eliashberg solver last temperature = ", self%Tend, " K"
       write(*, "(A, 1E16.8, A)") "Migdal-Eliashberg solver temperature difference = ", self%dT, " K"
       write(*, "(A, 1E16.8)") "Coulomb pseudopotential parameter = ", self%mustar
       write(*, "(A, L)") "Use isotropic Migdal-Eliashberg theory: ", self%isotropic
    end if

    sync all
  end subroutine initialize

  subroutine calculate_MAD_theory(self)
    !! Calculate the supercondting gap and transition temperature
    !! using the McMillan-Allen-Dynes (MAD) theory.
    !! P. B. Allen and R. C. Dynes Phys. Rev. B 12, 905 (1975).

    class(MigEl_sc), intent(inout) :: self
    
    !McMillan-Allen-Dynes Tc
    self%MAD_Tc = self%omegalog/1.2_dp*exp( -1.04_dp*(1.0_dp + self%iso_lambda0)/ &
         (self%iso_lambda0 - self%mustar*(1.0_dp + 0.62_dp*self%iso_lambda0)))/kB

    !BCS gap in the weak-coupling limit from the MAD Tc
    self%BCS_Delta = 1.72_dp*kB*self%MAD_Tc

    if(this_image() == 1) then
       write(*,"(A, (1E16.8, x), A)") ' McMillan-Allen-Dynes Tc =', &
            self%MAD_Tc, ' K'
       write(*,"(A, (1E16.8, x), A)") ' BCS gap =', &
            self%BCS_delta*1.0e3_dp, ' meV'
    end if
  end subroutine calculate_MAD_theory

  subroutine calculate_iso_MigEl_theory(self, el, wann, num, max_ph_en)
    !! Solver for the isotropic Migdal-Eliashberg equations.

    class(MigEl_sc), intent(inout) :: self
    type(electron), intent(in) :: el
    type(epw_wannier), intent(in) :: wann
    type(numerics), intent(in) :: num
    real(dp), intent(in) :: max_ph_en

    !Local variables
    logical :: in_T_bracket
    real(dp) :: T, norm_Z, norm_Zold, norm_Delta, norm_Deltaold
    real(dp), allocatable :: &
         iso_quasi_Delta(:), iso_quasi_Z(:), iso_quasi_dos(:), &
         iso_matsubara_lambda(:), iso_matsubara_Delta(:), iso_matsubara_Z(:)
    integer(k8) :: iter, nstates_irred
    
    !Total number of IBZ blocks states
    nstates_irred = el%nwv_irred*wann%numwannbands
    
    if(self%isotropic) then
       call print_message("Solving the isotropic Migdal-Eliashberg equations...")
       
       allocate(iso_quasi_Delta(self%numqp), &
            iso_quasi_Z(self%numqp), iso_quasi_dos(self%numqp))
    else
       !Place holder.
       !call print_message("Solving the anisotropic Migdal-Eliashberg equations...")
    end if

    !Calculate for all temperatures in the provided bracket
    T = self%Tstart
    in_T_bracket = .true.
    do while(in_T_bracket)
       if(this_image() == 1) then
          write(*, "(A, 1E16.8, A)") " Temperature = ", T, " K"
       end if
       
       !Generate the Matsubara meshes
       call self%generate_matsubara_meshes(T, max_ph_en)

       if(self%isotropic) then
          if(allocated(iso_matsubara_lambda)) deallocate(iso_matsubara_lambda)
          if(allocated(iso_matsubara_Delta)) deallocate(iso_matsubara_Delta)
          if(allocated(iso_matsubara_Z)) deallocate(iso_matsubara_Z)

          allocate(iso_matsubara_lambda(self%nummatsubara))
          allocate(iso_matsubara_Delta(self%nummatsubara))
          allocate(iso_matsubara_Z(self%nummatsubara))

          !Approximate initial gap to start off the iterator
          iso_matsubara_Delta(:) = self%BCS_delta
          norm_Delta = twonorm(iso_matsubara_Delta)

          iso_matsubara_Z = 1.0_dp !For the sake of non-zero initial norm
          norm_Z = twonorm(iso_matsubara_Z)

          !Calculate isotropic lambda function
          call calculate_isotropic_lambda(wann, num, self%omegas, self%bose_matsubara_ens, iso_matsubara_lambda)
       else
          !Place holder.
       end if

       if(this_image() == 1) write(*, "(A)") "   Iterating..."
       do iter = 1, num%maxiter
          if(self%isotropic) then
             call iterate_iso_matsubara_Z(iso_matsubara_lambda, self%fermi_matsubara_ens, &
                  iso_matsubara_Delta, iso_matsubara_Z, T)
             call iterate_iso_matsubara_Delta(iso_matsubara_lambda, self%fermi_matsubara_ens, &
                  iso_matsubara_Delta, iso_matsubara_Z, T, self%mustar)
          else
             !Place holder.
          end if

          norm_Zold = norm_Z
          norm_Deltaold = norm_Delta

          if(self%isotropic) then
             norm_Z = twonorm(iso_matsubara_Z)
             norm_Delta = twonorm(iso_matsubara_Delta)
          else
             !norm_Z = norm_tensor(Z)
             !norm_Delta = norm_tensor(Delta)
          end if

          !Output norms every 100 iterations
          if(this_image() == 1 .and. mod(iter, 100) == 0) then
             write(*, "(A, (1E16.8, x, 1E16.8))") "   ||Z||, ||Delta|| = ", norm_Z, norm_Delta
          end if

          !Iteration breaker
          if((abs(norm_Z - norm_Zold)/norm_Zold < num%conv_thres) .and. &
               (abs(norm_Delta - norm_Deltaold)/norm_Deltaold < num%conv_thres)) then
             if(this_image() == 1) write(*, "(A, (I5))") "   Convergence reached after iter", iter
             if(this_image() == 1) write(*, "(A, (1E16.8))") "   ||Delta|| = ", norm_Delta
             exit
          end if
       end do !iteration number
       
       !write(*, "(A)") "  Performing analytical continuation..."
       !TODO Perform analytic continuation
       
       !TODO Check phase
       
       !Next temperature
       T = T + self%dT

       !Check temperature bracket
       if(T > self%Tend .or. T < self%Tstart) in_T_bracket = .false.
    end do
  end subroutine calculate_iso_MigEl_theory

  subroutine iterate_iso_matsubara_Z(iso_matsubara_lambda, fermi_matsubara_ens, &
       Delta, Z, T)
    !! Iterate the isotropic Matsubara mass renormalization function for
    !! a given gap (Delta) at a given temperature (T) in K.

    real(dp), intent(in) :: iso_matsubara_lambda(:), Delta(:), &
         fermi_matsubara_ens(:), T
    real(dp), intent(out) :: Z(:)

    !Internal variables
    integer(k8) :: j, jp, nummatsubara, start, end, chunk, num_active_images
    real(dp) :: aux, pikBT

    nummatsubara = size(Delta)
    
    pikBT = pi*T*kB !eV

    call distribute_points(nummatsubara, chunk, start, end, num_active_images)
    
    Z = 0.0_dp
    !Run over Matsubara energies
    do j = start, end
       !Sum over Matsubara energies
       aux = 0.0_dp
       do jp = 1, nummatsubara
          aux = aux + iso_matsubara_lambda(abs(j - jp) + 1)* &
               fermi_matsubara_ens(jp)/ &
               sqrt(fermi_matsubara_ens(jp)**2 + Delta(jp)**2)
       end do
       Z(j) = Z(j) + aux/fermi_matsubara_ens(j) 
    end do

    !Reduce Z
    call co_sum(Z)

    Z = 1.0_dp + pikBT*Z
    sync all
  end subroutine iterate_iso_matsubara_Z

  subroutine iterate_iso_matsubara_Delta(iso_matsubara_lambda, fermi_matsubara_ens, &
       Delta, Z, T, mustar)
    !! Iterates the isotropic Matsubara gap function, Delta.
    !! This will take Delta^(tau) -> Delta^(tau + 1) for a given mass enhancement, Z^(tau),
    !! tau being the "time" in the iterator sense.

    real(dp), intent(in) :: iso_matsubara_lambda(:), fermi_matsubara_ens(:), &
         Z(:), T, mustar
    real(dp), intent(inout) :: Delta(:)

    !Internal variables
    integer(k8) :: j, jp, nummatsubara, start, end, chunk, num_active_images
    real(dp) :: aux, lambda_aux, pikBT
    real(dp), allocatable :: old_Delta(:)

    nummatsubara = size(Delta)

    allocate(old_Delta(nummatsubara))

    pikBT = pi*T*kB !eV

    call distribute_points(nummatsubara, chunk, start, end, num_active_images)
    
    old_Delta = Delta
    Delta = 0.0_dp
    !Run over Matsubara energies
    do j = start, end
       !Sum over Matsubara energies
       aux = 0.0_dp
       do jp = 1, nummatsubara
          lambda_aux = iso_matsubara_lambda(abs(j - jp) + 1)

          aux = aux + (lambda_aux - mustar)*old_Delta(jp)/ &
               sqrt(fermi_matsubara_ens(jp)**2 + old_Delta(jp)**2)
       end do
       Delta(j) = Delta(j) + aux
    end do

    !Reduce Delta
    call co_sum(Delta)

    Delta = Delta*pikBT/Z
    sync all
  end subroutine iterate_iso_matsubara_Delta
  
  subroutine generate_real_ens_meshes(self, max_ph_en)
    !! Create uniform mesh of phonon energies and quasiparticle energies

    class(MigEl_sc), intent(inout) :: self
    real(dp), intent(in) :: max_ph_en

    !Local variables
    integer(k8) :: i

    call print_message("Creating uniform phonon energy mesh...")

    !Number of phonon energy points in mesh
    self%numomega = ceiling((max_ph_en + 5.0e-3_dp)/self%domega)

    !Create uniform phonon energy mesh
    if(allocated(self%omegas)) deallocate(self%omegas)
    allocate(self%omegas(self%numomega))
    self%omegas(1) = 1.0e-5_dp !avoid zero phonon energy
    do i = 2, self%numomega
       self%omegas(i) = self%omegas(i - 1) + self%domega
    end do
    
    !Number of quasiparticle energy points in mesh
    self%numqp = self%qp_cutoff*ceiling(max_ph_en/self%domega)

    !Create uniform quasiparticle energy mesh
    allocate(self%qp_ens(self%numqp))
    self%qp_ens(1) = 0
    do i = 2, self%numqp
       self%qp_ens(i) = self%qp_ens(i - 1) + self%domega
    end do

    !Write energy meshes to file
    call write2file_rank1_real('omegas', self%omegas)
    call write2file_rank1_real('quasiparticle_ens', self%qp_ens)
  end subroutine generate_real_ens_meshes

  subroutine generate_matsubara_meshes(self, temp, max_ph_en)
    !! Calculate the Fermionic and Bosonic Matsubara frequency meshes
    !! for a given temperature in Kelvins. 

    class(MigEl_sc), intent(inout) :: self
    real(dp), intent(in) :: temp, max_ph_en
    
    !Local variables
    integer(k8) :: l, halfloc
    real(dp) :: dmatsubara, invbeta

    call print_message("   Creating Matsubara meshes...")

    !Temperature energy
    invbeta = kB*temp !eV
    
    !Matsubara energy mesh spacing
    dmatsubara = pi*invbeta !eV

    !The length of the Matsubara mesh is T-dependent.
    !Set number of Matsubara energy points.
    ! The ratio of the largest energy to the largest phonon energy = matsubara_cutoff
    self%nummatsubara = self%matsubara_cutoff*ceiling(max_ph_en/dmatsubara)
    if(mod(self%nummatsubara, 2) == 1) self%nummatsubara = self%nummatsubara + 1 !enforce evenness
    
    !Bose iomega_l = i2l.pi/beta, upper plane: l = 0, 1, 2, ...
    if(allocated(self%bose_matsubara_ens)) deallocate(self%bose_matsubara_ens)
    allocate(self%bose_matsubara_ens(self%nummatsubara))
    do l = 1, self%nummatsubara 
       self%bose_matsubara_ens(l) = 2*(l - 1)*dmatsubara 
    end do

    !Fermi iomega_j = i(2j + 1)pi/beta, j = ...-3, -2, -1, 0, 1, 2, ...
    if(allocated(self%fermi_matsubara_ens)) deallocate(self%fermi_matsubara_ens)
    allocate(self%fermi_matsubara_ens(self%nummatsubara))
    halfloc = self%nummatsubara/2+1
    do l = 1, self%nummatsubara
       self%fermi_matsubara_ens(l) = (2*(l-halfloc) + 1)*dmatsubara
    end do

    !Set number of fermionic Matsubara points on upper half plane
    !j = 0, 1, 2, ...
    self%nummatsubara_upper = self%nummatsubara/2 + 1

    !Write energy meshes to file
    !call write2file_rank1_real('bose_matsubara_ens', self%bose_matsubara_ens)
    !call write2file_rank1_real('fermi_matsubara_ens', self%fermi_matsubara_ens)
  end subroutine generate_matsubara_meshes
end module MigEl_sc_module
