! Copyright 2022 elphbolt contributors.
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
  
  use params, only: r64, i64, pi, kB, oneI
  use misc, only: subtitle, print_message, exit_with_message, write2file_rank1_real, &
       twonorm, distribute_points, Pade_continued, demux_state, mux_state
  use numerics_module, only: numerics
  use electron_module, only: electron
  use wannier_module, only: wannier
  use eliashberg, only: calculate_iso_Matsubara_lambda, calculate_aniso_Matsubara_lambda
  
  implicit none

  !external chdir

  private
  public MigEl_sc

  type MigEl_sc
     !! Data and procedures related to the Migdal-Eliashberg equations solver environment.
     
     integer(i64) :: numqp
     !! Number of point on quasiparticle energy grid
     real(r64), allocatable :: qp_ens(:)
     !! Uniform quasiparticle energy mesh
     integer(i64) :: qp_cutoff
     !! Quasiparticle energy cutoff (factor that multiplies the highest phonon energy)
     integer(i64) :: nummatsubara
     !! Number of points on Matsubara mesh
     integer(i64) :: nummatsubara_upper
     !! Number of points on upper plane Matsubara mesh
     integer(i64) :: matsubara_cutoff
     !! Matsubara energy cutoff (factor of highest phonon energy)
     real(r64), allocatable :: bose_matsubara_ens(:)
     !! Uniform Bosonic Matsubara mesh
     real(r64), allocatable :: fermi_matsubara_ens(:)
     !! Number of point on phonon energy grid
     real(r64), allocatable :: omegas(:)
     !! Uniform Fermionic Matsubara mesh
     integer(i64) :: numomega
     !! Uniform phonon energy mesh
     real(r64) :: omegalog
     !! Logarithmic average of phonon energy
     real(r64) :: iso_lambda0
     !! Standard, isotropic e-ph coupling
     real(r64) :: domega
     !! Uniform bosonic mesh energy difference
     real(r64) :: Tstart, Tend, dT
     !! Temperature sweep: start, end, difference
     real(r64) :: mustar
     !! Dimensionless Coulomb pseudopotential parameter
     real(r64) :: MAD_Tc
     !! Superconducting transition temperature in the McMillan-Allen-Dynes (MAD) theory
     real(r64) :: BCS_delta
     !! Superconducting gap from the BCS theory using the MAD Tc
     logical :: isotropic
     !! Use isotropic approximation?
     logical :: use_external_eps
     !! Use user generated |epsilon|^2 to screen a2F?

   contains
     procedure, public :: initialize, calculate_MAD_theory, calculate_MigEl_theory
     procedure, private :: generate_real_ens_meshes, generate_matsubara_meshes
  end type MigEl_sc

contains

  subroutine initialize(self, max_ph_en)
    !! Read input file and setup the T-independent part of the MigEl environment.
    
    class(MigEl_sc), intent(out) :: self
    real(r64), intent(in) :: max_ph_en

    !Local variables
    real(r64) :: domega, Tstart, Tend, &
         dT, mustar, qp_cutoff, matsubara_cutoff
    logical(r64) :: isotropic, use_external_eps
    
    namelist /superconductivity/ domega, matsubara_cutoff, qp_cutoff, &
         Tstart, Tend, dT, mustar, isotropic, use_external_eps

    call subtitle("Setting up Migdal-Eliashberg solver environment...")

    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Read superconductivity-related information
    qp_cutoff = 0_r64
    matsubara_cutoff = 0_r64
    domega = 0.0_r64
    Tstart = 0.0_r64
    Tend = 0.0_r64
    dT = 0.0_r64
    mustar = 0.0_r64
    isotropic = .false.
    use_external_eps = .false.
    read(1, nml = superconductivity)
    if(domega*qp_cutoff*matsubara_cutoff*Tstart* &
         Tend*dT == 0) then
       call exit_with_message('Bad input(s) provided for superconductivity solver.')
    end if
    self%qp_cutoff = nint(qp_cutoff)
    self%matsubara_cutoff = nint(matsubara_cutoff)
    self%domega = domega
    self%Tstart = Tstart
    self%Tend = Tend
    self%dT = dT
    self%mustar = mustar
    self%isotropic = isotropic
    self%use_external_eps = use_external_eps

    if(.not. self%isotropic .and. self%use_external_eps) &
         call exit_with_message('External screening for the anisotropic case is not supported. Exiting.')
    
    !Set up meshes
    call generate_real_ens_meshes(self, max_ph_en)
    
    !Print out information.
    if(this_image() == 1) then
       write(*, "(A, 1E16.8, A)") "Matsubara energy cut-off = ", self%matsubara_cutoff*max_ph_en, " eV"
       write(*, "(A, 1E16.8, A)") "Quasiparticle energy cut-off = ", self%qp_cutoff*max_ph_en, " eV"
       write(*, "(A, 1E16.8, A)") "Bosonic energy mesh spacing = ", self%domega, " eV"
       write(*, "(A, 1E16.8, A)") "First temperature = ", self%Tstart, " K"
       write(*, "(A, 1E16.8, A)") "Last temperature = ", self%Tend, " K"
       write(*, "(A, 1E16.8, A)") "Temperature step = ", self%dT, " K"
       write(*, "(A, 1E16.8)") "Coulomb pseudopotential parameter = ", self%mustar
       write(*, "(A, L)") "Use isotropic Migdal-Eliashberg theory: ", self%isotropic
       write(*, "(A, L)") "Use user generated |epsilon|^2 to screen a2F: ", self%use_external_eps
    end if

    sync all
  end subroutine initialize

  subroutine calculate_MAD_theory(self)
    !! Calculate the supercondting gap and transition temperature
    !! using the McMillan-Allen-Dynes (MAD) theory.
    !! P. B. Allen and R. C. Dynes Phys. Rev. B 12, 905 (1975).

    class(MigEl_sc), intent(inout) :: self
    
    !McMillan-Allen-Dynes Tc
    self%MAD_Tc = self%omegalog/1.2_r64*exp( -1.04_r64*(1.0_r64 + self%iso_lambda0)/ &
         (self%iso_lambda0 - self%mustar*(1.0_r64 + 0.62_r64*self%iso_lambda0)))/kB

    !BCS gap in the weak-coupling limit from the MAD Tc
    self%BCS_Delta = 1.72_r64*kB*self%MAD_Tc

    if(this_image() == 1) then
       write(*,"(A, (1E16.8, x), A)") 'McMillan-Allen-Dynes Tc =', &
            self%MAD_Tc, ' K'
       write(*,"(A, (1E16.8, x), A)") 'BCS gap =', &
            self%BCS_delta*1.0e3_r64, ' meV'
    end if
  end subroutine calculate_MAD_theory

  subroutine calculate_MigEl_theory(self, el, wann, num, max_ph_en)
    !! Solve the Migdal-Eliashberg equations.

    class(MigEl_sc), intent(inout) :: self
    type(electron), intent(in) :: el
    type(wannier), intent(in) :: wann
    type(numerics), intent(in) :: num
    real(r64), intent(in) :: max_ph_en

    !Local variables
    logical :: in_T_bracket
    real(r64) :: T, norm_Z, norm_Zold, norm_Delta, norm_Deltaold, aux
    real(r64), allocatable :: quasi_dos(:), &
         iso_matsubara_lambda(:), iso_matsubara_Delta(:), iso_matsubara_Z(:), &
         aniso_matsubara_Delta(:, :), aniso_matsubara_Z(:, :)
    complex(r64), allocatable :: iso_quasi_Delta(:), iso_quasi_Z(:), &
         aniso_quasi_Delta(:, :), aniso_quasi_Z(:, :)
    integer(i64) :: iter, nstates_irred, i, istate, m, ik
    character(len = 1024) :: filename, numcols
    
    !Total number of IBZ blocks states
    nstates_irred = el%nwv_irred*wann%numwannbands

    allocate(quasi_dos(self%numqp))
    
    if(self%isotropic) then
       call print_message("Solving the isotropic Migdal-Eliashberg equations...")
       
       allocate(iso_quasi_Delta(self%numqp), iso_quasi_Z(self%numqp))
    else
       call print_message("Solving the anisotropic Migdal-Eliashberg equations...")

       allocate(aniso_quasi_Delta(nstates_irred, self%numqp), &
            aniso_quasi_Z(nstates_irred, self%numqp))
    end if

    !Calculate for all temperatures in the provided bracket
    T = self%Tstart
    in_T_bracket = .true.
    do while(in_T_bracket)
       if(this_image() == 1) then
          write(*, "(A)") "---->"
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

          iso_matsubara_Z = 1.0_r64 !For the sake of non-zero initial norm
          norm_Z = twonorm(iso_matsubara_Z)

          !Calculate isotropic lambda function
          call calculate_iso_matsubara_lambda(wann, num, self%omegas, self%bose_matsubara_ens, iso_matsubara_lambda)
       else
          if(allocated(aniso_matsubara_Delta)) deallocate(aniso_matsubara_Delta)
          if(allocated(aniso_matsubara_Z)) deallocate(aniso_matsubara_Z)

          allocate(aniso_matsubara_Delta(nstates_irred, self%nummatsubara))
          allocate(aniso_matsubara_Z(nstates_irred, self%nummatsubara))

          !Approximate initial gap to start off the iterator
          aniso_matsubara_Delta(:, :) = self%BCS_delta
          norm_Delta = twonorm(aniso_matsubara_Delta)

          aniso_matsubara_Z = 1.0_r64 !For the sake of non-zero initial norm
          norm_Z = twonorm(aniso_matsubara_Z)

          !Calculate anisotropic lambda function
          call calculate_aniso_matsubara_lambda(wann, num, el, self%omegas, self%bose_matsubara_ens)
       end if

       if(this_image() == 1) write(*, "(A)") "   Iterating..."
       do iter = 1, max(500, num%maxiter) !At least 500 iterations.
          if(self%isotropic) then
             call iterate_iso_matsubara_Z(iso_matsubara_lambda, self%fermi_matsubara_ens, &
                  iso_matsubara_Delta, iso_matsubara_Z, T)
             call iterate_iso_matsubara_Delta(iso_matsubara_lambda, self%fermi_matsubara_ens, &
                  iso_matsubara_Delta, iso_matsubara_Z, T, self%mustar)
          else
             call iterate_aniso_matsubara_Z(el, wann, num, self%fermi_matsubara_ens, &
                  aniso_matsubara_Delta, aniso_matsubara_Z, T)
             call iterate_aniso_matsubara_Delta(el, wann, num, self%fermi_matsubara_ens, &
                  aniso_matsubara_Delta, aniso_matsubara_Z, T, self%mustar)
          end if

          norm_Zold = norm_Z
          norm_Deltaold = norm_Delta

          if(self%isotropic) then
             norm_Z = twonorm(iso_matsubara_Z)
             norm_Delta = twonorm(iso_matsubara_Delta)
          else
             norm_Z = twonorm(aniso_matsubara_Z)
             norm_Delta = twonorm(aniso_matsubara_Delta)
          end if

          !Zero out norm_Delata if it is smaller 1 micro eV 
          if(norm_Delta < 1.0e-6_r64) norm_Delta = 0.0_r64

          !Output norms every 100 iterations
          if(this_image() == 1 .and. mod(iter, 100) == 0) then
             write(*, "(A, (1E16.8, x, 1E16.8))") "   ||Z||, ||Delta|| = ", norm_Z, norm_Delta
          end if

          !Iteration breaker
          if((abs(norm_Z - norm_Zold)/norm_Zold < num%conv_thres) .and. &
               ((abs(norm_Delta - norm_Deltaold)/norm_Deltaold < num%conv_thres) &
               .or. norm_Delta == 0.0_r64)) then
             if(this_image() == 1) write(*, "(A, (I5))") "   Convergence reached after iter", iter
             if(this_image() == 1) write(*, "(A, (1E16.8))") "   ||Delta|| = ", norm_Delta
             exit
          end if
       end do !iteration number
       
       !Perform analytic continuation to positive real energies
       if(this_image() == 1) &
            write(*, "(A)") "  Performing analytic continuation..."

       if(self%isotropic) then
          iso_quasi_Delta = Pade_continued(&
               oneI*self%fermi_matsubara_ens(self%nummatsubara_upper:self%nummatsubara), &
               iso_matsubara_Delta(self%nummatsubara_upper:self%nummatsubara), self%qp_ens)

          !Reduced quasiparticle density of states
          !Eq. 11 of H.J. Choi et al. Physica C 385 (2003) 66–74
          quasi_dos = real((self%qp_ens + oneI*1.0e-6_r64)/ &
               sqrt((self%qp_ens + oneI*1.0e-6_r64)**2 - iso_quasi_Delta**2))

          !Write quasiparticle Delta and reduced DOS as text data to file
          if(this_image() == 1) then
             call chdir(num%cwd)
             write (filename, '(f10.3)') T
             filename = 'iso_quasiparticle_Delta.T' // trim(adjustl(filename))
             write(numcols, "(I0)") 3
             open(1,file = trim(filename), status = 'replace')
             do i = 1, self%numqp
                write(1, "("//trim(adjustl(numcols))//"E20.10)") self%qp_ens(i), &
                     iso_quasi_Delta(i)
             end do
             close(1)

             write (filename, '(f10.3)') T
             filename = 'iso_quasiparticle_DOS.T' // trim(adjustl(filename))
             write(numcols, "(I0)") 2
             open(1,file = trim(filename), status = 'replace')
             do i = 1, self%numqp
                write(1, "("//trim(adjustl(numcols))//"E20.10)") self%qp_ens(i), &
                     quasi_dos(i)
             end do
             close(1)
          end if
          sync all
       else
          do istate = 1, nstates_irred
             !Demux state index into band (m) and wave vector (ik) indices
             call demux_state(istate, wann%numwannbands, m, ik)

             !Apply energy window to IBZ blocks electron
             if(abs(el%ens_irred(ik, m) - el%enref) > el%fsthick) cycle

             aniso_quasi_Delta(istate, :) = Pade_continued(&
                  oneI*self%fermi_matsubara_ens(self%nummatsubara_upper:self%nummatsubara), &
                  aniso_matsubara_Delta(istate, self%nummatsubara_upper:self%nummatsubara), self%qp_ens)
          end do

          !TODO Need to decide if this should be printed out since it can be a pretty large file
          
          !Reduced quasiparticle density of states
          !Eq. 11 of H.J. Choi et al. Physica C 385 (2003) 66–74
          do i = 1, self%numqp
             aux = 0.0_r64
             do istate = 1, nstates_irred
                !Demux state index into band (m) and wave vector (ik) indices
                call demux_state(istate, wann%numwannbands, m, ik)
                
                aux = aux + el%nequiv(ik)*el%Ws_irred(ik, m)*&
                     real((self%qp_ens(i) + oneI*1.0e-6_r64)/ &
                     sqrt((self%qp_ens(i) + oneI*1.0e-6_r64)**2 - aniso_quasi_Delta(istate, i)**2))
             end do
             quasi_dos(i) = aux
          end do

          !Write quasiparticle reduced DOS as text data to file
          if(this_image() == 1) then
             call chdir(num%cwd)
             write (filename, '(f10.3)') T
             filename = 'aniso_quasiparticle_DOS.T' // trim(adjustl(filename))
             write(numcols, "(I0)") 2
             open(1,file = trim(filename), status = 'replace')
             do i = 1, self%numqp
                write(1, "("//trim(adjustl(numcols))//"E20.10)") self%qp_ens(i), &
                     quasi_dos(i)
             end do
             close(1)
          end if
          sync all
       end if
              
       !Next temperature
       T = T + self%dT

       !Check temperature bracket
       if(self%dT > 0) then
          if(T > self%Tend .or. T < self%Tstart) in_T_bracket = .false.
       else
          if(T < self%Tend .or. T > self%Tstart) in_T_bracket = .false.
       end if
    end do
  end subroutine calculate_MigEl_theory

  subroutine iterate_aniso_matsubara_Z(el, wann, num, fermi_matsubara_ens, &
       Delta, Z, T)
    !! Iterate the anisotropic Matsubara mass renormalization function for
    !! a given gap (Delta) at a given temperature (T) in K.

    type(electron), intent(in) :: el
    type(wannier), intent(in) :: wann
    type(numerics), intent(in) :: num
    real(r64), intent(in) :: Delta(:, :), fermi_matsubara_ens(:), T
    real(r64), intent(out) :: Z(:, :)

    !Internal variables
    integer(i64) :: j, jp, nstates_irred, nummatsubara, nprocs, start, end, chunk, &
         num_active_images, istate, count, m, n, s, ik, ikp, ikp_ibz, istatep_ibz
    real(r64) :: aux, pikBT
    real(r64), allocatable :: lambda(:), lambda_istate(:, :)
    character(len = 1024) :: filename

    nstates_irred = size(Delta(:, 1))
    nummatsubara = size(Delta(1, :))

    pikBT = pi*T*kB !eV

    allocate(lambda(nummatsubara))

    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    Z = 0.0_r64
    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over IBZ blocks states
       do istate = start, end
          !Initialize eligible process counter for this state
          count = 0

          !Demux state index into band (m) and wave vector (ik) indices
          call demux_state(istate, wann%numwannbands, m, ik)

          !Apply energy window to initial (IBZ blocks) electron
          if(abs(el%ens_irred(ik, m) - el%enref) > el%fsthick) cycle

          !Read anisotropic lambda_istate from file
          call chdir(trim(adjustl(num%scdir)))
          write (filename, '(I9)') istate
          filename = 'lambda.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'old', access = 'stream')
          read(1) nprocs
          if(allocated(lambda_istate)) deallocate(lambda_istate)
          allocate(lambda_istate(nprocs, nummatsubara))
          if(nprocs > 0) read(1) lambda_istate
          close(1)

          !Sum over final (FBZ blocks) electron wave vectors
          do ikp = 1, el%nwv
             !Grab IBZ vector index, ikp_ibz, corresponding to this FBZ vector
             ikp_ibz = el%fbz2ibz_map(ikp)

             !Sum over final electron bands
             do n = 1, wann%numwannbands
                !Apply energy window to final electron
                if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle

                !Mux primed state index
                istatep_ibz = mux_state(wann%numwannbands, n, ikp_ibz)

                !Branch summed lambda
                lambda = 0.0_r64
                do s = 1, wann%numbranches
                   !Increment anisotropic lambda_istate processes counter
                   count = count + 1
                   lambda(:) = lambda(:) + lambda_istate(count, :)
                end do

                !Run over Matsubara axis
                do j = 1, nummatsubara
                   !Sum over Matsubara axis
                   aux = 0.0_r64
                   do jp = 1, nummatsubara
                      aux = aux + el%Ws(ikp, n)* &
                           lambda(abs(j - jp) + 1)*fermi_matsubara_ens(jp)/ &
                           sqrt(fermi_matsubara_ens(jp)**2 + Delta(istatep_ibz, jp)**2)
                   end do
                   Z(istate, j) = Z(istate, j) + aux/fermi_matsubara_ens(j) 
                end do
             end do
          end do
       end do
    end if

    sync all
    call co_sum(Z)
    sync all
    
    Z = 1.0_r64 + Z*pikBT
    sync all
  end subroutine iterate_aniso_matsubara_Z

  subroutine iterate_aniso_matsubara_Delta(el, wann, num, fermi_matsubara_ens, &
       Delta, Z, T, mustar)
    !! Iterate the isotropic Matsubara gap function, Delta.
    !! This will take Delta^(tau) -> Delta^(tau + 1) for a given mass enhancement, Z^(tau),
    !! tau being the "time" in the iterator sense.

    type(electron), intent(in) :: el
    type(wannier), intent(in) :: wann
    type(numerics), intent(in) :: num
    real(r64), intent(in) :: Z(:, :), fermi_matsubara_ens(:), T, mustar
    real(r64), intent(inout) :: Delta(:, :)

    !Internal variables
    integer(i64) :: j, jp, nstates_irred, nummatsubara, nprocs, start, end, chunk, &
         num_active_images, istate, count, m, n, s, ik, ikp, ikp_ibz, istatep_ibz
    real(r64) :: aux, pikBT
    real(r64), allocatable :: lambda(:), lambda_istate(:, :), old_Delta(:, :)
    character(len = 1024) :: filename

    nstates_irred = size(Delta(:, 1))
    nummatsubara = size(Delta(1, :))

    pikBT = pi*T*kB !eV

    allocate(lambda(nummatsubara))

    allocate(old_Delta(nstates_irred, nummatsubara))
    
    call distribute_points(nstates_irred, chunk, start, end, num_active_images)

    old_Delta = Delta
    Delta = 0.0_r64

    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over IBZ blocks states
       do istate = start, end
          !Initialize eligible process counter for this state
          count = 0

          !Demux state index into band (m) and wave vector (ik) indices
          call demux_state(istate, wann%numwannbands, m, ik)

          !Apply energy window to initial (IBZ blocks) electron
          if(abs(el%ens_irred(ik, m) - el%enref) > el%fsthick) cycle

          !Read anisotropic lambda_istate from file
          call chdir(trim(adjustl(num%scdir)))
          write (filename, '(I9)') istate
          filename = 'lambda.istate'//trim(adjustl(filename))
          open(1, file = trim(filename), status = 'old', access = 'stream')
          read(1) nprocs
          if(allocated(lambda_istate)) deallocate(lambda_istate)
          allocate(lambda_istate(nprocs, nummatsubara))
          if(nprocs > 0) read(1) lambda_istate
          close(1)

          !Sum over final (FBZ blocks) electron wave vectors
          do ikp = 1, el%nwv
             !Grab IBZ vector index, ikp_ibz, corresponding to this FBZ vector
             ikp_ibz = el%fbz2ibz_map(ikp)

             !Sum over final electron bands
             do n = 1, wann%numwannbands
                !Apply energy window to final electron
                if(abs(el%ens(ikp, n) - el%enref) > el%fsthick) cycle

                !Mux primed state index
                istatep_ibz = mux_state(wann%numwannbands, n, ikp_ibz)

                !Branch summed lambda
                lambda = 0.0_r64
                do s = 1, wann%numbranches
                   !Increment anisotropic lambda_istate processes counter
                   count = count + 1
                   lambda(:) = lambda(:) + lambda_istate(count, :)
                end do

                !Run over Matsubara axis
                do j = 1, nummatsubara
                   !Sum over Matsubara axis
                   aux = 0.0_r64
                   do jp = 1, nummatsubara
                      aux = aux + el%Ws(ikp, n)* &
                           (lambda(abs(j - jp) + 1) - mustar)*old_Delta(istatep_ibz, jp)/ &
                           sqrt(fermi_matsubara_ens(jp)**2 + old_Delta(istatep_ibz, jp)**2)
                   end do
                   Delta(istate, j) = Delta(istate, j) + aux 
                end do
             end do
          end do
       end do
    end if

    sync all
    call co_sum(Delta)
    sync all
    
    Delta = pikBT*Delta/Z
    sync all
  end subroutine iterate_aniso_matsubara_Delta
  
  subroutine iterate_iso_matsubara_Z(iso_matsubara_lambda, fermi_matsubara_ens, &
       Delta, Z, T)
    !! Iterate the isotropic Matsubara mass renormalization function for
    !! a given gap (Delta) at a given temperature (T) in K.
    
    real(r64), intent(in) :: iso_matsubara_lambda(:), Delta(:), &
         fermi_matsubara_ens(:), T
    real(r64), intent(out) :: Z(:)

    !Internal variables
    integer(i64) :: j, jp, nummatsubara, start, end, chunk, num_active_images
    real(r64) :: aux, pikBT

    nummatsubara = size(Delta)
    
    pikBT = pi*T*kB !eV

    call distribute_points(nummatsubara, chunk, start, end, num_active_images)

    Z = 0.0_r64

    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over Matsubara energies
       do j = start, end
          !Sum over Matsubara energies
          aux = 0.0_r64
          do jp = 1, nummatsubara
             aux = aux + iso_matsubara_lambda(abs(j - jp) + 1)* &
                  fermi_matsubara_ens(jp)/ &
                  sqrt(fermi_matsubara_ens(jp)**2 + Delta(jp)**2)
          end do
          Z(j) = Z(j) + aux/fermi_matsubara_ens(j) 
       end do
    end if

    !Reduce Z
    sync all
    call co_sum(Z)
    sync all
    
    Z = 1.0_r64 + pikBT*Z
    sync all
  end subroutine iterate_iso_matsubara_Z

  subroutine iterate_iso_matsubara_Delta(iso_matsubara_lambda, fermi_matsubara_ens, &
       Delta, Z, T, mustar)
    !! Iterate the isotropic Matsubara gap function, Delta.
    !! This will take Delta^(tau) -> Delta^(tau + 1) for a given mass enhancement, Z^(tau),
    !! tau being the "time" in the iterator sense.

    real(r64), intent(in) :: iso_matsubara_lambda(:), fermi_matsubara_ens(:), &
         Z(:), T, mustar
    real(r64), intent(inout) :: Delta(:)

    !Internal variables
    integer(i64) :: j, jp, nummatsubara, start, end, chunk, num_active_images
    real(r64) :: aux, lambda_aux, pikBT
    real(r64), allocatable :: old_Delta(:)

    nummatsubara = size(Delta)

    allocate(old_Delta(nummatsubara))

    pikBT = pi*T*kB !eV

    call distribute_points(nummatsubara, chunk, start, end, num_active_images)
    
    old_Delta = Delta
    Delta = 0.0_r64

    !Only work with the active images
    if(this_image() <= num_active_images) then
       !Run over Matsubara energies
       do j = start, end
          !Sum over Matsubara energies
          aux = 0.0_r64
          do jp = 1, nummatsubara
             lambda_aux = iso_matsubara_lambda(abs(j - jp) + 1)

             aux = aux + (lambda_aux - mustar)*old_Delta(jp)/ &
                  sqrt(fermi_matsubara_ens(jp)**2 + old_Delta(jp)**2)
          end do
          Delta(j) = Delta(j) + aux
       end do
    end if

    !Reduce Delta
    sync all
    call co_sum(Delta)
    sync all
    
    Delta = Delta*pikBT/Z
    sync all
  end subroutine iterate_iso_matsubara_Delta
  
  subroutine generate_real_ens_meshes(self, max_ph_en)
    !! Create uniform mesh of phonon energies and quasiparticle energies

    class(MigEl_sc), intent(inout) :: self
    real(r64), intent(in) :: max_ph_en

    !Local variables
    integer(i64) :: i

    call print_message("Creating uniform phonon energy mesh...")

    !Number of phonon energy points in mesh
    self%numomega = ceiling((max_ph_en + 5.0e-3_r64)/self%domega)

    !Create uniform phonon energy mesh
    if(allocated(self%omegas)) deallocate(self%omegas)
    allocate(self%omegas(self%numomega))
    self%omegas(1) = 1.0e-5_r64 !avoid zero phonon energy
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
    real(r64), intent(in) :: temp, max_ph_en
    
    !Local variables
    integer(i64) :: l, halfloc
    real(r64) :: dmatsubara, invbeta

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
  end subroutine generate_matsubara_meshes
end module MigEl_sc_module
