! Copyright 2020 elphbolt contributors.
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

module optimization
  !! Module for finding correction matrices M_T enforcing Kelvin-Onsager relations by solving an optimization problem with SLSQP
  !! The optimization problem is defined in the subroutine fun, which calculates the objective function and constraints for given correction matrices M_T (reshaped as a vector x)
  !! The main subroutine is find_correction, which initializes the optimization variable, sets up the SLSQP solver, and calls the optimization routine. 
  use slsqp_module, only: slsqp_solver
  use symmetry_module, only: symmetrize_3x3_tensor
  use precision, only: r128, r64, i64
  use misc, only: eye, compute_eigenvalues, F_norm

  implicit none

  public :: find_correction

  type(slsqp_solver) :: solver

  type :: coeffs_for_optimization
     !! Data type for storing the coefficients needed for the optimization problem. Includes the transport coefficients, temperature and symmetry elements of the system.
     
     real(r64) :: T, el_sigma(3,3), el_kappa0(3,3), el_sigmaS(3,3), el_alpha(3,3), ph_kappa(3,3), &
          ph_alpha(3,3), tot_kappa(3,3), tot_alpha(3,3)
     real(r64), allocatable :: crotations(:,:,:)
  end type coeffs_for_optimization
  
  ! Global variables for storing the coefficients needed for the optimization problem. These are set in bte.f90 before calling find_correction 
  ! and used in the subroutine fun to calculate the objective function and constraints for given correction matrices M_T. 
  ! global_coeffs_opp_B is used when magnetic field is on. 
  type(coeffs_for_optimization) :: global_coeffs, global_coeffs_opp_B

contains

  subroutine fun(me, x, f, c)
    !! Subroutine for calculating the objective function and constraints for given correction matrices M_T (reshaped as a vector x).
    !! The objective function is defined as the deviation from the Kelvin-Onsager relations, which should be minimized.
    !! The constraints are defined as the positive semidefiniteness of a matrix L constructed from the transport coefficients, 
    !! which ensures that the corrected transport coefficients satisfy the second law of thermodynamics.

    class(slsqp_solver), intent(inout) :: me
    real(r64),  intent(in) :: x(:)
    real(r64), intent(out) :: f, c(:)

    real(r64), allocatable :: M_T(:, :,:), A(:, :, :)
    real(r64) :: L(6,6), eigenvalues_real(6), eigenvalues_imag(6)
    integer :: i, n

    n = size(x)  

    ! Check that n is divisible by 9
    if (mod(n, 9) /= 0) then
       print *, "Error: size(x) must be multiple of 9"
       return
    end if

    allocate(M_T(n/9, 3, 3), A(n/9 + 1, 3, 3))
    do i=1,n/9
       M_T(i, :, :) = reshape(x((i-1)*9 + 1 : i*9), [3, 3])
       !! Symmetrize M_T to ensure that the correction matrices are symmetric, which is a physical requirement for the transport coefficients.
       call symmetrize_3x3_tensor(M_T(i, :, :), global_coeffs%crotations)
    end do
    A(1, :, :) = matmul(global_coeffs%el_sigmaS, M_T(1, :, :)) ! el_sigmaS(B)@M_I^T
    A(2, :, :) = matmul(global_coeffs%el_kappa0, M_T(1, :, :)) ! kappa_el(B)@M_I^T

    if(n == 9) then
       ! Non-magnetic case: M_T has only M_I^T

       ! Calculate tot_alpha(B) =  el_alpha(B) + ph_alpha(B) to shorten code further
       global_coeffs%tot_alpha = global_coeffs%ph_alpha+global_coeffs%el_alpha

       ! Objective function: ||el_sigmaS @ M_T - alpha^T/T||_F/||alpha/T||_F
       f = F_norm(A(1, :, :)*global_coeffs%T - transpose(global_coeffs%tot_alpha))/F_norm(global_coeffs%tot_alpha)

       ! Constraint: L matrix must be positive semidefinite
       ! L = [[sigma, el_sigmaS @ M_T],
       !      [alpha, el_kappa0 @ M_T + ph_kappa]]

       ! Top-left block: sigma  
       L(1:3, 1:3) = global_coeffs%el_sigma

       ! Top-right block: el_sigmaS @ M_T * T 
       L(1:3, 4:6) = A(1, :, :)*global_coeffs%T

       ! Bottom-left block: total_alpha  
       L(4:6, 1:3) = global_coeffs%tot_alpha

       ! Bottom-right block: (el_kappa0 @ M_T + ph_kappa)*T
       L(4:6, 4:6) = (A(2, :, :) + global_coeffs%ph_kappa)*global_coeffs%T

    elseif(n == 27) then
       ! Magnetic case with M_T = [M_I^T, M_J^T, M_F^T] 

       ! Calculate tot_alpha(-B) =  el_alpha(-B) + ph_alpha(-B) and 
       ! tot_kappa(-B) = el_kappa0(-B) + ph_kappa(-B) to shorten code further
       global_coeffs_opp_B%tot_alpha = global_coeffs_opp_B%ph_alpha+global_coeffs_opp_B%el_alpha
       global_coeffs_opp_B%tot_kappa = global_coeffs_opp_B%el_kappa0 + global_coeffs_opp_B%ph_kappa

       ! Objective function: ||el_sigmaS(B)@M_I^T - alpha(-B)^T/T||_F/||alpha(-B)/T||_F + 
       !                      + ||sigma(B)@M_J^T - sigma(-B)^T||_F/||sigma(-B)|| + 
       !                      + ||kappa_ph(B)@M_F^T + kappa_el(B)@M_I^T - kappa(-B)^T||_F/||kappa(-B)||_F
       A(3, :, :) = matmul(global_coeffs%el_sigma, M_T(2, :, :)) ! sigma(B)@M_J^T
       A(4, :, :) = matmul(global_coeffs%ph_kappa, M_T(3, :, :)) ! kappa_ph(B)@M_F^T

       f = F_norm(A(1, :, :)*global_coeffs%T - transpose(global_coeffs_opp_B%tot_alpha))/&
            F_norm(global_coeffs_opp_B%tot_alpha)
       f = f + F_norm(A(3, :, :) - transpose(global_coeffs_opp_B%el_sigma))/F_norm(global_coeffs_opp_B%el_sigma)
       f = f + F_norm(A(4, :, :) + A(2, :, :) - transpose(global_coeffs_opp_B%tot_kappa))/F_norm(global_coeffs_opp_B%tot_kappa)

       ! Constraint: L matrix must be positive semidefinite
       ! L = [[sigma @ M_J^T, el_sigmaS @ M_I^T],
       !      [el_alpha @ M_J^T + ph_alpha , el_kappa0 @ M_I^T + ph_kappa@M_F^T]]

       ! Top-left block: sigma @ M_J^T 
       L(1:3, 1:3) = A(3, :, :)

       ! Top-right block: el_sigmaS @ M_I^T 
       L(1:3, 4:6) = A(1, :, :)*global_coeffs%T

       ! Bottom-left block: el_alpha @ M_J^T + ph_alpha 
       L(4:6, 1:3) = matmul(global_coeffs%el_alpha, M_T(2, :, :)) + global_coeffs%ph_alpha

       ! Bottom-right block: el_kappa0 @ M_I^T + ph_kappa @ M_F^T 
       L(4:6, 4:6) = (A(2, :, :) + A(4, :, :))*global_coeffs%T
    else 
       print *, "Error: Unexpected size of optimization variable x: ", n
       return
    end if

    ! Compute eigenvalues of the symmetric part of L 
    call compute_eigenvalues((L + transpose(L))/2, eigenvalues_real, eigenvalues_imag)

    ! Constraint: min(real(eigenvalue)) >= 0
    c(1) = minval(eigenvalues_real)
    ! Constraint: max(|imag(eigenvalue)|) = 0
    c(2) = 1e-6 - maxval(abs(eigenvalues_imag)) 

    deallocate(M_T, A)

  end subroutine fun

  subroutine dummy_grad(me, x, g, a)
    !! This subroutine is never called because gradient_mode=1
    !! Just initialize outputs to avoid uninitialized variable warnings
    class(slsqp_solver), intent(inout) :: me
    real(r64), dimension(:), intent(in) :: x
    real(r64), dimension(:), intent(out) :: g
    real(r64), dimension(:,:), intent(out) :: a

    g = 0.0_r64
    a = 0.0_r64
  end subroutine dummy_grad

  subroutine find_correction(result_x, dev, corr, corr_threshold)
    !! Subroutine for finding the correction matrices M_T enforcing Kelvin-Onsager relations by solving the optimization problem with SLSQP
    !! result_x: output optimization variable containing the correction matrices M_T (reshaped as a vector)
    !! dev: output deviation of the objective function from zero (deviation from Kelvin-Onsager relations)
    !! corr: output maximum deviation of the correction matrices M_T from identity (in %)
    !! corr_threshold: optional input threshold for maximum allowed deviation of M_T from identity (in %), default is 100%
    
    real(r64), intent(out) :: result_x(:)
    real(r64), intent(out), optional :: dev, corr
    real(r64), intent(in), optional :: corr_threshold
    integer :: n, m = 2, meq = 0, i ! Two inequality constraint (positive semidefinite, eigenvalues are real)
    integer :: maxit = 1000, exit_code
    real(r64) :: feastol = 1.0e-6_r64,  c(2), corr_threshold_actual
    real(r64), allocatable :: bl(:), bu(:)
    logical :: status_ok

    n = size(result_x)
    allocate(bl(n), bu(n))

    ! Set bounds for M_T elements
    bl = -100.0_r64
    bu = 100.0_r64

    if (present(corr_threshold)) then
       corr_threshold_actual = corr_threshold
    else
       corr_threshold_actual = 100.0_r64   ! Default value
    end if

    ! Initial guess: identity matrix (matrices)
    result_x = reshape([(reshape(eye(3_i64) * 1.0_r64, [9]), i=1, n/9)], [n])   

    ! Initialize solver with correct parameters
    call solver%initialize(n = n, m = m, meq = meq, max_iter = maxit, acc = feastol, f = fun, &
         g = dummy_grad, gradient_mode = 3, xl = bl, xu = bu, status_ok = status_ok, &
         gradient_delta = 1.0e-4_r64, toldx = 1.0e-6_r64, iprint = 0)

    if(.not. status_ok) then
       exit_code = -1
       return 
    end if

    ! Solve the optimization problem
    call solver%optimize(result_x, exit_code)

    if(exit_code /= 0) then
       if(this_image() == 1) print *, "Warning: KO correction did not converge, skipping correction."
       result_x = reshape([(reshape(eye(3_i64) * 1.0_r64, [9]), i=1, n/9)], [n])
    end if

    !Calculate deviation of the correction matrices from identity, corr = 100*maxval|M_T - I|
    corr = 100.0*maxval(abs(result_x - reshape([(reshape(eye(3_i64) * 1.0_r64, [9]), i=1, n/9)], [n])))
    
    ! Check if correction is bigger then a threshold
    if(corr > corr_threshold_actual) print *, "Warning: KO correction is too large, KO_corr[%] = ", corr, "%  &
            You might better use a finer mesh." 

    !Calculate deviation from identity,  dev = objective function f(x = result_x) * 100  
    call fun(solver, result_x, dev, c)
    dev = 100.0*dev

  end subroutine find_correction

end module optimization