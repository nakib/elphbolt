module optimization
    use slsqp_module, only: slsqp_solver
    use symmetry_module, only: symmetrize_3x3_tensor
    use precision, only: r128, r64, i64
    use misc, only: eye, compute_eigenvalues, F_norm

    implicit none

    public :: find_correction

    type(slsqp_solver) :: solver
    
    type :: coeffs_for_optimization
        real(r64) :: T, el_sigma(3,3), el_kappa0(3,3), el_sigmaS(3,3), el_alpha(3,3), ph_kappa(3,3), &
        ph_alpha(3,3), tot_kappa(3,3), tot_alpha(3,3)
        real(r64), allocatable :: crotations(:,:,:)
    end type coeffs_for_optimization

    type(coeffs_for_optimization) :: global_coeffs, global_coeffs_opp_B

contains

    ! Objective and constraint function for SLSQP
    ! Must match the interface expected by slsqp_solver
    subroutine fun(me, x, f, c)
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

            ! Top-right block: el_sigmaS @ M_T  
            L(1:3, 4:6) = A(1, :, :)
            
            ! Bottom-left block: total_alpha  
            L(4:6, 1:3) = global_coeffs%tot_alpha
            
            ! Bottom-right block: el_kappa0 @ M_T + ph_kappa 
            L(4:6, 4:6) = A(2, :, :) + global_coeffs%ph_kappa 

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
            L(1:3, 4:6) = A(1, :, :)
            
            ! Bottom-left block: el_alpha @ M_J^T + ph_alpha 
            L(4:6, 1:3) = matmul(global_coeffs%el_alpha, M_T(2, :, :)) + global_coeffs%ph_alpha
            
            ! Bottom-right block: el_kappa0 @ M_I^T + ph_kappa @ M_F^T 
            L(4:6, 4:6) = A(2, :, :) + A(4, :, :)
        else 
            print *, "Error: Unexpected size of optimization variable x: ", n
            return
        end if

        ! Compute eigenvalues
        call compute_eigenvalues(L, eigenvalues_real, eigenvalues_imag)

        ! Constraint: min(real(eigenvalue)) >= 0
        c(1) = minval(eigenvalues_real)
        ! Constraint: max(|imag(eigenvalue)|) = 0
        c(2) = 1e-6 - maxval(abs(eigenvalues_imag)) 

        deallocate(M_T, A)

    end subroutine fun

    subroutine dummy_grad(me, x, g, a)
        class(slsqp_solver), intent(inout) :: me
        real(r64), dimension(:), intent(in) :: x
        real(r64), dimension(:), intent(out) :: g
        real(r64), dimension(:,:), intent(out) :: a

        ! This subroutine is never called because gradient_mode=1
        ! Just initialize outputs to avoid uninitialized variable warnings
        g = 0.0_r64
        a = 0.0_r64
    end subroutine dummy_grad

    subroutine find_correction(result_x, dev, corr, corr_threshold)
        ! Input matrices
        ! type(coeffs_for_optimization), intent(in) :: coeffs_local
        real(r64), intent(out) :: result_x(:)
        real(r64), intent(out), optional :: dev, corr
        real(r64), intent(in), optional :: corr_threshold
        ! dev = objective function f(x = result_x) * 100 
        ! corr = 100*maxval|M_T - I|
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
               if(this_image() == 1) then
                print *, "Warning: KO correction did not converge, skipping correction."

               end if 
               result_x = reshape([(reshape(eye(3_i64) * 1.0_r64, [9]), i=1, n/9)], [n])
        end if

        !Calculate deviation from identity
        corr = 100.0*maxval(abs(result_x - reshape([(reshape(eye(3_i64) * 1.0_r64, [9]), i=1, n/9)], [n])))
        if(corr > corr_threshold_actual) then
            print *, "Warning: KO correction is too large, KO_corr[%] = ", corr, "%  &
            You might better use a finer mesh." 

            ! corr = 0.0_r64
            ! result_x = reshape([(reshape(eye(3_i64) * 1.0_r64, [9]), i=1, n/9)], [n])
        end if

        call fun(solver, result_x, dev, c)
        dev = 100.0*dev

    end subroutine find_correction
    
end module optimization