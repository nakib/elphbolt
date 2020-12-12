program elphBolt
  !! Author: Nakib H. Protik
  !! Summary: Main driver program.
  !!
  !! elphBolt is a program for solving the coupled electron-phonon Boltzmann
  !! transport equations (e-ph BTEs) as formulated in Phys. Rev. B 101, 075202 (2020)
  !! and arXiv:2008.08722 (2020) with both the electron-phonon and phonon-phonon
  !! interactions computed ab initio.

  use config, only: initialize_system
  use mesh, only: calculate_electrons, calculate_phonons
  use wannier, only: calculate_g_mixed, calculate_g_bloch, gmixed_epw_gamma
  
  implicit none
  
  if(this_image() == 1) then
     print*, 'Number of images = ', num_images()
  end if

  !Read inputs, find symmetry operations, and set up calculation environment.
  call initialize_system

  !TODO Test electron bands, phonon dispersions, and g along path.
  if(this_image() == 1) then
     call test_wannier
  end if
  
  !Calculate electrons
  call calculate_electrons

  !Calculate phonons
  call calculate_phonons

  !Calculate mixed Bloch-Wannier space e-ph vertex
  call calculate_g_mixed

  !Calculate Bloch space e-ph vertex
  call calculate_g_bloch

  !Calculate ph-ph vertex
  !TODO call calculate_V

  !Calculate transition probabilities

  !Iterate BTEs

  contains

  !Unit tester for Wannier interpolation with EPW inputs.
  subroutine test_wannier
    use wannier, only: ph_wann_epw, el_wann_epw, g2_epw, gmixed_epw_gamma
    use params
    use data

    integer(k4) :: i, nqpath, m, n, s
    real(dp) :: k(3), kp(3)
    real(dp), allocatable :: qpathvecs(:,:), ph_ens_path(:,:), &
         el_ens_path(:,:), ph_vels_path(:,:,:), el_ens_kp(:), &
         el_vels_kp(:,:), g2_qpath(:), el_ens_k(:), el_vels_k(:,:)
    complex(dp), allocatable :: ph_evecs_path(:,:,:), el_evecs_kp(:,:), &
         el_evecs_k(:,:), gmixed_gamma(:,:,:,:) 
    character(len = 1024) :: filename
    character(len=8) :: saux

    !Read list of wavevectors in crystal coordinates
    nqpath = 601
    allocate(qpathvecs(nqpath,3))
    open(1,file=trim('highsymqpts.txt'),status='old')
    read(1,*) !skip first line
    do i = 1, nqpath
       read(1,*) qpathvecs(i,:)
    end do

    !Test phonon Wannier interpolation
    print*, 'Doing phonon Wannier test...'

    !Calculate phonon dispersions
    allocate(ph_ens_path(nqpath,numbranches), ph_vels_path(nqpath,numbranches,3),&
         ph_evecs_path(nqpath,numbranches,numbranches))
    call ph_wann_epw(nqpath, qpathvecs, ph_ens_path, ph_vels_path, ph_evecs_path)

    !Output phonon dispersions
    write(saux,"(I0)") numbranches
    open(1,file="ph_ens_path",status="replace")
    do i = 1, nqpath
       write(1,"("//trim(adjustl(saux))//"E20.10)") ph_ens_path(i,:)
    end do
    close(1)

    !Test electron Wannier interpolation
    print*, 'Doing electron Wannier test...'

    !Calculate electron bands
    allocate(el_ens_path(nqpath,numwannbands))
    call el_wann_epw(nqpath, qpathvecs, el_ens_path)

    !Output electron dispersions
    write(saux,"(I0)") numwannbands
    open(1,file="el_ens_path",status="replace")
    do i=1,nqpath
       write(1,"("//trim(adjustl(saux))//"E20.10)") el_ens_path(i,:)
    end do
    close(1)

    !Test matrix elements Wannier interpolations
    print*, 'Doing g(k,q) Wannier test...'

    !Initial electron at Gamma
    k = (/0.0_dp, 0.0_dp, 0.0_dp/)
    allocate(el_ens_k(numwannbands), el_vels_k(numwannbands,3),&
         el_evecs_k(numwannbands,numwannbands))
    call el_wann_epw(1_k4, k, el_ens_k, el_vels_k, el_evecs_k)

    !Calculate g(k, Rp)
    call gmixed_epw_gamma

    !Load gmixed from file
    !Change to data output directory
    call chdir(trim(adjustl(g2dir)))
    allocate(gmixed_gamma(numwannbands, numwannbands, numbranches, nwsq))
    filename = 'gmixed.gamma' !//trim(adjustl(filename))
    open(1,file=filename,status="old",access='stream')
    read(1) gmixed_gamma
    close(1)
    !Change back to working directory
    call chdir(cwd)

    !print*, 'Rydberg2eV = ', Ryd2eV
    !print*, 'Rydberg2amu = ', Ryd2amu

    !Calculate g(k, k') where k' = (k + q) modulo G
    !q is along the path loaded from the file
    !phonon branch, s = 3 (LA)
    !initial (final) electron band, m(n) = 1(1) 
    m = 1
    n = 2
    s = 3
    allocate(el_ens_kp(numwannbands), el_vels_kp(numwannbands,3),&
         el_evecs_kp(numwannbands,numwannbands))
    allocate(g2_qpath(nqpath))
    do i = 1, nqpath
       kp = qpathvecs(i, :) !for k = Gamma 

       !Calculate electrons at path point
       call el_wann_epw(1_k4, kp, el_ens_kp, el_vels_kp, el_evecs_kp)

       !Calculate |g(k,k')|^2
       !(q is the same as k')
       g2_qpath(i) = g2_epw(kp, el_evecs_k(m, :), &
            el_evecs_kp(n, :), ph_evecs_path(i, s, :), ph_ens_path(i, s), &
            gmixed_gamma)

       !print*, el_ens_kp(n), ph_ens_path(i,s)*1e3, sqrt(g2_qpath(i))*1.0d3
    end do

    print*, 'done calculating g'
    open(1,file='g_qpath_123',status="replace")
    write(1,*) '#|g_SE| [eV]'
    do i=1,nqpath
       write(1,"(E20.10)") sqrt(g2_qpath(i))
    end do
    close(1)

  end subroutine test_wannier
end program elphBolt
