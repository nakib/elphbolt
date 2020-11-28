module wannier
  !! Subroutines for reading Wannier information and Wannier->Bloch transformations

  use params, only: dp, k4 , Ryd2eV, Ryd2radTHz, oneI, twopi, Ryd2amu
  use misc, only: exit_with_message, print_message, expi
  use data, only: numwannbands, numbranches, nwsk, nwsq, nwsg, &
       rcells_k, rcells_q, rcells_g, elwsdeg, phwsdeg, atomtypes, &
       gwsdeg, Hwann, Dphwann, gwann, lattvecs, numatoms, qmesh, masses
  
  implicit none

  public

  !EPW File names
  character(len=*), parameter :: filename_epwdata = "epwdata.fmt"
  character(len=*), parameter :: filename_epwgwann = "epmatwp1"
  character(len=*), parameter :: filename_elwscells = "rcells_k"
  character(len=*), parameter :: filename_phwscells = "rcells_q"
  character(len=*), parameter :: filename_gwscells = "rcells_g"
  character(len=*), parameter :: filename_elwsdeg = "wsdeg_k"
  character(len=*), parameter :: filename_phwsdeg = "wsdeg_q"
  character(len=*), parameter :: filename_gwsdeg = "wsdeg_g"

  !Unit conversion constant
  !real(dp), parameter :: g2unitfactor = Rydberg2eV**3*Rydberg2amu

contains

  subroutine read_EPW_Wannier
    !! Read Wannier representation of the hamiltonian, dynamical matrix, and the
    !! e-ph matrix elements from file epwdata.fmt.

    integer(k4) :: iuc, ib, jb, numbranches_dummy

    open(1,file=filename_epwdata,status='old')
    read(1,*) !ef
    read(1,*) numwannbands, nwsk, numbranches_dummy, nwsq, nwsg
    read(1,*) !zstar, epsil: non-zero only for polar materials

    !Read real space hamiltonian
    call print_message("Reading Wannier rep. Hamiltonian...")
    allocate(Hwann(nwsk,numwannbands,numwannbands))
    do ib = 1,numwannbands
       do jb = 1,numwannbands
          do iuc = 1,nwsk !Number of real space electron cells
             read (1, *) Hwann(iuc,ib,jb)
          end do
       end do
    end do

    !Read real space dynamical matrix
    call print_message("Reading Wannier rep. dynamical matrix...")
    allocate(Dphwann(nwsq,numbranches,numbranches))
    do ib = 1,numbranches
       do jb = 1,numbranches
          do iuc = 1,nwsq !Number of real space phonon cells
             read (1, *) Dphwann(iuc,ib,jb)
          end do
       end do
    end do
    close(1)

    !Read real space matrix elements
    call print_message("Reading Wannier rep. e-ph vertex...")
    open(1, file = filename_epwgwann, status = 'old', access = 'stream')
    allocate(gwann(numwannbands,numwannbands,nwsk,numbranches,nwsg))
    gwann = 0.0_dp
    read(1) gwann
    close(1)

    !Read cell maps of q, k, g meshes.
    call print_message("Reading Wannier cells and multiplicities...")
    allocate(rcells_k(nwsk,3))
    allocate(elwsdeg(nwsk))
    open(1, file = filename_elwscells, status = "old")
    open(2, file = filename_elwsdeg, status = "old")
    do iuc = 1,nwsk
       read(1, *) rcells_k(iuc, :)
       read(2, *) elwsdeg(iuc)
    end do
    close(1)
    close(2)

    allocate(rcells_q(nwsq, 3))
    allocate(phwsdeg(nwsq))
    open(1, file = filename_phwscells, status = "old")
    open(2, file = filename_phwsdeg, status = "old")
    do iuc = 1,nwsq
       read(1, *) rcells_q(iuc, :)
       read(2, *) phwsdeg(iuc)
    end do
    close(1)
    close(2)

    allocate(rcells_g(nwsg, 3))
    allocate(gwsdeg(nwsg))
    open(1, file = filename_gwscells, status = "old")
    open(2, file = filename_gwsdeg, status = "old")
    do iuc = 1,nwsg
       read(1, *) rcells_g(iuc, :)
       read(2, *) gwsdeg(iuc)
    end do
    close(1)
    close(2)
  end subroutine read_EPW_Wannier
  
  subroutine el_wann_epw(nk, kvecs, energies, velocities, evecs)
    !! Wannier interpolate electrons on list of arb. k-vecs
    
    integer(k4), intent(in) :: nk
    real(dp), intent(in) :: kvecs(nk,3) !Crystal coordinates
    real(dp), intent(out) :: energies(nk,numwannbands)
    real(dp), optional, intent(out) :: velocities(nk,numwannbands,3)
    complex(dp), optional, intent(out) :: evecs(nk,numwannbands,numwannbands)

    integer(k4) :: iuc, ib, jb, ipol, ik, nwork, tmp
    real(dp) :: rcart(3)
    real(dp),  allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)
    complex(dp) :: caux, H(numwannbands,numwannbands), &
         dH(3,numwannbands,numwannbands)

    !Catch error for optional velocity calculation
    if(present(velocities) .and. .not. present(evecs)) &
         call exit_with_message("In Wannier, velocity is present but not eigenvecs.")
    
    nwork = 1
    allocate(work(nwork))
    allocate(rwork(max(1,7*numwannbands)))
    
    do ik = 1,nk
       !Form Hamiltonian (H) and k-derivative of H (dH) 
       !from Hwann, rcells_k, and elwsdeg
       H = 0
       dH = 0
       do iuc = 1,nwsk
          caux = expi(twopi*dot_product(kvecs(ik,:),rcells_k(iuc,:)))&
               /elwsdeg(iuc)
          H = H + caux*Hwann(iuc,:,:)

          if(present(velocities)) then
             rcart = matmul(lattvecs,rcells_k(iuc,:))
             do ipol = 1,3
                dH(ipol,:,:) = dH(ipol,:,:) + &
                     oneI*rcart(ipol)*caux*Hwann(iuc,:,:)
             end do
          end if
       end do

       !Force Hermiticity
       do ib = 1,numwannbands
          do jb = ib + 1,numwannbands
             H(ib,jb) = (H(ib,jb) + conjg(H(jb,ib)))*0.5_dp
             H(jb,ib) = H(ib,jb)
          end do
       end do

       !Diagonalize H
       call zheev("V", "U", numwannbands, H(:,:), numwannbands, energies(ik,:), &
            work, -1, rwork, tmp)
       if(real(work(1)) > nwork) then
          nwork = nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V", "U", numwannbands, H(:,:), numwannbands, energies(ik,:), &
            work, nwork, rwork, tmp)

       if(present(evecs)) then
          evecs(ik,:,:)=transpose(H(:,:))
       end if

       if(present(velocities)) then
          !Calculate velocities using Feynman-Hellmann thm
          do ib = 1,numwannbands
             do ipol = 1,3
                velocities(ik,ib,ipol)=real(dot_product(evecs(ik,ib,:), &
                     matmul(dH(ipol,:,:), evecs(ik,ib,:))))
             end do
          end do
       end if

       !energies(ik,:) = energies(ik,:)*Rydberg2radTHz !2piTHz
       energies(ik,:) = energies(ik,:)*Ryd2eV !eV

       if(present(velocities)) then
          velocities(ik,:,:) = velocities(ik,:,:)*Ryd2radTHz !nmTHz = Km/s
       end if
    end do !ik
  end subroutine el_wann_epw

  subroutine ph_wann_epw(nq, qvecs, energies, velocities, evecs)
    !! Wannier interpolate phonons on list of arb. q-vec
  
    integer(k4), intent(in) :: nq
    real(dp), intent(in) :: qvecs(nq, 3) !Crystal coordinates
    real(dp), intent(out) :: energies(nq, numbranches)
    real(dp), intent(out), optional :: velocities(nq, numbranches, 3)
    complex(dp), intent(out), optional :: evecs(nq, numbranches, numbranches)

    integer(k4) :: iuc, ib, jb, ipol, iq, na, nb, nwork, aux
    complex(dp) :: caux
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)
    real(dp) :: omega2(numbranches), rcart(3), massnorm
    complex(dp) :: dynmat(numbranches, numbranches), ddynmat(3, numbranches, numbranches)

    !Catch error for optional velocity calculation
    if(present(velocities) .and. .not. present(evecs)) &
         call exit_with_message("In Wannier, velocity is present but not eigenvecs.")
    
    nwork = 1
    allocate(work(nwork))
    allocate(rwork(max(1, 9*numatoms-2)))
    
    do iq = 1, nq
       !Form dynamical matrix (dynmat) and q-derivative of dynmat (ddynmat) 
       !from Dphwann, rcells_q, and phwsdeg
       dynmat = 0
       ddynmat = 0
       do iuc = 1, nwsq
          caux = expi(twopi*dot_product(qvecs(iq, :), rcells_q(iuc, :)))&
               /phwsdeg(iuc)
          dynmat = dynmat + caux*Dphwann(iuc, :, :)

          if(present(velocities)) then
             rcart = matmul(lattvecs, rcells_q(iuc, :))
             do ipol=1, 3
                ddynmat(ipol, :, :) = ddynmat(ipol, :, :) + &
                     oneI*rcart(ipol)*caux*Dphwann(iuc, :, :)
             end do
          end if
       end do

       !Non-analytic correction
       !if(polar) then
       !   call dyn_nonanalytic(matmul(rlattvec,qvecs(iq,:))*bohr2nm, dynmat, ddynmat)
       !end if

       !Force Hermiticity
       do ib = 1, numbranches
          do jb = ib + 1, numbranches
             dynmat(ib, jb) = (dynmat(ib, jb) + conjg(dynmat(jb, ib)))*0.5_dp
             dynmat(jb, ib) = dynmat(ib, jb)
          end do
       end do
       
       !Mass normalize
       do na = 1, numatoms
          do nb = 1, numatoms
             massnorm = 1.d0/sqrt(masses(atomtypes(na))*masses(atomtypes(nb)))*Ryd2amu
             dynmat(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
                  dynmat(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb)*massnorm
             if(present(velocities)) then
                do ipol=1, 3
                   ddynmat(ipol, 3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
                        ddynmat(ipol, 3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb)*massnorm
                end do
             end if
          end do
       end do
       
       !Diagonalize dynmat
       call zheev("V", "U", numbranches, dynmat(:, :), numbranches, omega2, work, -1, rwork, aux)
       if(real(work(1)) > nwork) then
          nwork = nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V", "U", numbranches, dynmat(:, :), numbranches, omega2, work, nwork, rwork, aux)
       
       energies(iq, :) = sign(sqrt(abs(omega2)), omega2)
       if(present(evecs)) then
          evecs(iq, :, :) = transpose(dynmat(:, :))
       end if

       if(present(velocities)) then
          !Calculate velocities using Feynman-Hellmann thm
          do ib = 1, numbranches
             do ipol = 1, 3
                velocities(iq, ib, ipol) = real(dot_product(dynmat(:, ib), &
                     matmul(ddynmat(ipol, :, :), dynmat(:, ib))))
             end do
             velocities(iq, ib, :) = velocities(iq, ib, :)/(2.0_dp*energies(iq, ib))
          end do
       end if

       !energies(iq, :) = energies(iq, :)*Rydberg2radTHz !2piTHz
       !energies(iq, :) = energies(iq, :)*Rydberg2eV*1.0e3_dp !meV
       energies(iq, :) = energies(iq, :)*Ryd2eV !eV
       if(present(velocities)) then
          velocities(iq, :, :) = velocities(iq, :, :)*Ryd2radTHz !nmTHz = Km/s
       end if
       
       !Take care of gamma point.
       if(nq == product(qmesh) .or. all(qvecs(1,:) == 0)) then
          energies(1, 1:3) = 0
          if(present(velocities)) then
             velocities(1, :, :) = 0
          end if
       end if

       !Handle negative energy phonons
       do ib = 1, numbranches
          if(energies(iq,ib) < -0.005_dp) then
             print*, 'Oh no!'
             call exit_with_message('Large negative phonon energy found! Stopping!')
             
          else if(energies(iq,ib) < 0 .and. energies(iq,ib) > -0.005_dp) then
             print*, 'Warning! Small negative phonon energy for q-point, branch =', &
                  iq, ib
             print*, 'Setting it to zero. Please check Wannier parameters.'
             energies(iq,ib) = 0
          end if
       end do
    end do !iq
  end subroutine ph_wann_epw

!!$  subroutine dyn_nonanalytic(q, dyn, ddyn)
!!$    !! Calculate the long-range correction to the
!!$    !! dynamical matrix and its derivative for a given phonon mode.
!!$    !!
!!$    !! q is the phonon wvec in Cartesian coords., Bohr^-1
!!$    !! (d)dyn is the (derivative of)dynamical matrix
!!$    
!!$    real(dp), intent(in) :: q(3) !Cartesian
!!$    complex(dp), intent(inout) :: dyn(nbranches,nbranches)
!!$    complex(dp), intent(inout) :: ddyn(3,nbranches,nbranches)
!!$
!!$    complex(dp) :: dyn_l(nbranches,nbranches)
!!$    complex(dp) :: ddyn_l(3,nbranches,nbranches)
!!$    real(dp) :: qeq,     &! <q+g| epsilon |q+g>
!!$         arg, zag(3),zbg(3), g(3), gmax, alph, geg,&
!!$         tpiba,dgeg(3),fnat(3),rr(natoms,natoms,3)
!!$    integer(kind=4) :: na,nb,i,idim,jdim,ipol,jpol,m1,m2,m3,nq1,nq2,nq3
!!$    complex(dp) :: fac, facqd, facq
!!$
!!$    tpiba = twopi/normtwo(lattvec(:,1))*bohr2nm
!!$
!!$    !Recall that the phonon supercell in elphBolt is the
!!$    !same as the EPW coarse phonon mesh.
!!$    nq1 = scell(1)
!!$    nq2 = scell(2)
!!$    nq3 = scell(3)
!!$
!!$    gmax= 14.d0 !dimensionless
!!$    alph= tpiba**2 !bohr^-2
!!$    geg = gmax*alph*4.0d0
!!$    !In Ry units, qe = sqrt(2.0)
!!$    fac = 8.d0*pi/(V/bohr2nm**3)
!!$
!!$    dyn_l = 0.d0
!!$    ddyn_l = 0.d0
!!$    do m1 = -nq1,nq1
!!$       do m2 = -nq2,nq2
!!$          do m3 = -nq3,nq3
!!$             g(:) = (m1*rlattvec(:,1)+m2*rlattvec(:,2)+m3*rlattvec(:,3))*bohr2nm
!!$             qeq = dot_product(g,matmul(epsilon,g))
!!$             
!!$             if (qeq > 0.d0 .and. qeq/alph/4.d0 < gmax ) then
!!$                facqd = exp(-qeq/alph/4.0d0)/qeq
!!$
!!$                do na = 1,natoms
!!$                   zag(:)=matmul(g,born(:,:,na))
!!$                   fnat(:)=0.d0
!!$                   do nb = 1,natoms
!!$                      rr(na,nb,:) = (cartesian(:,na)-cartesian(:,nb))/bohr2nm
!!$                      arg = dot_product(g,rr(na,nb,:))
!!$                      zbg(:) = matmul(g,born(:,:,nb))
!!$                      fnat(:) = fnat(:) + zbg(:)*phexp(arg)
!!$                   end do
!!$                   do jpol=1,3
!!$                      jdim=(na-1)*3+jpol
!!$                      do ipol=1,3
!!$                         idim=(na-1)*3+ipol
!!$                         dyn_l(idim,jdim) = dyn_l(idim,jdim) - facqd * &
!!$                              zag(ipol) * fnat(jpol)
!!$                      end do
!!$                   end do
!!$                end do
!!$             end if
!!$
!!$             g = g + q
!!$             qeq = dot_product(g,matmul(epsilon,g))
!!$             if (qeq > 0.d0 .and. qeq/alph/4.d0 < gmax ) then
!!$                facqd = exp(-qeq/alph/4.0d0)/qeq
!!$                dgeg=matmul(epsilon+transpose(epsilon),g)
!!$                do nb = 1,natoms
!!$                   zbg(:)=matmul(g,born(:,:,nb))
!!$                   do na = 1,natoms
!!$                      rr(na,nb,:) = (cartesian(:,na)-cartesian(:,nb))/bohr2nm
!!$                      zag(:)=matmul(g,born(:,:,na))
!!$                      arg = dot_product(g,rr(na,nb,:))
!!$                      facq = facqd*phexp(arg)
!!$                      do jpol=1,3
!!$                         jdim=(nb-1)*3+jpol
!!$                         do ipol=1,3
!!$                            idim=(na-1)*3+ipol
!!$                            dyn_l(idim,jdim) = dyn_l(idim,jdim) + facq * &
!!$                                 zag(ipol)*zbg(jpol)
!!$                            !Correction to derivative of dynmat
!!$                            ddyn_l(:,idim,jdim)=ddyn_l(:,idim,jdim)+&
!!$                                 facq*&
!!$                                 ( zbg(jpol)*born(:,ipol,na)+zag(ipol)*born(:,jpol,nb)+&
!!$                                 zag(ipol)*zbg(jpol)*(iunit*rr(na,nb,:)-&
!!$                                 dgeg(:)/alph/4.0-dgeg(:)/qeq) )
!!$                         end do
!!$                      end do
!!$                   end do
!!$                end do
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$    dyn = dyn + dyn_l*fac
!!$    ddyn = ddyn + ddyn_l*fac
!!$  end subroutine dyn_nonanalytic
  
!!$
!!$  function g2_epw(qvec, el_evec_k, el_evec_kp, ph_evec_q, ph_en, gmixed_ik)
!!$    !! Function to calculate |g|^2.
!!$    !! This works with EPW real space data
!!$    !! qvec: phonon wave vector in crystal coords
!!$    !! el_evec_k(kp): initial(final) electron eigenvector in bands m(n) 
!!$    !! ph_evec_q: phonon eigenvector branch s 
!!$    !! ph_en: phonon energy in mode (s,qvec)
!!$    !! gmixed_ik: e-ph matrix element
!!$    use configuration, only: numbranches, numwannbands, nwsq, nwsg, masses, & 
!!$         atomtypes, rcells_g, gwsdeg
!!$
!!$    implicit none
!!$
!!$    real(dp),intent(in) :: qvec(3), ph_en
!!$    complex(dp),intent(in) :: el_evec_k(numwannbands),&
!!$         el_evec_kp(numwannbands), ph_evec_q(numbranches), &
!!$         gmixed_ik(numwannbands, numwannbands, numbranches, nwsq)
!!$    integer(kind=4) :: ip, ig, np, mp, sp, mtype
!!$    complex(dp) :: caux, u(numbranches), UkpgkUk(numbranches, nwsq), &
!!$         UkpgkUkuq(nwsq), gbloch, overlap(numwannbands,numwannbands)
!!$    real(dp) :: g2_epw
!!$
!!$    !Mass normalize the phonon matrix
!!$    do ip = 1, numbranches ! d.o.f of basis atoms
!!$       !demux atom type from d.o.f
!!$       mtype = (ip - 1)/3 + 1 
!!$       !normalize and conjugate eigenvector
!!$       u(ip) = ph_evec_q(ip)/sqrt(masses(atomtypes(mtype)))
!!$    end do
!!$    
!!$    if(ph_en == 0) then !zero out matrix elements at Gamma
!!$       g2_epw = 0
!!$    else
!!$       UkpgkUk = 0 !g(k,Rp) rotated by the electron U(k), U(k') matrices
!!$       UkpgkUkuq = 0 !above quantity rotated by the phonon u(q) matrix
!!$       gbloch = 0
!!$
!!$       !Copy data to accelerator device
!!$       !$acc data copyin(UkpgkUk(:,:), UkpgkUkuq(:), el_evec_k(:), el_evec_kp(:), u(:), overlap(:,:))
!!$       !$acc kernels present(gmixed_ik(:,:,:,:), gwsdeg(:), rcells_g(:,:))
!!$
!!$       !Create the <n'|m'> overlap matrix
!!$       do np = 1, numwannbands !over final electron band
!!$          do mp = 1, numwannbands !over initial electron band
!!$             overlap(mp,np) = conjg(el_evec_kp(np))*el_evec_k(mp)
!!$          end do
!!$       end do
!!$
!!$       do ig = 1, nwsg !over matrix elements WS cell
!!$          !Apply electron rotations
!!$          do sp = 1, numbranches
!!$             caux = 0
!!$             do np = 1, numwannbands !over final electron band
!!$                do mp = 1, numwannbands !over initial electron band
!!$                   caux = caux + overlap(mp,np)*gmixed_ik(np, mp, sp, ig)
!!$                end do
!!$             end do
!!$             UkpgkUk(sp, ig) = UkpgkUk(sp, ig) + caux
!!$          end do
!!$       end do
!!$       
!!$       do ig = 1, nwsg !over matrix elements WS cell
!!$          !Apply phonon rotation
!!$          UkpgkUkuq(ig) = UkpgkUkuq(ig) + dot_product(conjg(u),UkpgkUk(:, ig))
!!$       end do
!!$
!!$       do ig = 1, nwsg !over matrix elements WS cell
!!$          !Fourier transform to q-space
!!$          caux = exp(twopiI*dot_product(qvec, rcells_g(ig, :)))&
!!$               /gwsdeg(ig)
!!$          gbloch = gbloch + caux*UkpgkUkuq(ig)
!!$       end do
!!$       
!!$       !$acc end kernels
!!$       !$acc end data
!!$       
!!$       g2_epw = 0.5_dp*real(gbloch*conjg(gbloch))/ &
!!$            ph_en*g2unitfactor !eV^2
!!$    end if
!!$  end function g2_epw
!!$
!!$  subroutine calc_all_gk2
!!$    !! MPI parallelizer of g2_epw over IBZ electron states within the Fermi window.
!!$    !In the FBZ and IBZ blocks a wave vector was retained when at least one
!!$    !band belonged within the energy window. Here the bands outside energy window
!!$    !will be skipped in the calculation as they are irrelevant for transport.
!!$    !This subroutine will calculate the full Bloch rep. matrix elements for
!!$    !all the energy window restricted electron-phonon processes for a given
!!$    !irreducible initial electron state = (band, wave vector). 
!!$    !This list will be written to disk in files tagged with the muxed state index. 
!!$    use configuration, only: nk, nk_irred, procid, numwannbands, &
!!$         numbranches, el_ens_irred, enref, fsthick, kmesh, el_wavevecs_irred, &
!!$         el_wavevecs, el_ens, el_evecs, el_evecs_irred, ph_evecs, ph_ens, gdir, &
!!$         cwd, nwsq, el_indexlist_irred, nstates_inwindow, gwsdeg, rcells_g
!!$
!!$    implicit none
!!$
!!$    integer(k4) :: nstates_irred, istate, m, ik, ik_muxed, n, ikp, s, &
!!$         iq, start, end, chunk, ierr, k_indvec(3), kp_indvec(3), &
!!$         q_indvec(3), count, g2size
!!$    real(dp) :: k(3), kp(3), q(3)
!!$    real(dp), allocatable :: g2_istate(:)
!!$    complex(dp) :: gmixed_ik(numwannbands, numwannbands, numbranches, nwsq)
!!$    character(len = 1024) :: filename
!!$
!!$    call print_message("Doing g(k,Rp) -> |g(k,q)|^2 for all IBZ states...")
!!$
!!$    !Length of g2_istate
!!$    g2size = nstates_inwindow*numbranches
!!$    allocate(g2_istate(g2size))
!!$
!!$    !Total number of IBZ blocks states
!!$    nstates_irred = nk_irred*numwannbands
!!$    
!!$    call mpi_chunk(nstates_irred, chunk, start, end)
!!$
!!$    if(procid == 0) then
!!$       print*, "   #states = ", nstates_irred
!!$       print*, "   #states/process = ", chunk
!!$    end if
!!$
!!$    !Copy gwsdeg to accelerator
!!$    !$acc data copyin(gwsdeg(:),rcells_g(:,:))
!!$
!!$    count = 0
!!$    do istate = start, end !over IBZ blocks states
!!$       !Initialize eligible process counter for this state
!!$       count = 0
!!$
!!$       !Demux state index into band (m) and wave vector (ik) indices
!!$       call demux_state(istate, numwannbands, m, ik)
!!$       
!!$       !Get the muxed index of wave vector from the IBZ blocks index list
!!$       ik_muxed = el_indexlist_irred(ik)
!!$
!!$       !Apply energy window to initial (IBZ blocks) electron
!!$       if(abs(el_ens_irred(ik, m) - enref) > fsthick) cycle
!!$
!!$       !Load gmixed(ik) here for use inside the loops below
!!$       call chdir(trim(adjustl(gdir)))
!!$       write (filename, '(I6)') ik
!!$       filename = 'gmixed.ik'//trim(adjustl(filename))
!!$       open(1,file=filename,status="old",access='stream')
!!$       read(1) gmixed_ik
!!$       close(1)
!!$       call chdir(cwd)
!!$
!!$       !Copy gmixed.ik to accelerator
!!$       !$acc data copyin(gmixed_ik(:,:,:,:))
!!$       
!!$       !Initial (IBZ blocks) wave vector (crystal coords.)
!!$       k = el_wavevecs_irred(ik, :)
!!$
!!$       !Convert from crystal to 0-based index vector
!!$       k_indvec = nint(k*kmesh)
!!$       
!!$       !Run over final (FBZ blocks) electron wave vectors
!!$       do ikp = 1, nk
!!$          !Final wave vector (crystal coords.)
!!$          kp = el_wavevecs(ikp, :)
!!$
!!$          !Convert from crystal to 0-based index vector
!!$          kp_indvec = nint(kp*kmesh)
!!$          
!!$          !Run over final electron bands
!!$          do n = 1, numwannbands
!!$             !Apply energy window to final electron
!!$             if(abs(el_ens(ikp, n) - enref) > fsthick) cycle
!!$
!!$             !Find interacting phonon wave vector
!!$             !Note that q, k, and k' are all on the same mesh
!!$             q_indvec = modulo(kp_indvec - k_indvec, kmesh) !0-based index vector
!!$             q = q_indvec/dble(kmesh) !crystal coords.
!!$
!!$             !Muxed index of q
!!$             iq = mux_vector(q_indvec, kmesh, 0)
!!$             
!!$             !Run over phonon branches
!!$             do s = 1, numbranches
!!$                !Increment g2 processes counter
!!$                count = count + 1
!!$
!!$                !Calculate |g_mns(<k>,q)|^2
!!$                g2_istate(count) = g2_epw(q, el_evecs_irred(ik, m, :), &
!!$                     el_evecs(ikp, n, :), ph_evecs(iq, s, :), &
!!$                     ph_ens(iq, s), gmixed_ik)
!!$             end do !s
!!$          end do !n
!!$       end do !ikp
!!$
!!$       !$acc end data
!!$
!!$       !Change to data output directory
!!$       call chdir(trim(adjustl(gdir)))
!!$
!!$       !Write data in binary format
!!$       !Note: this will overwrite existing data!
!!$       write (filename, '(I6)') istate
!!$       filename = 'g2.istate'//trim(adjustl(filename))
!!$       open(1, file = trim(filename), status = 'replace', access = 'stream')
!!$       write(1) g2_istate
!!$       close(1)
!!$
!!$       !Change back to working directory
!!$       call chdir(cwd)
!!$    end do
!!$
!!$    !$acc end data
!!$
!!$    !TODO at this point we can delete the gmixed disk data
!!$    !call clear_gmixed_cache
!!$
!!$    call mpi_barrier(mpi_comm_world, ierr)
!!$  end subroutine calc_all_gk2
!!$
!!$  subroutine gmixed_epw(ik)
!!$    !! Calculate the Bloch-Wannier mixed rep. e-ph matrix elements g(k,Rp),
!!$    !! where k is an IBZ electron wave vector and Rp is a phonon unit cell.
!!$    !! Note: this step *DOES NOT* perform the rotation over the Wannier bands space.
!!$    !The result will be saved to disk.
!!$    use configuration, only: numwannbands, numbranches, nwsk, nwsq, &
!!$         rcells_k, elwsdeg, el_wavevecs_irred, gdir, gwann, cwd
!!$
!!$    implicit none
!!$    
!!$    integer(k4), intent(in) :: ik
!!$    integer(k4) :: iuc
!!$    complex(dp) :: caux
!!$    complex(dp), allocatable:: gmixed(:, :, :, :)
!!$    real(dp) :: kvec(3)
!!$    character(len = 1024) :: filename
!!$
!!$    allocate(gmixed(numwannbands, numwannbands, numbranches, nwsq))
!!$    
!!$    !Electron wave vector (crystal coords.) in IBZ blocks 
!!$    kvec = el_wavevecs_irred(ik, :)
!!$
!!$    !Fourier transform to k-space
!!$    gmixed = 0
!!$    do iuc = 1, nwsk
!!$       caux = expi(twopi*dot_product(kvec, rcells_k(iuc, :))) &
!!$            /elwsdeg(iuc)
!!$       gmixed(:, :, :, :) = gmixed(:, :, :, :) + &
!!$            caux*gwann(:, :, iuc, :, :)
!!$    end do
!!$
!!$    !Change to data output directory
!!$    call chdir(trim(adjustl(gdir)))
!!$
!!$    !Write data in binary format
!!$    !Note: this will overwrite existing data!
!!$    write (filename, '(I6)') ik
!!$    filename = 'gmixed.ik'//trim(adjustl(filename))
!!$    open(1, file = trim(filename), status = 'replace', access = 'stream')
!!$    write(1) gmixed
!!$    close(1)
!!$
!!$    !Change back to working directory
!!$    call chdir(cwd)
!!$  end subroutine gmixed_epw
!!$
!!$  subroutine gmixed_epw_driver
!!$    !! MPI parallizer of gmixed_epw over IBZ electron wave vectors
!!$    use configuration, only: nk_irred
!!$
!!$    implicit none
!!$
!!$    integer(k4) :: ik, ikstart, ikend, chunk, ierr
!!$
!!$    call print_message("Doing g(Re,Rp) -> g(k,Rp) for all IBZ k...")
!!$
!!$    call mpi_chunk(nk_irred, chunk, ikstart, ikend)
!!$
!!$    if(procid == 0) then
!!$       print*, "   #k = ", nk_irred
!!$       print*, "   #k/process = ", chunk
!!$    end if
!!$
!!$    do ik = ikstart, ikend
!!$       call gmixed_epw(ik)
!!$    end do
!!$
!!$    call mpi_barrier(mpi_comm_world, ierr)
!!$  end subroutine gmixed_epw_driver
!!$
!!$  !For testing and debugging:
!!$
!!$  !Calculate the Bloch-Wannier mixed rep. e-ph matrix elements g(k,Rp),
!!$  !where k is an IBZ electron wave vector and Rp is a phonon unit cell.
!!$  !This is for the electron at gamma.
!!$  subroutine gmixed_epw_gamma
!!$    use configuration, only: numwannbands, numbranches, nwsk, nwsq, &
!!$         elwsdeg, gdir, gwann, cwd
!!$
!!$    implicit none
!!$    
!!$    integer(k4) :: iuc
!!$    complex(dp) :: caux
!!$    complex(dp), allocatable:: gmixed(:, :, :, :)
!!$    character(len = 1024) :: filename
!!$
!!$    allocate(gmixed(numwannbands, numwannbands, numbranches, nwsq))
!!$    
!!$    !Fourier transform to k-space
!!$    gmixed = 0
!!$    do iuc = 1, nwsk
!!$       caux = (1.0_dp,0.0_dp)/elwsdeg(iuc)           
!!$       gmixed(:, :, :, :) = gmixed(:, :, :, :) + &
!!$            caux*gwann(:, :, iuc, :, :)
!!$    end do
!!$
!!$    !Change to data output directory
!!$    call chdir(trim(adjustl(gdir)))
!!$
!!$    !Write data in binary format
!!$    !Note: this will overwrite existing data!
!!$    filename = 'gmixed.gamma'//trim(adjustl(filename))
!!$    open(1, file = trim(filename), status = 'replace', access = 'stream')
!!$    write(1) gmixed
!!$    close(1)
!!$
!!$    !Change back to working directory
!!$    call chdir(cwd)
!!$  end subroutine gmixed_epw_gamma
end module wannier
