module wannier
  !! Subroutines for reading Wannier information and Wannier->Bloch transformations

  use params, only: dp, k4 , Ryd2eV, Ryd2radTHz, oneI, pi, twopi, twopiI, &
       Ryd2amu, bohr2nm
  use misc, only: exit_with_message, print_message, expi, twonorm, &
       distribute_points, demux_state, mux_vector
  use data, only: atomtypes, lattvecs, numatoms, qmesh, masses, &
       reclattvecs, born, epsilon, basis_cart, volume, polar, nk_irred, &
       el_wavevecs_irred, cwd, g2dir, el_wavevecs, enref, fsthick, kmesh, nk, &
       nstates_inwindow, el_indexlist_irred, el_ens_irred, el_ens, el_evecs, &
       el_evecs_irred, ph_evecs, ph_ens

  implicit none

  private
  public epw_wannier 

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
  real(dp), parameter :: g2unitfactor = Ryd2eV**3*Ryd2amu

  type epw_wannier
     !! Data and procedures related to Wannierization.

     integer(k4) :: numwannbands
     !! Number of Wannier bands.
     integer(k4) :: numbranches
     !! Number of phonon branches.
     integer(k4) :: nwsk
     !! Number of real space cells for electrons.
     integer(k4) :: nwsq
     !! Number of real space cells for phonons.
     integer(k4) :: nwsg
     !! Number of real space cells for electron-phonon vertex.
     integer(k4), allocatable :: rcells_k(:, :)
     !! Real space cell locations for electrons.
     integer(k4), allocatable :: rcells_q(:, :)
     !! Real space cell locations for phonons.
     integer(k4), allocatable :: rcells_g(:, :)
     !! Real space cell locations for electron-phonon vertex.     
     integer(k4), allocatable :: elwsdeg(:)
     !! Real space cell multiplicity for electrons.
     integer(k4), allocatable :: phwsdeg(:)
     !! Real space cell multiplicity for phonons.
     integer(k4), allocatable :: gwsdeg(:)
     !! Real space cell multiplicity for electron-phonon vertex.
     complex(dp), allocatable :: Hwann(:, :, :)
     !! Hamiltonian in Wannier representation.
     complex(dp), allocatable :: Dphwann(:, :, :)
     !! Dynamical matrix in Wannier representation.
     complex(dp), allocatable :: gwann(:, :, :, :, :)
     !! e-ph vertex in Wannier representation.

   contains

     procedure :: read=>read_EPW_Wannier, el_wann_epw, ph_wann_epw, calculate_g_bloch, &
          gmixed_epw, gmixed_epw_gamma, g2_epw, test_wannier

  end type epw_wannier

contains

  subroutine read_EPW_Wannier(w)
    !! Read Wannier representation of the hamiltonian, dynamical matrix, and the
    !! e-ph matrix elements from file epwdata.fmt.

    class(epw_wannier), intent(out) :: w
    integer(k4) :: iuc, ib, jb !, numbranches

    open(1,file=filename_epwdata,status='old')
    read(1,*) !ef
    read(1,*) w%numwannbands, w%nwsk, w%numbranches, w%nwsq, w%nwsg
    read(1,*) !zstar, epsil: non-zero only for polar materials

    !Read real space hamiltonian
    call print_message("Reading Wannier rep. Hamiltonian...")
    allocate(w%Hwann(w%nwsk,w%numwannbands,w%numwannbands))
    do ib = 1,w%numwannbands
       do jb = 1,w%numwannbands
          do iuc = 1,w%nwsk !Number of real space electron cells
             read (1, *) w%Hwann(iuc,ib,jb)
          end do
       end do
    end do

    !Read real space dynamical matrix
    call print_message("Reading Wannier rep. dynamical matrix...")
    allocate(w%Dphwann(w%nwsq,w%numbranches,w%numbranches))
    do ib = 1,w%numbranches
       do jb = 1,w%numbranches
          do iuc = 1,w%nwsq !Number of real space phonon cells
             read (1, *) w%Dphwann(iuc,ib,jb)
          end do
       end do
    end do
    close(1)

    !Read real space matrix elements
    call print_message("Reading Wannier rep. e-ph vertex...")
    open(1, file = filename_epwgwann, status = 'old', access = 'stream')
    allocate(w%gwann(w%numwannbands,w%numwannbands,w%nwsk,w%numbranches,w%nwsg))
    w%gwann = 0.0_dp
    read(1) w%gwann
    close(1)

    !Read cell maps of q, k, g meshes.
    call print_message("Reading Wannier cells and multiplicities...")
    allocate(w%rcells_k(w%nwsk,3))
    allocate(w%elwsdeg(w%nwsk))
    open(1, file = filename_elwscells, status = "old")
    open(2, file = filename_elwsdeg, status = "old")
    do iuc = 1,w%nwsk
       read(1, *) w%rcells_k(iuc, :)
       read(2, *) w%elwsdeg(iuc)
    end do
    close(1)
    close(2)

    allocate(w%rcells_q(w%nwsq, 3))
    allocate(w%phwsdeg(w%nwsq))
    open(1, file = filename_phwscells, status = "old")
    open(2, file = filename_phwsdeg, status = "old")
    do iuc = 1,w%nwsq
       read(1, *) w%rcells_q(iuc, :)
       read(2, *) w%phwsdeg(iuc)
    end do
    close(1)
    close(2)

    allocate(w%rcells_g(w%nwsg, 3))
    allocate(w%gwsdeg(w%nwsg))
    open(1, file = filename_gwscells, status = "old")
    open(2, file = filename_gwsdeg, status = "old")
    do iuc = 1,w%nwsg
       read(1, *) w%rcells_g(iuc, :)
       read(2, *) w%gwsdeg(iuc)
    end do
    close(1)
    close(2)
  end subroutine read_EPW_Wannier

  subroutine el_wann_epw(w, nk, kvecs, energies, velocities, evecs)
    !! Wannier interpolate electrons on list of arb. k-vecs

    class(epw_wannier), intent(in) :: w
    integer(k4), intent(in) :: nk
    real(dp), intent(in) :: kvecs(nk,3) !Crystal coordinates
    real(dp), intent(out) :: energies(nk,w%numwannbands)
    real(dp), optional, intent(out) :: velocities(nk,w%numwannbands,3)
    complex(dp), optional, intent(out) :: evecs(nk,w%numwannbands,w%numwannbands)

    integer(k4) :: iuc, ib, jb, ipol, ik, nwork, tmp
    real(dp) :: rcart(3)
    real(dp),  allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)
    complex(dp) :: caux, H(w%numwannbands,w%numwannbands), &
         dH(3,w%numwannbands,w%numwannbands)

    !Catch error for optional velocity calculation
    if(present(velocities) .and. .not. present(evecs)) &
         call exit_with_message("In Wannier, velocity is present but not eigenvecs.")

    nwork = 1
    allocate(work(nwork))
    allocate(rwork(max(1,7*w%numwannbands)))

    do ik = 1,nk
       !Form Hamiltonian (H) and k-derivative of H (dH) 
       !from Hwann, rcells_k, and elwsdeg
       H = 0
       dH = 0
       do iuc = 1,w%nwsk
          caux = expi(twopi*dot_product(kvecs(ik,:),w%rcells_k(iuc,:)))&
               /w%elwsdeg(iuc)
          H = H + caux*w%Hwann(iuc,:,:)

          if(present(velocities)) then
             rcart = matmul(lattvecs,w%rcells_k(iuc,:))
             do ipol = 1,3
                dH(ipol,:,:) = dH(ipol,:,:) + &
                     oneI*rcart(ipol)*caux*w%Hwann(iuc,:,:)
             end do
          end if
       end do

       !Force Hermiticity
       do ib = 1, w%numwannbands
          do jb = ib + 1, w%numwannbands
             H(ib,jb) = (H(ib,jb) + conjg(H(jb,ib)))*0.5_dp
             H(jb,ib) = H(ib,jb)
          end do
       end do

       !Diagonalize H
       call zheev("V", "U", w%numwannbands, H(:,:), w%numwannbands, energies(ik,:), &
            work, -1, rwork, tmp)
       if(real(work(1)) > nwork) then
          nwork = nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V", "U", w%numwannbands, H(:,:), w%numwannbands, energies(ik,:), &
            work, nwork, rwork, tmp)

       if(present(evecs)) then
          evecs(ik,:,:)=transpose(H(:,:))
       end if

       if(present(velocities)) then
          !Calculate velocities using Feynman-Hellmann thm
          do ib = 1,w%numwannbands
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

  subroutine ph_wann_epw(w, nq, qvecs, energies, velocities, evecs)
    !! Wannier interpolate phonons on list of arb. q-vec

    class(epw_wannier), intent(in) :: w
    integer(k4), intent(in) :: nq
    real(dp), intent(in) :: qvecs(nq, 3) !Crystal coordinates
    real(dp), intent(out) :: energies(nq, w%numbranches)
    real(dp), intent(out), optional :: velocities(nq, w%numbranches, 3)
    complex(dp), intent(out), optional :: evecs(nq, w%numbranches, w%numbranches)

    integer(k4) :: iuc, ib, jb, ipol, iq, na, nb, nwork, aux
    complex(dp) :: caux
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)
    real(dp) :: omega2(w%numbranches), rcart(3), massnorm
    complex(dp) :: dynmat(w%numbranches, w%numbranches), ddynmat(3, w%numbranches, w%numbranches)

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
       do iuc = 1, w%nwsq
          caux = expi(twopi*dot_product(qvecs(iq, :), w%rcells_q(iuc, :)))&
               /w%phwsdeg(iuc)
          dynmat = dynmat + caux*w%Dphwann(iuc, :, :)

          if(present(velocities)) then
             rcart = matmul(lattvecs, w%rcells_q(iuc, :))
             do ipol=1, 3
                ddynmat(ipol, :, :) = ddynmat(ipol, :, :) + &
                     oneI*rcart(ipol)*caux*w%Dphwann(iuc, :, :)
             end do
          end if
       end do

       !Non-analytic correction
       if(polar) then
          call dyn_nonanalytic(w, matmul(reclattvecs,qvecs(iq,:))*bohr2nm, dynmat, ddynmat)
       end if

       !Force Hermiticity
       do ib = 1, w%numbranches
          do jb = ib + 1, w%numbranches
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
       call zheev("V", "U", w%numbranches, dynmat(:, :), w%numbranches, omega2, work, -1, rwork, aux)
       if(real(work(1)) > nwork) then
          nwork = nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V", "U", w%numbranches, dynmat(:, :), w%numbranches, omega2, work, nwork, rwork, aux)

       energies(iq, :) = sign(sqrt(abs(omega2)), omega2)
       if(present(evecs)) then
          evecs(iq, :, :) = transpose(dynmat(:, :))
       end if

       if(present(velocities)) then
          !Calculate velocities using Feynman-Hellmann thm
          do ib = 1, w%numbranches
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
       do ib = 1, w%numbranches
          if(energies(iq,ib) < -0.005_dp) then
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

  subroutine dyn_nonanalytic(w, q, dyn, ddyn)
    !! Calculate the long-range correction to the
    !! dynamical matrix and its derivative for a given phonon mode.
    !!
    !! q: the phonon wave vector in Cartesian coords., Bohr^-1
    !! (d)dyn: the (derivative of) dynamical matrix

    class(epw_wannier), intent(in) :: w
    real(dp), intent(in) :: q(3) !Cartesian
    complex(dp), intent(inout) :: dyn(w%numbranches,w%numbranches)
    complex(dp), intent(inout) :: ddyn(3,w%numbranches,w%numbranches)

    complex(dp) :: dyn_l(w%numbranches,w%numbranches)
    complex(dp) :: ddyn_l(3,w%numbranches,w%numbranches)
    real(dp) :: qeq,     &! <q+g| epsilon |q+g>
         arg, zag(3),zbg(3), g(3), gmax, alph, geg,&
         tpiba,dgeg(3),fnat(3),rr(numatoms,numatoms,3)
    integer(kind=4) :: na,nb,i,idim,jdim,ipol,jpol,m1,m2,m3,nq1,nq2,nq3
    complex(dp) :: fac, facqd, facq

    tpiba = twopi/twonorm(lattvecs(:,1))*bohr2nm

    !Recall that the phonon supercell in elphBolt is the
    !same as the EPW coarse phonon mesh.
    nq1 = qmesh(1)
    nq2 = qmesh(2)
    nq3 = qmesh(3)

    gmax= 14.d0 !dimensionless
    alph= tpiba**2 !bohr^-2
    geg = gmax*alph*4.0d0
    !In Ry units, qe = sqrt(2.0)
    fac = 8.d0*pi/(volume/bohr2nm**3)

    dyn_l = 0.d0
    ddyn_l = 0.d0
    do m1 = -nq1,nq1
       do m2 = -nq2,nq2
          do m3 = -nq3,nq3
             g(:) = (m1*reclattvecs(:,1)+m2*reclattvecs(:,2)+m3*reclattvecs(:,3))*bohr2nm
             qeq = dot_product(g,matmul(epsilon,g))

             if (qeq > 0.d0 .and. qeq/alph/4.d0 < gmax ) then
                facqd = exp(-qeq/alph/4.0d0)/qeq

                do na = 1,numatoms
                   zag(:)=matmul(g,born(:,:,na))
                   fnat(:)=0.d0
                   do nb = 1,numatoms
                      rr(na,nb,:) = (basis_cart(:,na)-basis_cart(:,nb))/bohr2nm
                      arg = dot_product(g,rr(na,nb,:))
                      zbg(:) = matmul(g,born(:,:,nb))
                      fnat(:) = fnat(:) + zbg(:)*expi(arg)
                   end do
                   do jpol=1,3
                      jdim=(na-1)*3+jpol
                      do ipol=1,3
                         idim=(na-1)*3+ipol
                         dyn_l(idim,jdim) = dyn_l(idim,jdim) - facqd * &
                              zag(ipol) * fnat(jpol)
                      end do
                   end do
                end do
             end if

             g = g + q
             qeq = dot_product(g,matmul(epsilon,g))
             if (qeq > 0.d0 .and. qeq/alph/4.d0 < gmax ) then
                facqd = exp(-qeq/alph/4.0d0)/qeq
                dgeg=matmul(epsilon+transpose(epsilon),g)
                do nb = 1,numatoms
                   zbg(:)=matmul(g,born(:,:,nb))
                   do na = 1,numatoms
                      rr(na,nb,:) = (basis_cart(:,na)-basis_cart(:,nb))/bohr2nm
                      zag(:)=matmul(g,born(:,:,na))
                      arg = dot_product(g,rr(na,nb,:))
                      facq = facqd*expi(arg)
                      do jpol=1,3
                         jdim=(nb-1)*3+jpol
                         do ipol=1,3
                            idim=(na-1)*3+ipol
                            dyn_l(idim,jdim) = dyn_l(idim,jdim) + facq * &
                                 zag(ipol)*zbg(jpol)
                            !Correction to derivative of dynmat
                            ddyn_l(:,idim,jdim)=ddyn_l(:,idim,jdim)+&
                                 facq*&
                                 ( zbg(jpol)*born(:,ipol,na)+zag(ipol)*born(:,jpol,nb)+&
                                 zag(ipol)*zbg(jpol)*(oneI*rr(na,nb,:)-&
                                 dgeg(:)/alph/4.0-dgeg(:)/qeq) )
                         end do
                      end do
                   end do
                end do
             end if
          end do
       end do
    end do
    dyn = dyn + dyn_l*fac
    ddyn = ddyn + ddyn_l*fac
  end subroutine dyn_nonanalytic

  subroutine long_range_prefac(w, q, uqs, glprefac)
    !! Calculate the long-range correction prefactor of
    !! the e-ph matrix element for a given phonon mode.
    !! q: phonon wvec in Cartesian coords., Bohr^-1
    !! uqs: phonon eigenfn for mode (s,q)
    !! glprefac: is the output in Ry units (EPW/QE)

    class(epw_wannier), intent(in) :: w
    real(dp), intent(in) :: q(3) !Cartesian
    complex(dp), intent(in) :: uqs(w%numbranches)
    complex(dp), intent(inout) :: glprefac

    real(dp) :: qeq,     &! <q+g| epsilon |q+g>
         arg, zaq, g(3), gmax, alph, geg,tpiba
    integer(k4) :: na,ipol, m1,m2,m3,nq1,nq2,nq3
    complex(dp) :: fac, facqd, facq

    tpiba = twopi/twonorm(lattvecs(:,1))*bohr2nm

    !Recall that the phonon supercell in elphBolt is the
    !same as the EPW coarse phonon mesh.
    nq1 = qmesh(1)
    nq2 = qmesh(2)
    nq3 = qmesh(3)

    gmax= 14.d0 !dimensionless
    alph= tpiba**2 !bohr^-2
    geg = gmax*alph*4.0d0
    !In Ry units, qe = sqrt(2.0)
    fac = 8.d0*pi/(volume/bohr2nm**3)*oneI
    glprefac = (0.d0,0.d0)

    do m1 = -nq1,nq1
       do m2 = -nq2,nq2
          do m3 = -nq3,nq3
             g(:) = (m1*reclattvecs(:,1)+m2*reclattvecs(:,2)+m3*reclattvecs(:,3))*bohr2nm + q
             qeq = dot_product(g,matmul(epsilon,g))

             if (qeq > 0.d0 .and. qeq/alph/4.d0 < gmax ) then
                facqd = exp(-qeq/alph/4.0d0)/qeq

                do na = 1,numatoms
                   arg = -dot_product(g,basis_cart(:,na))/bohr2nm
                   facq = facqd*expi(arg)
                   do ipol=1,3
                      zaq = dot_product(g,born(:,ipol,na))
                      glprefac = glprefac + facq*zaq*uqs(3*(na-1)+ipol)
                   end do
                end do
             end if
          end do
       end do
    end do
    glprefac = glprefac*fac
  end subroutine long_range_prefac

  function g2_epw(w, qvec, el_evec_k, el_evec_kp, ph_evec_q, ph_en, gmixed_ik)
    !! Function to calculate |g|^2.
    !! This works with EPW real space data
    !! qvec: phonon wave vector in crystal coords
    !! el_evec_k(kp): initial(final) electron eigenvector in bands m(n) 
    !! ph_evec_q: phonon eigenvector branch s 
    !! ph_en: phonon energy in mode (s,qvec)
    !! gmixed_ik: e-ph matrix element

    class(epw_wannier), intent(in) :: w
    real(dp),intent(in) :: qvec(3), ph_en
    complex(dp),intent(in) :: el_evec_k(w%numwannbands),&
         el_evec_kp(w%numwannbands), ph_evec_q(w%numbranches), &
         gmixed_ik(w%numwannbands, w%numwannbands, w%numbranches, w%nwsq)
    integer(kind=4) :: ip, ig, np, mp, sp, mtype
    complex(dp) :: caux, u(w%numbranches), UkpgkUk(w%numbranches, w%nwsq), &
         UkpgkUkuq(w%nwsq), gbloch, overlap(w%numwannbands,w%numwannbands), glprefac
    real(dp) :: g2_epw, unm

    !Mass normalize the phonon matrix
    do ip = 1, w%numbranches ! d.o.f of basis atoms
       !demux atom type from d.o.f
       mtype = (ip - 1)/3 + 1 
       !normalize and conjugate eigenvector
       u(ip) = ph_evec_q(ip)/sqrt(masses(atomtypes(mtype)))
    end do

    if(ph_en == 0) then !zero out matrix elements at Gamma
       g2_epw = 0
    else
       UkpgkUk = 0 !g(k,Rp) rotated by the electron U(k), U(k') matrices
       UkpgkUkuq = 0 !above quantity rotated by the phonon u(q) matrix
       gbloch = 0

       !Create the <n'|m'> overlap matrix
       do np = 1, w%numwannbands !over final electron band
          do mp = 1, w%numwannbands !over initial electron band
             overlap(mp,np) = conjg(el_evec_kp(np))*el_evec_k(mp)
          end do
       end do

       do ig = 1, w%nwsg !over matrix elements WS cell
          !Apply electron rotations
          do sp = 1, w%numbranches
             caux = 0
             do np = 1, w%numwannbands !over final electron band
                do mp = 1, w%numwannbands !over initial electron band
                   caux = caux + overlap(mp,np)*gmixed_ik(np, mp, sp, ig)
                end do
             end do
             UkpgkUk(sp, ig) = UkpgkUk(sp, ig) + caux
          end do
       end do

       do ig = 1, w%nwsg !over matrix elements WS cell
          !Apply phonon rotation
          UkpgkUkuq(ig) = UkpgkUkuq(ig) + dot_product(conjg(u),UkpgkUk(:, ig))
       end do

       do ig = 1, w%nwsg !over matrix elements WS cell
          !Fourier transform to q-space
          caux = exp(twopiI*dot_product(qvec, w%rcells_g(ig, :)))&
               /w%gwsdeg(ig)
          gbloch = gbloch + caux*UkpgkUkuq(ig)
       end do

       if(polar) then !Long-range correction
          unm = dot_product(conjg(el_evec_k),el_evec_kp)
          call long_range_prefac(w, matmul(reclattvecs,qvec)*bohr2nm,u,glprefac)
          gbloch = gbloch + glprefac*unm
       end if

       g2_epw = 0.5_dp*real(gbloch*conjg(gbloch))/ &
            ph_en*g2unitfactor !eV^2
    end if
  end function g2_epw

  subroutine calculate_g_bloch(w)
    !! MPI parallelizer of g2_epw over IBZ electron states within the Fermi window.
    !
    !In the FBZ and IBZ blocks a wave vector was retained when at least one
    !band belonged within the energy window. Here the bands outside energy window
    !will be skipped in the calculation as they are irrelevant for transport.
    !This subroutine will calculate the full Bloch rep. matrix elements for
    !all the energy window restricted electron-phonon processes for a given
    !irreducible initial electron state = (band, wave vector). 
    !This list will be written to disk in files tagged with the muxed state index.

    class(epw_wannier), intent(in) :: w
    integer(k4) :: nstates_irred, istate, m, ik, ik_muxed, n, ikp, s, &
         iq, start, end, chunk, ierr, k_indvec(3), kp_indvec(3), &
         q_indvec(3), count, g2size
    real(dp) :: k(3), kp(3), q(3)
    real(dp), allocatable :: g2_istate(:)
    complex(dp) :: gmixed_ik(w%numwannbands,w%numwannbands,w%numbranches,w%nwsq)
    character(len = 1024) :: filename

    call print_message("Doing g(k,Rp) -> |g(k,q)|^2 for all IBZ states...")

    !Length of g2_istate
    g2size = nstates_inwindow*w%numbranches
    allocate(g2_istate(g2size))

    !Total number of IBZ blocks states
    nstates_irred = nk_irred*w%numwannbands

    call distribute_points(nstates_irred, chunk, start, end)

    if(this_image() == 1) then
       print*, "   #states = ", nstates_irred
       print*, "   #states/process = ", chunk
    end if

    count = 0
    do istate = start, end !over IBZ blocks states
       !Initialize eligible process counter for this state
       count = 0

       !Demux state index into band (m) and wave vector (ik) indices
       call demux_state(istate,w%numwannbands,m,ik)

       !Get the muxed index of wave vector from the IBZ blocks index list
       ik_muxed = el_indexlist_irred(ik)

       !Apply energy window to initial (IBZ blocks) electron
       if(abs(el_ens_irred(ik, m) - enref) > fsthick) cycle

       !Load gmixed(ik) here for use inside the loops below
       call chdir(trim(adjustl(g2dir)))
       write (filename, '(I6)') ik
       filename = 'gmixed.ik'//trim(adjustl(filename))
       open(1,file=filename,status="old",access='stream')
       read(1) gmixed_ik
       close(1)
       call chdir(cwd)

       !Initial (IBZ blocks) wave vector (crystal coords.)
       k = el_wavevecs_irred(ik, :)

       !Convert from crystal to 0-based index vector
       k_indvec = nint(k*kmesh)

       !Run over final (FBZ blocks) electron wave vectors
       do ikp = 1, nk
          !Final wave vector (crystal coords.)
          kp = el_wavevecs(ikp, :)

          !Convert from crystal to 0-based index vector
          kp_indvec = nint(kp*kmesh)

          !Run over final electron bands
          do n = 1, w%numwannbands
             !Apply energy window to final electron
             if(abs(el_ens(ikp, n) - enref) > fsthick) cycle

             !Find interacting phonon wave vector
             !Note that q, k, and k' are all on the same mesh
             q_indvec = modulo(kp_indvec - k_indvec, kmesh) !0-based index vector
             q = q_indvec/dble(kmesh) !crystal coords.

             !Muxed index of q
             iq = mux_vector(q_indvec, kmesh, 0_dp)

             !Run over phonon branches
             do s = 1, w%numbranches
                !Increment g2 processes counter
                count = count + 1

                !Calculate |g_mns(<k>,q)|^2
                g2_istate(count) = g2_epw(w, q, el_evecs_irred(ik, m, :), &
                     el_evecs(ikp, n, :), ph_evecs(iq, s, :), &
                     ph_ens(iq, s), gmixed_ik)
             end do !s
          end do !n
       end do !ikp

       !Change to data output directory
       call chdir(trim(adjustl(g2dir)))

       !Write data in binary format
       !Note: this will overwrite existing data!
       write (filename, '(I6)') istate
       filename = 'g2.istate'//trim(adjustl(filename))
       open(1, file = trim(filename), status = 'replace', access = 'stream')
       write(1) g2_istate
       close(1)

       !Change back to working directory
       call chdir(cwd)
    end do

    !TODO at this point we can delete the gmixed disk data
    !call clear_gmixed_cache

    sync all
  end subroutine calculate_g_bloch

  subroutine gmixed_epw(w,ik)
    !! Calculate the Bloch-Wannier mixed rep. e-ph matrix elements g(k,Rp),
    !! where k is an IBZ electron wave vector and Rp is a phonon unit cell.
    !! Note: this step *DOES NOT* perform the rotation over the Wannier bands space.
    !
    !The result will be saved to disk tagged with k-index.

    class(epw_wannier), intent(in) :: w
    integer(k4), intent(in) :: ik
    integer(k4) :: iuc
    complex(dp) :: caux
    complex(dp), allocatable:: gmixed(:,:,:,:)
    real(dp) :: kvec(3)
    character(len = 1024) :: filename

    allocate(gmixed(w%numwannbands, w%numwannbands, w%numbranches, w%nwsq))

    !Electron wave vector (crystal coords.) in IBZ blocks 
    kvec = el_wavevecs_irred(ik, :)

    !Fourier transform to k-space
    gmixed = 0
    do iuc = 1,w%nwsk
       caux = expi(twopi*dot_product(kvec, w%rcells_k(iuc,:)))/w%elwsdeg(iuc)
       gmixed(:,:,:,:) = gmixed(:,:,:,:) + caux*w%gwann(:,:,iuc,:,:)
    end do

    !Change to data output directory
    call chdir(trim(adjustl(g2dir)))

    !Write data in binary format
    !Note: this will overwrite existing data!
    write (filename, '(I6)') ik
    filename = 'gmixed.ik'//trim(adjustl(filename))
    open(1, file = trim(filename), status = 'replace', access = 'stream')
    write(1) gmixed
    close(1)

    !Change back to working directory
    call chdir(cwd)
  end subroutine gmixed_epw

  subroutine calculate_g_mixed(w)
    !! Parallel driver of gmixed_epw over IBZ electron wave vectors

    class(epw_wannier), intent(in) :: w
    integer(k4) :: ik, ikstart, ikend, chunk, ierr

    call print_message("Doing g(Re,Rp) -> g(k,Rp) for all IBZ k...")

    call distribute_points(nk_irred, chunk, ikstart, ikend)

    if(this_image() == 1) then
       print*, "   #k = ", nk_irred
       print*, "   #k/process = ", chunk
    end if

    do ik = ikstart, ikend
       call gmixed_epw(w, ik)
    end do

    sync all
  end subroutine calculate_g_mixed

  !For testing and debugging:
  subroutine gmixed_epw_gamma(w)
    !! Calculate the Bloch-Wannier mixed rep. e-ph matrix elements g(k,Rp),
    !! where k is an IBZ electron wave vector and Rp is a phonon unit cell.
    !! This is for the electron at gamma.

    class(epw_wannier), intent(in) :: w
    integer(k4) :: iuc
    complex(dp) :: caux
    complex(dp), allocatable:: gmixed(:, :, :, :)
    character(len = 1024) :: filename

    allocate(gmixed(w%numwannbands, w%numwannbands, w%numbranches, w%nwsq))

    !Fourier transform to k-space
    gmixed = 0
    do iuc = 1, w%nwsk
       caux = (1.0_dp,0.0_dp)/w%elwsdeg(iuc)           
       gmixed(:, :, :, :) = gmixed(:, :, :, :) + &
            caux*w%gwann(:, :, iuc, :, :)
    end do

    !Change to data output directory
    call chdir(trim(adjustl(g2dir)))

    !Write data in binary format
    !Note: this will overwrite existing data!
    filename = 'gmixed.gamma' !//trim(adjustl(filename))
    open(1, file = trim(filename), status = 'replace', access = 'stream')
    write(1) gmixed
    close(1)

    !Change back to working directory
    call chdir(cwd)
  end subroutine gmixed_epw_gamma

  !Unit tester for Wannier interpolation with EPW inputs.
  subroutine test_wannier(w)
    use params
    use data

    class(epw_wannier), intent(in) :: w
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
    allocate(ph_ens_path(nqpath,w%numbranches), ph_vels_path(nqpath,w%numbranches,3),&
         ph_evecs_path(nqpath,w%numbranches,w%numbranches))
    call ph_wann_epw(w, nqpath, qpathvecs, ph_ens_path, ph_vels_path, ph_evecs_path)

    !Output phonon dispersions
    write(saux,"(I0)") w%numbranches
    open(1,file="ph_ens_path",status="replace")
    do i = 1, nqpath
       write(1,"("//trim(adjustl(saux))//"E20.10)") ph_ens_path(i,:)
    end do
    close(1)

    !Test electron Wannier interpolation
    print*, 'Doing electron Wannier test...'

    !Calculate electron bands
    allocate(el_ens_path(nqpath,w%numwannbands))
    call el_wann_epw(w, nqpath, qpathvecs, el_ens_path)

    !Output electron dispersions
    write(saux,"(I0)") w%numwannbands
    open(1,file="el_ens_path",status="replace")
    do i=1,nqpath
       write(1,"("//trim(adjustl(saux))//"E20.10)") el_ens_path(i,:)
    end do
    close(1)

    !Test matrix elements Wannier interpolations
    print*, 'Doing g(k,q) Wannier test...'

    !Initial electron at Gamma
    k = (/0.0_dp, 0.0_dp, 0.0_dp/)
    allocate(el_ens_k(w%numwannbands), el_vels_k(w%numwannbands,3),&
         el_evecs_k(w%numwannbands,w%numwannbands))
    call el_wann_epw(w, 1_k4, k, el_ens_k, el_vels_k, el_evecs_k)

    !Calculate g(k, Rp)
    call gmixed_epw_gamma(w)

    !Load gmixed from file
    !Change to data output directory
    call chdir(trim(adjustl(g2dir)))
    allocate(gmixed_gamma(w%numwannbands, w%numwannbands, w%numbranches, w%nwsq))
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
    allocate(el_ens_kp(w%numwannbands), el_vels_kp(w%numwannbands,3),&
         el_evecs_kp(w%numwannbands,w%numwannbands))
    allocate(g2_qpath(nqpath))
    do i = 1, nqpath
       kp = qpathvecs(i, :) !for k = Gamma 

       !Calculate electrons at path point
       call el_wann_epw(w, 1_k4, kp, el_ens_kp, el_vels_kp, el_evecs_kp)

       !Calculate |g(k,k')|^2
       !(q is the same as k')
       g2_qpath(i) = g2_epw(w, kp, el_evecs_k(m, :), &
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

end module wannier
