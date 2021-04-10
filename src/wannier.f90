module wannier_module
  !! Module containing type and procedures related to Wannierization.

  use params, only: dp, k4, Ryd2eV, Ryd2radTHz, oneI, pi, twopi, twopiI, &
       Ryd2amu, bohr2nm
  use misc, only: exit_with_message, print_message, expi, twonorm, &
       distribute_points, demux_state, mux_vector
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  
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
     integer(k4) :: coarse_qmesh(3)
     !! Coarse phonon wave vector mesh in Wannier calculation.
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

     procedure :: read=>read_EPW_Wannier, el_wann_epw, ph_wann_epw, &
          gmixed_epw, g2_epw
     procedure :: test_wannier

  end type epw_wannier

contains

  subroutine read_EPW_Wannier(wann)
    !! Read Wannier representation of the hamiltonian, dynamical matrix, and the
    !! e-ph matrix elements from file epwdata.fmt.

    class(epw_wannier), intent(out) :: wann

    !Local variables
    integer(k4) :: iuc, ib, jb
    integer(k4) :: coarse_qmesh(3)

    namelist /wannier/ coarse_qmesh

    !Open input file
    open(1, file = 'input.nml', status = 'old')

    wann%coarse_qmesh = (/1, 1, 1/)
    read(1, nml = wannier)
    if(any(coarse_qmesh <= 0)) then
       call exit_with_message('Bad input(s) in wannier.')
    end if
    wann%coarse_qmesh = coarse_qmesh
    
    !Close input file
    close(1)
    
    open(1,file=filename_epwdata,status='old')
    read(1,*) !ef
    read(1,*) wann%numwannbands, wann%nwsk, wann%numbranches, wann%nwsq, wann%nwsg
    read(1,*) !zstar, epsil: non-zero only for polar materials

    !Read real space hamiltonian
    call print_message("Reading Wannier rep. Hamiltonian...")
    allocate(wann%Hwann(wann%nwsk,wann%numwannbands,wann%numwannbands))
    do ib = 1,wann%numwannbands
       do jb = 1,wann%numwannbands
          do iuc = 1,wann%nwsk !Number of real space electron cells
             read (1, *) wann%Hwann(iuc,ib,jb)
          end do
       end do
    end do

    !Read real space dynamical matrix
    call print_message("Reading Wannier rep. dynamical matrix...")
    allocate(wann%Dphwann(wann%nwsq,wann%numbranches,wann%numbranches))
    do ib = 1,wann%numbranches
       do jb = 1,wann%numbranches
          do iuc = 1,wann%nwsq !Number of real space phonon cells
             read (1, *) wann%Dphwann(iuc,ib,jb)
          end do
       end do
    end do
    close(1)

    !Read real space matrix elements
    call print_message("Reading Wannier rep. e-ph vertex...")
    open(1, file = filename_epwgwann, status = 'old', access = 'stream')
    allocate(wann%gwann(wann%numwannbands,wann%numwannbands,wann%nwsk,wann%numbranches,wann%nwsg))
    wann%gwann = 0.0_dp
    read(1) wann%gwann
    close(1)

    !Read cell maps of q, k, g meshes.
    call print_message("Reading Wannier cells and multiplicities...")
    allocate(wann%rcells_k(wann%nwsk,3))
    allocate(wann%elwsdeg(wann%nwsk))
    open(1, file = filename_elwscells, status = "old")
    open(2, file = filename_elwsdeg, status = "old")
    do iuc = 1,wann%nwsk
       read(1, *) wann%rcells_k(iuc, :)
       read(2, *) wann%elwsdeg(iuc)
    end do
    close(1)
    close(2)

    allocate(wann%rcells_q(wann%nwsq, 3))
    allocate(wann%phwsdeg(wann%nwsq))
    open(1, file = filename_phwscells, status = "old")
    open(2, file = filename_phwsdeg, status = "old")
    do iuc = 1,wann%nwsq
       read(1, *) wann%rcells_q(iuc, :)
       read(2, *) wann%phwsdeg(iuc)
    end do
    close(1)
    close(2)

    allocate(wann%rcells_g(wann%nwsg, 3))
    allocate(wann%gwsdeg(wann%nwsg))
    open(1, file = filename_gwscells, status = "old")
    open(2, file = filename_gwsdeg, status = "old")
    do iuc = 1,wann%nwsg
       read(1, *) wann%rcells_g(iuc, :)
       read(2, *) wann%gwsdeg(iuc)
    end do
    close(1)
    close(2)
  end subroutine read_EPW_Wannier

  subroutine el_wann_epw(wann, crys, nk, kvecs, energies, velocities, evecs)
    !! Wannier interpolate electrons on list of arb. k-vecs

    class(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    integer(k4), intent(in) :: nk
    real(dp), intent(in) :: kvecs(nk,3) !Crystal coordinates
    real(dp), intent(out) :: energies(nk,wann%numwannbands)
    real(dp), optional, intent(out) :: velocities(nk,wann%numwannbands,3)
    complex(dp), optional, intent(out) :: evecs(nk,wann%numwannbands,wann%numwannbands)

    !Local variables
    integer(k4) :: iuc, ib, jb, ipol, ik, nwork, tmp
    real(dp) :: rcart(3)
    real(dp),  allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)
    complex(dp) :: caux, H(wann%numwannbands,wann%numwannbands), &
         dH(3,wann%numwannbands,wann%numwannbands)
    
    !Catch error for optional velocity calculation
    if(present(velocities) .and. .not. present(evecs)) &
         call exit_with_message("In Wannier, velocity is present but not eigenvecs.")

    nwork = 1
    allocate(work(nwork))
    allocate(rwork(max(1,7*wann%numwannbands)))
    
    do ik = 1,nk
       !Form Hamiltonian (H) and k-derivative of H (dH) 
       !from Hwann, rcells_k, and elwsdeg
       H = 0
       dH = 0
       do iuc = 1,wann%nwsk
          caux = expi(twopi*dot_product(kvecs(ik,:),wann%rcells_k(iuc,:)))&
               /wann%elwsdeg(iuc)
          H = H + caux*wann%Hwann(iuc,:,:)

          if(present(velocities)) then
             rcart = matmul(crys%lattvecs,wann%rcells_k(iuc,:))
             do ipol = 1,3
                dH(ipol,:,:) = dH(ipol,:,:) + &
                     oneI*rcart(ipol)*caux*wann%Hwann(iuc,:,:)
             end do
          end if
       end do

       !Force Hermiticity
       do ib = 1, wann%numwannbands
          do jb = ib + 1, wann%numwannbands
             H(ib,jb) = (H(ib,jb) + conjg(H(jb,ib)))*0.5_dp
             H(jb,ib) = H(ib,jb)
          end do
       end do

       !Diagonalize H
       call zheev("V", "U", wann%numwannbands, H(:,:), wann%numwannbands, energies(ik,:), &
            work, -1, rwork, tmp)
       if(real(work(1)) > nwork) then
          nwork = nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V", "U", wann%numwannbands, H(:,:), wann%numwannbands, energies(ik,:), &
            work, nwork, rwork, tmp)

       if(present(evecs)) then
          evecs(ik,:,:)=transpose(H(:,:))
       end if

       if(present(velocities)) then
          !Calculate velocities using Feynman-Hellmann thm
          do ib = 1,wann%numwannbands
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

  subroutine ph_wann_epw(wann, crys, nq, qvecs, energies, velocities, evecs)
    !! Wannier interpolate phonons on list of arb. q-vec

    class(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys

    !Local variables
    integer(k4), intent(in) :: nq
    real(dp), intent(in) :: qvecs(nq, 3) !Crystal coordinates
    real(dp), intent(out) :: energies(nq, wann%numbranches)
    real(dp), intent(out), optional :: velocities(nq, wann%numbranches, 3)
    complex(dp), intent(out), optional :: evecs(nq, wann%numbranches, wann%numbranches)

    integer(k4) :: iuc, ib, jb, ipol, iq, na, nb, nwork, aux
    complex(dp) :: caux
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)
    real(dp) :: omega2(wann%numbranches), rcart(3), massnorm
    complex(dp) :: dynmat(wann%numbranches, wann%numbranches), &
         ddynmat(3, wann%numbranches, wann%numbranches)

    !Catch error for optional velocity calculation
    if(present(velocities) .and. .not. present(evecs)) &
         call exit_with_message("In Wannier, velocity is present but not eigenvecs.")
    
    nwork = 1
    allocate(work(nwork))
    allocate(rwork(max(1, 9*crys%numatoms-2)))

    do iq = 1, nq
       !Form dynamical matrix (dynmat) and q-derivative of dynmat (ddynmat) 
       !from Dphwann, rcells_q, and phwsdeg
       dynmat = 0
       ddynmat = 0
       do iuc = 1, wann%nwsq
          caux = expi(twopi*dot_product(qvecs(iq, :), wann%rcells_q(iuc, :)))&
               /wann%phwsdeg(iuc)
          dynmat = dynmat + caux*wann%Dphwann(iuc, :, :)

          if(present(velocities)) then
             rcart = matmul(crys%lattvecs, wann%rcells_q(iuc, :))
             do ipol=1, 3
                ddynmat(ipol, :, :) = ddynmat(ipol, :, :) + &
                     oneI*rcart(ipol)*caux*wann%Dphwann(iuc, :, :)
             end do
          end if
       end do

       !Non-analytic correction
       if(crys%polar) then
          call dyn_nonanalytic(wann, crys, matmul(crys%reclattvecs,qvecs(iq,:))*bohr2nm, dynmat, ddynmat)
       end if

       !Force Hermiticity
       do ib = 1, wann%numbranches
          do jb = ib + 1, wann%numbranches
             dynmat(ib, jb) = (dynmat(ib, jb) + conjg(dynmat(jb, ib)))*0.5_dp
             dynmat(jb, ib) = dynmat(ib, jb)
          end do
       end do

       !Mass normalize
       do na = 1, crys%numatoms
          do nb = 1, crys%numatoms
             massnorm = 1.d0/sqrt(crys%masses(crys%atomtypes(na))*&
                  crys%masses(crys%atomtypes(nb)))*Ryd2amu
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
       call zheev("V", "U", wann%numbranches, dynmat(:, :), wann%numbranches, omega2, work, -1, rwork, aux)
       if(real(work(1)) > nwork) then
          nwork = nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V", "U", wann%numbranches, dynmat(:, :), wann%numbranches, omega2, work, nwork, rwork, aux)

       energies(iq, :) = sign(sqrt(abs(omega2)), omega2)
       if(present(evecs)) then
          evecs(iq, :, :) = transpose(dynmat(:, :))
       end if

       if(present(velocities)) then
          !Calculate velocities using Feynman-Hellmann thm
          do ib = 1, wann%numbranches
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
       if(all(qvecs(iq,:) == 0)) then
          energies(iq, 1:3) = 0
          if(present(velocities)) then
             velocities(iq, :, :) = 0
          end if
       end if

       !Handle negative energy phonons
       do ib = 1, wann%numbranches
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

  subroutine dyn_nonanalytic(wann, crys, q, dyn, ddyn)
    !! Calculate the long-range correction to the
    !! dynamical matrix and its derivative for a given phonon mode.
    !!
    !! q: the phonon wave vector in Cartesian coords., Bohr^-1
    !! (d)dyn: the (derivative of) dynamical matrix

    class(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    
    !Local variables
    real(dp), intent(in) :: q(3) !Cartesian
    complex(dp), intent(inout) :: dyn(wann%numbranches,wann%numbranches)
    complex(dp), intent(inout) :: ddyn(3,wann%numbranches,wann%numbranches)

    complex(dp) :: dyn_l(wann%numbranches,wann%numbranches)
    complex(dp) :: ddyn_l(3,wann%numbranches,wann%numbranches)
    real(dp) :: qeq,     &! <q+g| epsilon |q+g>
         arg, zag(3), zbg(3), g(3), gmax, alph, geg,&
         tpiba, dgeg(3), fnat(3), rr(crys%numatoms,crys%numatoms,3)
    integer(k4) :: na,nb,i,idim,jdim,ipol,jpol,m1,m2,m3,nq1,nq2,nq3
    complex(dp) :: fac, facqd, facq

    tpiba = twopi/twonorm(crys%lattvecs(:,1))*bohr2nm

    !Recall that the phonon supercell in elphbolt is the
    !same as the EPW coarse phonon mesh.
    nq1 = wann%coarse_qmesh(1)
    nq2 = wann%coarse_qmesh(2)
    nq3 = wann%coarse_qmesh(3)

    gmax= 14.0_dp !dimensionless
    alph= tpiba**2 !bohr^-2
    geg = gmax*alph*4.0d0
    !In Ry units, qe = sqrt(2.0)
    fac = 8.0_dp*pi/(crys%volume/bohr2nm**3)

    dyn_l = 0.0_dp
    ddyn_l = 0.0_dp
    do m1 = -nq1,nq1
       do m2 = -nq2,nq2
          do m3 = -nq3,nq3
             g(:) = (m1*crys%reclattvecs(:,1)+m2*crys%reclattvecs(:,2)+m3*crys%reclattvecs(:,3))*bohr2nm
             qeq = dot_product(g,matmul(crys%epsilon,g))

             if (qeq > 0.0_dp .and. qeq/alph/4.0_dp < gmax ) then
                facqd = exp(-qeq/alph/4.0_dp)/qeq

                do na = 1,crys%numatoms
                   zag(:)=matmul(g,crys%born(:,:,na))
                   fnat(:)=0.0_dp
                   do nb = 1,crys%numatoms
                      rr(na,nb,:) = (crys%basis_cart(:,na)-crys%basis_cart(:,nb))/bohr2nm
                      arg = dot_product(g,rr(na,nb,:))
                      zbg(:) = matmul(g,crys%born(:,:,nb))
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
             qeq = dot_product(g,matmul(crys%epsilon,g))
             if (qeq > 0.0_dp .and. qeq/alph/4.0_dp < gmax ) then
                facqd = exp(-qeq/alph/4.0d0)/qeq
                dgeg=matmul(crys%epsilon+transpose(crys%epsilon),g)
                do nb = 1,crys%numatoms
                   zbg(:)=matmul(g,crys%born(:,:,nb))
                   do na = 1,crys%numatoms
                      rr(na,nb,:) = (crys%basis_cart(:,na)-crys%basis_cart(:,nb))/bohr2nm
                      zag(:)=matmul(g,crys%born(:,:,na))
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
                                 ( zbg(jpol)*crys%born(:,ipol,na)+zag(ipol)*crys%born(:,jpol,nb)+&
                                 zag(ipol)*zbg(jpol)*(oneI*rr(na,nb,:)-&
                                 dgeg(:)/alph/4.0_dp-dgeg(:)/qeq) )
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

  function g2_epw(wann, crys, qvec, el_evec_k, el_evec_kp, ph_evec_q, ph_en, gmixed_ik)
    !! Function to calculate |g|^2.
    !! This works with EPW real space data
    !! qvec: phonon wave vector in crystal coords
    !! el_evec_k(kp): initial(final) electron eigenvector in bands m(n) 
    !! ph_evec_q: phonon eigenvector branchs 
    !! ph_en: phonon energy in mode (s,qvec)
    !! gmixed_ik: e-ph matrix element

    class(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    
    real(dp),intent(in) :: qvec(3), ph_en
    complex(dp),intent(in) :: el_evec_k(wann%numwannbands),&
         el_evec_kp(wann%numwannbands), ph_evec_q(wann%numbranches), &
         gmixed_ik(wann%numwannbands, wann%numwannbands, wann%numbranches, wann%nwsq)
    integer(kind=4) :: ip, ig, np, mp, sp, mtype
    complex(dp) :: caux, u(wann%numbranches), UkpgkUk(wann%numbranches, wann%nwsq), &
         UkpgkUkuq(wann%nwsq), gbloch, overlap(wann%numwannbands,wann%numwannbands), glprefac
    real(dp) :: g2_epw, unm

    !Mass normalize the phonon matrix
    do ip = 1, wann%numbranches ! d.o.f of basis atoms
       !demux atom type from d.o.f
       mtype = (ip - 1)/3 + 1 
       !normalize and conjugate eigenvector
       u(ip) = ph_evec_q(ip)/sqrt(crys%masses(crys%atomtypes(mtype)))
    end do

    if(ph_en == 0) then !zero out matrix elements at Gamma
       g2_epw = 0
    else
       UkpgkUk = 0 !g(k,Rp) rotated by the electron U(k), U(k') matrices
       UkpgkUkuq = 0 !above quantity rotated by the phonon u(q) matrix
       gbloch = 0

       !Create the <n'|m'> overlap matrix
       do np = 1, wann%numwannbands !over final electron band
          do mp = 1, wann%numwannbands !over initial electron band
             overlap(mp,np) = conjg(el_evec_kp(np))*el_evec_k(mp)
          end do
       end do

       do ig = 1, wann%nwsg !over matrix elements WS cell
          !Apply electron rotations
          do sp = 1, wann%numbranches
             caux = 0
             do np = 1, wann%numwannbands !over final electron band
                do mp = 1, wann%numwannbands !over initial electron band
                   caux = caux + overlap(mp,np)*gmixed_ik(np, mp, sp, ig)
                end do
             end do
             UkpgkUk(sp, ig) = UkpgkUk(sp, ig) + caux
          end do
       end do

       do ig = 1, wann%nwsg !over matrix elements WS cell
          !Apply phonon rotation
          UkpgkUkuq(ig) = UkpgkUkuq(ig) + dot_product(conjg(u),UkpgkUk(:, ig))
       end do

       do ig = 1, wann%nwsg !over matrix elements WS cell
          !Fourier transform to q-space
          caux = exp(twopiI*dot_product(qvec, wann%rcells_g(ig, :)))&
               /wann%gwsdeg(ig)
          gbloch = gbloch + caux*UkpgkUkuq(ig)
       end do

       if(crys%polar) then !Long-range correction
          unm = dot_product(conjg(el_evec_k),el_evec_kp)
          call long_range_prefac(wann, crys, &
               matmul(crys%reclattvecs,qvec)*bohr2nm,u,glprefac)
          gbloch = gbloch + glprefac*unm
       end if

       g2_epw = 0.5_dp*real(gbloch*conjg(gbloch))/ &
            ph_en*g2unitfactor !eV^2
    end if
  end function g2_epw
  
  subroutine long_range_prefac(wann, crys, q, uqs, glprefac)
    !! Calculate the long-range correction prefactor of
    !! the e-ph matrix element for a given phonon mode.
    !! q: phonon wvec in Cartesian coords., Bohr^-1
    !! uqs: phonon eigenfn for mode (s,q)
    !! glprefac: is the output in Ry units (EPW/QE)

    class(epw_wannier), intent(in) :: wann
    type(crystal), intent(in) :: crys
    
    real(dp), intent(in) :: q(3) !Cartesian
    complex(dp), intent(in) :: uqs(wann%numbranches)
    complex(dp), intent(inout) :: glprefac

    real(dp) :: qeq,     &! <q+g| epsilon |q+g>
         arg, zaq, g(3), gmax, alph, geg,tpiba
    integer(k4) :: na,ipol, m1,m2,m3,nq1,nq2,nq3
    complex(dp) :: fac, facqd, facq

    tpiba = twopi/twonorm(crys%lattvecs(:,1))*bohr2nm

    !Recall that the phonon supercell in elphbolt is the
    !same as the EPW coarse phonon mesh.
    nq1 = wann%coarse_qmesh(1)
    nq2 = wann%coarse_qmesh(2)
    nq3 = wann%coarse_qmesh(3)

    gmax= 14.d0 !dimensionless
    alph= tpiba**2 !bohr^-2
    geg = gmax*alph*4.0d0
    !In Ry units, qe = sqrt(2.0)
    fac = 8.d0*pi/(crys%volume/bohr2nm**3)*oneI
    glprefac = (0.d0,0.d0)

    do m1 = -nq1,nq1
       do m2 = -nq2,nq2
          do m3 = -nq3,nq3
             g(:) = (m1*crys%reclattvecs(:,1)+m2*crys%reclattvecs(:,2)+m3*crys%reclattvecs(:,3))*bohr2nm + q
             qeq = dot_product(g,matmul(crys%epsilon,g))

             if (qeq > 0.d0 .and. qeq/alph/4.d0 < gmax ) then
                facqd = exp(-qeq/alph/4.0d0)/qeq

                do na = 1,crys%numatoms
                   arg = -dot_product(g,crys%basis_cart(:,na))/bohr2nm
                   facq = facqd*expi(arg)
                   do ipol=1,3
                      zaq = dot_product(g,crys%born(:,ipol,na))
                      glprefac = glprefac + facq*zaq*uqs(3*(na-1)+ipol)
                   end do
                end do
             end if
          end do
       end do
    end do
    glprefac = glprefac*fac
  end subroutine long_range_prefac

  subroutine gmixed_epw(wann, num, ik, kvec)
    !! Calculate the Bloch-Wannier mixed rep. e-ph matrix elements g(k,Rp),
    !! where kvec is an IBZ electron wave vector and Rp is a phonon unit cell.
    !! Note: this step *DOES NOT* perform the rotation over the Wannier bands space.
    !!
    !! The result will be saved to disk tagged with k-index.

    class(epw_wannier), intent(in) :: wann
    type(numerics), intent(in) :: num
    integer(k4), intent(in) :: ik
    real(dp), intent(in) :: kvec(3)

    !Local variables
    integer(k4) :: iuc
    complex(dp) :: caux
    complex(dp), allocatable:: gmixed(:,:,:,:)

    character(len = 1024) :: filename

    allocate(gmixed(wann%numwannbands, wann%numwannbands, wann%numbranches, wann%nwsq))

    !Fourier transform to k-space
    gmixed = 0
    do iuc = 1,wann%nwsk
       caux = expi(twopi*dot_product(kvec, wann%rcells_k(iuc,:)))/wann%elwsdeg(iuc)
       gmixed(:,:,:,:) = gmixed(:,:,:,:) + caux*wann%gwann(:,:,iuc,:,:)
    end do

    !Change to data output directory
    call chdir(trim(adjustl(num%g2dir)))

    !Write data in binary format
    !Note: this will overwrite existing data!
    write (filename, '(I9)') ik
    filename = 'gmixed.ik'//trim(adjustl(filename))
    open(1, file = trim(filename), status = 'replace', access = 'stream')
    write(1) gmixed
    close(1)

    !Change back to working directory
    call chdir(num%cwd)
  end subroutine gmixed_epw

  !For testing and debugging:
  subroutine gmixed_epw_gamma(wann, num)
    !! Calculate the Bloch-Wannier mixed rep. e-ph matrix elements g(k,Rp),
    !! where k is an IBZ electron wave vector and Rp is a phonon unit cell.
    !! This is for the electron at gamma.

    class(epw_wannier), intent(in) :: wann
    type(numerics), intent(in) :: num
    
    integer(k4) :: iuc
    complex(dp) :: caux
    complex(dp), allocatable:: gmixed(:, :, :, :)
    character(len = 1024) :: filename

    allocate(gmixed(wann%numwannbands, wann%numwannbands, wann%numbranches, wann%nwsq))

    !Fourier transform to k-space
    gmixed = 0
    do iuc = 1, wann%nwsk
       caux = (1.0_dp,0.0_dp)/wann%elwsdeg(iuc)           
       gmixed(:, :, :, :) = gmixed(:, :, :, :) + &
            caux*wann%gwann(:, :, iuc, :, :)
    end do

    !Change to data output directory
    call chdir(trim(adjustl(num%g2dir)))

    !Write data in binary format
    !Note: this will overwrite existing data!
    filename = 'gmixed.gamma'
    open(1, file = trim(filename), status = 'replace', access = 'stream')
    write(1) gmixed
    close(1)

    !Change back to working directory
    call chdir(num%cwd)
  end subroutine gmixed_epw_gamma
  
  subroutine test_wannier(wann, crys, num)
    !! Unit tester for Wannier interpolation with EPW inputs.
    
    class(epw_wannier), intent(in) :: wann
    type(numerics), intent(in) :: num
    type(crystal), intent(in) :: crys

    !Local variables
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
    !nqpath = 601 !251 !601
    nqpath = 251
    allocate(qpathvecs(nqpath,3))
    open(1,file=trim('highsymqpts.txt'),status='old')
    read(1,*) !skip first line
    do i = 1, nqpath
       read(1,*) qpathvecs(i,:)
    end do

    !Test phonon Wannier interpolation
    print*, 'Doing phonon Wannier test...'

    !Calculate phonon dispersions
    allocate(ph_ens_path(nqpath,wann%numbranches), ph_vels_path(nqpath,wann%numbranches,3),&
         ph_evecs_path(nqpath,wann%numbranches,wann%numbranches))
    call ph_wann_epw(wann, crys, nqpath, qpathvecs, ph_ens_path, ph_vels_path, ph_evecs_path)

    !Output phonon dispersions
    write(saux,"(I0)") wann%numbranches
    open(1,file="ph_ens_path",status="replace")
    do i = 1, nqpath
       write(1,"("//trim(adjustl(saux))//"E20.10)") ph_ens_path(i,:)
    end do
    close(1)

    !Test electron Wannier interpolation
    print*, 'Doing electron Wannier test...'

    !Calculate electron bands
    allocate(el_ens_path(nqpath,wann%numwannbands))
    call el_wann_epw(wann, crys, nqpath, qpathvecs, el_ens_path)

    !Output electron dispersions
    write(saux,"(I0)") wann%numwannbands
    open(1,file="el_ens_path",status="replace")
    do i=1,nqpath
       write(1,"("//trim(adjustl(saux))//"E20.10)") el_ens_path(i,:)
    end do
    close(1)

    !Test matrix elements Wannier interpolations
    print*, 'Doing g(k,q) Wannier test...'

    !Initial electron at Gamma
    k = (/0.0_dp, 0.0_dp, 0.0_dp/)
    allocate(el_ens_k(wann%numwannbands), el_vels_k(wann%numwannbands,3),&
         el_evecs_k(wann%numwannbands,wann%numwannbands))
    call el_wann_epw(wann, crys, 1_k4, k, el_ens_k, el_vels_k, el_evecs_k)

    !Calculate g(k, Rp)
    call gmixed_epw_gamma(wann, num)

    !Load gmixed from file
    !Change to data output directory
    call chdir(trim(adjustl(num%g2dir)))
    allocate(gmixed_gamma(wann%numwannbands, wann%numwannbands, wann%numbranches, wann%nwsq))
    filename = 'gmixed.gamma' !//trim(adjustl(filename))
    open(1,file=filename,status="old",access='stream')
    read(1) gmixed_gamma
    close(1)
    !Change back to working directory
    call chdir(num%cwd)

    !print*, 'Rydberg2eV = ', Ryd2eV
    !print*, 'Rydberg2amu = ', Ryd2amu

    !Calculate g(k, k') where k' = (k + q) modulo G
    !q is along the path loaded from the file
    !phonon branch, s = 3 (LA)
    !initial (final) electron band, m(n) = 1(1) 
    m = 1
    n = 1
    s = 6
    allocate(el_ens_kp(wann%numwannbands), el_vels_kp(wann%numwannbands,3),&
         el_evecs_kp(wann%numwannbands,wann%numwannbands))
    allocate(g2_qpath(nqpath))
    do i = 1, nqpath
       kp = qpathvecs(i, :) !for k = Gamma 

       !Calculate electrons at path point
       call el_wann_epw(wann, crys, 1_k4, kp, el_ens_kp, el_vels_kp, el_evecs_kp)

       !Calculate |g(k,k')|^2
       !(q is the same as k')
       g2_qpath(i) = g2_epw(wann, crys, kp, el_evecs_k(m, :), &
            el_evecs_kp(n, :), ph_evecs_path(i, s, :), ph_ens_path(i, s), &
            gmixed_gamma)

       !print*, el_ens_kp(n), ph_ens_path(i,s)*1e3, sqrt(g2_qpath(i))*1.0d3
    end do

!!$    print*, 'done calculating g'
    !open(1,file='g_qpath_123',status="replace")
    open(1,file='g_qpath_116',status="replace")
    write(1,*) '#|g_SE| [eV]'
    do i=1,nqpath
       write(1,"(E20.10)") sqrt(g2_qpath(i))
    end do
    close(1)
  end subroutine test_wannier
end module wannier_module
