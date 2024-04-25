! Copyright 2020-2024 elphbolt contributors.
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

module nano_module
  !! Module containing type and procedures related to nanostructures.
  !! The effective quantities are provided through suppression factors
  !! which come from solving the averaged (over unbounded direction) of 
  !! the time independent non-homogeneous coupled electron-phonon 
  !! Boltzmann transport equation under the appropiate boundary conditions
  !! (i.e. fully diffusive boundaries, so that it reemits all incident 
  !! particles at reference equilibrium of the absorbing wall).
  !! See the subroutine compute_suppression for more information
  !!
  use precision, only: r64, i64
  use params, only: pi
  use misc, only: print_message, exit_with_message
  use symmetry_module, only: symmetry
  use phonon_module, only: phonon
  use electron_module, only: electron
  
  implicit none
  
  private
  public nanostructure
  
  !! This type contains informnation about nanostructuration
  type nanostructure
    
    character(len = 2), allocatable :: tag(:)
    !! Nanosystem identifier
    real(r64), allocatable :: limit(:)
    !! limiting length of the nanostructure
    real(r64), allocatable :: taxis(:,:)
    real(r64), allocatable :: tnorm(:)
    !! unbounded axis and its norm2
    real(r64), allocatable :: naxis(:,:)
    real(r64), allocatable :: nnorm(:)
    !! normal axis for planar structures and its norm2
    integer(i64) :: nsys 
    !! number of nanostructures
    real(r64), allocatable :: Sph(:,:,:) , Sel(:,:,:)
    !! the suppression factors for phonons and electrons
    logical, allocatable :: conv_ph(:), conv_el(:)
    !! convergence for each distribution for all system
    real(r64), allocatable :: vg_ph(:,:,:) 
    !! phonon group velocities (iq,ib, isystem)
    real(r64), allocatable :: vg_el(:,:,:)
    !! phonon group velocities (ik,ib, isystem)
    
   contains
   
   procedure, public :: initialize=>read_nanostructure, compute_suppression, &
                 compute_transport_vg, print_nanogeominfo, clean
  
  end type nanostructure
  
  contains

  subroutine clean(self)
    !! Subroutine to clean up the object
    !!
    !! self nansotructure object

    class(nanostructure), intent(inout) :: self
    
    if (allocated(self%tag))     deallocate(self%tag)
    if (allocated(self%limit))   deallocate(self%limit)
    if (allocated(self%taxis))   deallocate(self%taxis)
    if (allocated(self%tnorm))   deallocate(self%tnorm)
    if (allocated(self%naxis))   deallocate(self%naxis)
    if (allocated(self%nnorm))   deallocate(self%nnorm)
    if (allocated(self%Sph))     deallocate(self%Sph)
    if (allocated(self%Sel))     deallocate(self%Sel)
    if (allocated(self%conv_ph)) deallocate(self%conv_ph)
    if (allocated(self%conv_el)) deallocate(self%conv_el)
    if (allocated(self%vg_ph))   deallocate(self%vg_ph)
    if (allocated(self%vg_el))   deallocate(self%vg_el)

  end subroutine clean
  
  subroutine read_nanostructure(self)
    !! Subroutine to read nanostructure information
    !!
    !! self nanostructure object
    
    class(nanostructure), intent(inout) :: self
    
    character(len = 2), allocatable :: tag(:)
    !! Nanosystem identifier
    real(r64), allocatable :: limit(:)
    !! limiting length of the nanostructure
    real(r64), allocatable :: taxis(:,:)
    !! unbounded axis
    real(r64), allocatable :: naxis(:,:)
    !! normal axis for thin films
    integer(i64) :: nsys 
    !! Number of systems
    
    integer(i64) :: i
    !! loop counter
    integer(i64) :: nmax = 1000_i64
    !! max number of accepted structures 1000, should leave enough room (note than more than 25 makes input file difficult)
    
    namelist /nano/ tag, &
             nsys, taxis, naxis, limit 
    
    !Open input file
    open(1, file = 'input.nml', status = 'old')
    
    
    allocate(tag(nmax), limit(nmax), taxis(nmax,3), &
             naxis(nmax,3))
    
    !Defaults
    nsys  = -1
    tag   = 'xx'
    limit = -1.0_r64
    taxis = 0.0_r64
    naxis = 0.0_r64
    
    !Read nanostructures
    read(1, nml = nano)
    
    self%nsys = nsys
    
    if (self%nsys < 1 .or. nsys > nmax) then
        call exit_with_message('Bad input(s) in nanostructure. nsys must be [1,1000]')
    else
        allocate(self%tag(self%nsys), self%limit(self%nsys), self%taxis(self%nsys,3), &
             self%naxis(self%nsys,3), self%tnorm(self%nsys), self%nnorm(self%nsys))
        
        self%tag(:nsys) = tag(:nsys)
        self%limit(:nsys) = limit(:nsys)
        self%taxis(:nsys,:) = taxis(:nsys,:)
        self%naxis(:nsys,:) = naxis(:nsys,:)
        
        do i = 1, nsys
          self%tnorm(i) = norm2(self%taxis(i,:))
          self%nnorm(i) = norm2(self%naxis(i,:))
        end do 
    end if
    
    !Doing some checks
    do i = 1, nsys
        ! Check tags
        if ( (self%tag(i) .ne. 'nw') .and. &
             (self%tag(i) .ne. 'nr') .and. &
             (self%tag(i) .ne. 'tf')) then 
          call exit_with_message('Bad input(s) in nano: One tag is not valid')
        end if
        
        ! Check taxis not null
        if ( self%tnorm(i) < 1.0e-8_r64 ) then
          call exit_with_message(&
          'Bad input(s) in nano: Transport axis must have non-zero norm')
        end if
        
        ! Check normal axis for nr and thin films is not null
        if ( (self%tag(i) .ne. 'nw') .and. (self%nnorm(i) < 1.0e-8_r64 )) then
          call exit_with_message(&
          'Bad input(s) in nano: normal axis for nr and thin films must have non-zero norm.')
        end if
        
        ! Check normal and transport are perpedicular
        if ( (self%tag(i) .ne. 'nw') .and.  (abs(dot_product(self%naxis(i,:), & 
             self%taxis(i,:))) > 1.0e-8_r64 )) then
          call exit_with_message(& 
           'Bad input(s) in nano: normal axis must be perpedicular to transport aixs for nr and thin films.')
        end if
    end do
    
    deallocate(tag, limit, taxis, naxis)
    
    close(1)

  end subroutine read_nanostructure
  
  subroutine print_nanogeominfo(self)
    !! Prints nanostructure geometrical information to a file
    !!
    !! self nanostructure object
    !!
    class(nanostructure), intent(inout) :: self
    
    integer(i64) :: i
    
    if (this_image() == 1) then
        
        ! Open a file
        open(1, file = "nano_info", status = 'replace')
        
        do i = 1, self%nsys
            write(1, "(A , 7(1E20.10))") self%tag(i), self%limit(i), self%taxis(i,:), self%naxis(i,:)
        end do
        ! Close file
        close(1)
        
    end if
    
    sync all    
    
  end subroutine print_nanogeominfo
  
  subroutine compute_transport_vg(self, species_prefix, ph, el)
    !! Obtains the group velocity along the transport direction
    !! it rewrites the el and ph object velocities
    !!
    !! self nanostructure object
    !! species_prefix kind of particle
    !! ph phonon information
    !! el electron information
    !! 
    class(nanostructure) , intent(inout)  :: self
    character(len = 2), intent(in)        :: species_prefix
    type(phonon), intent(in), optional    :: ph
    type(electron), intent(in), optional  :: el
    
    !! Local
    integer(i64) :: i, j, k
    
    if (species_prefix == 'ph') then

      if (.not. present(ph)) then
        call exit_with_message("Error asked for ph vg in nanostructures but no ph object is provided")
      end if
      
      allocate(self%vg_ph(size(ph%vels,1),size(ph%vels,2),self%nsys))
      
      do i = 1, size(ph%vels,1)
        do j = 1, size(ph%vels,2)
          do k = 1, self%nsys
            self%vg_ph(i,j,k) = dot_product(ph%vels(i,j,:), &
                              self%taxis(k,:)/self%tnorm(k))
          end do
        end do
      end do
      
    else if (species_prefix == 'el') then

      if (.not. present(el)) then
        call exit_with_message("Error asked for el vg in nansotructures but no el object is provided")
      end if
      
      allocate(self%vg_el(size(el%vels,1),size(el%vels,2),self%nsys))
      
      do i = 1, size(el%vels,1)
        do j = 1, size(el%vels,2)
          do k = 1, self%nsys
            self%vg_el(i,j,k) = dot_product(el%vels(i,j,:), &
                              self%taxis(k,:)/self%tnorm(k))
          end do
        end do
      end do

    else
      call exit_with_message("Unknown particle in compute_suppression. Exiting.")
    end if
    
  end subroutine compute_transport_vg
  
  subroutine compute_suppression(self, species_prefix, sym, rta_rates_ibz, ph, el)
    !! Compute the suppression factor for carriers
    !! This function is based on the paper 10.1103/PhysRevB.85.195436
    !! getting the NW formula from such work. The formula for NR is 
    !! obtained from 10.1016/j.cpc.2022.108504 and the one for TF is 
    !! from 10.1016/j.cpc.2017.06.023.
    !! The model is based on the Chalmers solution of an 
    !! averaged non-homogeneous electron phonon Boltzmann transport
    !! Equation. It is averaged in the sense that RTA-deviations are 
    !! averaged over the cross section normal to the transport direction
    !! See 10.1016/j.ijheatmasstransfer.2024.125385 for more information regarding the methodology
    !! for the coupled electron-phonon BTE
    !!
    !! self Nanostructuration object
    !! species_prefix kind of particle
    !! sym Symmertry object
    !! rta_rates_ibz RTA scattering rates for the involved particle
    !! ph phonon information
    !! el electron information
    
    class(nanostructure) , intent(inout) :: self
    character(len = 2), intent(in) :: species_prefix
    type(symmetry), intent(in) :: sym
    real(r64), intent(in) :: rta_rates_ibz(:,:)
    type(phonon), intent(in)   :: ph
    type(electron), intent(in), optional :: el
    
    !! Local
    integer(i64) :: isys, i , l, il, s, ib
    real(r64) :: theta, phi, M, tau_ibz, vnr(2)
    real(r64) :: vg(3), u(3), rotnrinv(2,2), rot(3,3)
    
    !! Allocate 
    if (.not. allocated(self%Sph)) then
      allocate(self%Sph(size(ph%vels(:,1,1)),ph%numbands,self%nsys))
    end if
    if (present(el) .and. .not. allocated(self%Sel)) then
      allocate(self%Sel(el%nwv,el%numbands,self%nsys))
    end if
    
    do isys = 1 , self%nsys
      !! First procede with a rotation so that 
      !! for thin films and nanoribbons normal vector is aligned with z
      !! axis
      if (self%tag(isys) == 'nw') then
        rot(1,:) = (/1.0_r64, 0.0_r64, 0.0_r64/)
        rot(2,:) = (/0.0_r64, 1.0_r64, 0.0_r64/)
        rot(3,:) = (/0.0_r64, 0.0_r64, 1.0_r64/)
      else
        theta = acos(self%naxis(isys,3)/self%nnorm(isys))
        phi   = 0.0_r64
        if (abs(self%naxis(isys,3)/self%nnorm(isys)) < 1.0) then
          phi = atan2(self%naxis(isys,2)/self%nnorm(isys) , & 
                      self%naxis(isys,1)/self%nnorm(isys))
        end if
        rot(1,:) = (/cos(phi) * cos(theta), sin(phi) * cos(theta), -sin(theta)/)
        rot(2,:) = (/-sin(phi) , cos(phi), 0.0_r64/)
        rot(3,:) = (/cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)/)
      end if
      
      u = matmul(rot,self%taxis(isys,:)/self%tnorm(isys))
      !Some extra rotations for the nanoribbon case
      rotnrinv(1,:) = 1.0_r64 /(u(2)*u(2) + u(1)*u(1)) * (/ u(2) , -u(1) /)
      rotnrinv(2,:) = 1.0_r64 /(u(2)*u(2) + u(1)*u(1)) * (/ u(1) ,  u(2) /)
      
      if (species_prefix == 'ph') then
        do i = 1, ph%nwv_irred !an irreducible point
          do l = 1, ph%nequiv(i) !number of equivalent points of i
            il = ph%ibz2fbz_map(l, i, 2) ! (i, l) -> il
            s = ph%ibz2fbz_map(l, i, 1) ! mapping rotation
            !ik_ibz = ph%fbz2ibz_map(il)
            do ib = 1, ph%numbands
              !Get group velocity in the rotated coordinates
              vg = matmul(rot,ph%vels(il, ib, :))
              !Get RTA lifetime
              tau_ibz = 0.0_r64
              if(rta_rates_ibz(i, ib) /= 0.0_r64) then
                  tau_ibz = 1.0_r64/rta_rates_ibz(i, ib)
              end if
              !!! Here we compute the suppression ratio
              if (self%tag(isys) == 'nw') then
                M = norm2(vg - dot_product(vg,u) * u ) * tau_ibz
                if (M == 0.0_r64)  then
                  self%Sph(il,ib,isys) = 1.0_r64
                else
                  ! We call function as it requires computational integration
                  self%Sph(il,ib,isys) = supress_nw(self%limit(isys),M)
                end if
              else if (self%tag(isys) == 'nr') then
                vnr = matmul(rotnrinv,vg(:2))
                M = abs(vnr(1)) * tau_ibz
                if (M == 0.0_r64)  then
                  self%Sph(il,ib,isys) = 1.0_r64
                else
                  self%Sph(il,ib,isys) = 1.0_r64 + (M/self%limit(isys)) * &
                  ( exp(-self%limit(isys)/M) - 1.0_r64)
                end if
              else if (self%tag(isys) == 'tf') then
                M = abs(vg(3))/norm2(vg) * norm2(tau_ibz * vg) / self%limit(isys)
                if (M == 0.0_r64 .or. norm2(vg) < 1.0e-8_r64) then
                  self%Sph(il,ib,isys) = 1.0_r64
                else
                  self%Sph(il,ib,isys) = 1.0_r64 - M * (1.0_r64 - exp(-1.0_r64 / M))
                end if
              else 
                call exit_with_message("Unknown system in compute_suppression. Exiting. Exiting.")
              end if
            end do ! bands
          end do ! equivalences
        end do ! irreductible 
      else if (species_prefix == 'el') then
        !!! Do the same for electrons (be aware that it is done
        !!! on the restricted section )
        do i = 1, el%nwv !Iterate FBZ 
          il = el%fbz2ibz_map(i) ! Get the id of the point in the IBZ
          do ib = 1, el%numbands
            !Get group velocity in the rotated coordinates
            vg = matmul(rot,el%vels(i, ib, :))
            !Get RTA lifetime
            tau_ibz = 0.0_r64
            if(rta_rates_ibz(il, ib) /= 0.0_r64) then
              tau_ibz = 1.0_r64/rta_rates_ibz(il, ib)
            end if
            !!! Here compute the suppression factor
            if (self%tag(isys) == 'nw') then
              M = norm2(vg - dot_product(vg,u) * u ) * tau_ibz
              if (M == 0.0_r64)  then
                self%Sel(i,ib,isys) = 1.0_r64
              else
                ! We call function as it requires computational integration
                self%Sel(i,ib,isys) = supress_nw(self%limit(isys),M)
              end if
            else if (self%tag(isys) == 'nr') then
              vnr = matmul(rotnrinv,vg(:2))
              M = abs(vnr(1)) * tau_ibz
              if (M == 0.0_r64)  then
                self%Sel(i,ib,isys) = 1.0_r64
              else
                self%Sel(i,ib,isys) = 1.0_r64 + (M/self%limit(isys)) * &
                  ( exp(-self%limit(isys)/M) - 1.0_r64)
              end if
            else if (self%tag(isys) == 'tf') then
              M = abs(vg(3))/norm2(vg) * norm2(tau_ibz * vg) / self%limit(isys)
              if (M == 0.0_r64 .or. norm2(vg) < 1.0e-8_r64) then
                self%Sel(i,ib,isys) = 1.0_r64
              else
                self%Sel(i,ib,isys) = 1.0_r64 - M * (1.0_r64 - exp(-1.0_r64 / M))
              end if
            else 
              call exit_with_message("Unknown system in compute_suppression. Exiting. Exiting.")
            end if
          end do ! bands
        end do
      else 
        call exit_with_message("Unknown particle in compute_suppression. Exiting.")
      end if
    end do
    
    ! Saving information to file so that it can be plotted
    ! or used for other purposes
    if (this_image() == 1) then
        if (species_prefix == 'ph') then
            open(1, file = "Sph", status = 'replace')
            do isys = 1, self%nsys
                write(1,*) "# isys : ", isys
                do i = 1, size(self%Sph,1)
                    do ib = 1, ph%numbands
                        write(1,*) ph%ens(i, ib), ph%vels(i, ib, :), self%Sph(i,ib,isys)
                    end do
                end do
            end do
            close(1)
        else if (species_prefix == 'el') then
            open(1, file = "Sel", status = 'replace')
            do isys = 1, self%nsys
                write(1,*) "# isys : ", isys
                do i = 1, size(self%Sel,1)
                    do ib = 1, el%numbands
                        write(1,*) el%ens(i, ib), el%vels(i, ib, :), self%Sel(i,ib,isys)
                    end do
                end do
            end do
            close(1)
        else 
            call exit_with_message("Unknown particle in compute_suppression. Exiting.")
        end if
    end if
    sync all
    
    contains !
    
    !We implement as it is need for computation the trapezoidal rule
    !an analytical formula using Bessel and Struve modified functions exists
    !but has too much numerical noise. So numerical integration is the way to go.
    pure function supress_nw(R,mfp) result(sf)
      implicit none
      real(r64), intent(in) :: R, mfp
      real(r64) :: sf
      integer(i64) :: ir
      real(r64) :: dr
      real(r64) :: rs(5000)

      dr = R / size(rs)
      do ir = 1, size(rs)
        rs(ir) = (ir - 0.5_r64) * dr
      end do

      !First analytical part of the integral
      sf = -2.0_r64 * mfp / ( pi * R ) + 1.0_r64

      !The missing part \mathrm{\frac{M_{NW}}{\pi R^2}\int_{-R}^{R}e^\frac{-2\sqrt{R^2-x^2}}{M_{NW}}dx =
      !\frac{M_{NW}}{R}\left[I_1\left(\frac{2R}{M_{NW}}\right)-L_{-1}\left(\frac{2R}{M_{NW}}\right)\right ]}
      !with I_n and L_n being the n-th---1 and -1, respectively---order modified Bessel and Struve functions.
      !Numerical implementation based on such functions is full of numerical noise for 2*R/M_{NW} > 25
      !thus is better to evaluate it numerically
      do ir = 1, size(rs)
         sf = sf + 2.0_r64 * mfp/(pi * R**2 ) * exp( -2.0_r64*sqrt(R**2 - rs(ir)**2)/mfp) * dr
      end do
    end function supress_nw
    
  end subroutine compute_suppression
  
end module nano_module

